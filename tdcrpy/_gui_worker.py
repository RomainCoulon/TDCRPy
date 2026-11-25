# -*- coding: utf-8 -*-
"""Compute engine behind :mod:`tdcrpy.gui`, run as a separate process.

The window never imports ``tdcrpy`` itself. Every calculation is handed to a
fresh copy of this module over a pipe, which buys four things that matter
more here than the cost of starting a process:

* the window stays responsive during a run that can take minutes;
* Cancel is a kill, not a flag the Monte-Carlo loop has to check;
* a native crash in the numba kernel -- this package has produced SIGSEGV in
  long sweeps -- is reported as a failed run instead of taking the GUI with
  it;
* the configuration is applied before anything is imported, so no stale
  module global can survive into the calculation.

Protocol
--------
``--describe``
    Print one JSON object giving the version, the cocktail list, the nuclide
    list and the current configuration, then exit. The GUI populates its
    widgets from this rather than keeping its own copy of any of it.
``--run``
    Read one JSON job from stdin, emit JSON lines on stdout:
    ``{"event": "stage", ...}`` for progress, then exactly one of
    ``{"event": "result", ...}`` or ``{"event": "error", ...}``.

Anything the model prints, and anything ``warnings`` emits, goes to stderr,
which the GUI shows in its log pane. stdout carries only the protocol: the
real stdout is duplicated to a private file descriptor at start-up and
``sys.stdout`` is pointed at stderr, so a stray ``print`` inside the model
cannot corrupt a JSON line.
"""
import argparse
import json
import os
import sys
import traceback
import warnings

#: Efficiency channels in the order TDCRPy returns them, with the label the
#: GUI shows. Kept here, next to the unpacking code, so the two cannot drift.
CHANNELS = (
    ('S', 'Single'),
    ('D', 'Double (logical sum)'),
    ('T', 'Triple'),
    ('AB', 'Double AB'),
    ('BC', 'Double BC'),
    ('AC', 'Double AC'),
    ('D2', 'Double C/N'),
)


def _protocol_stream():
    """Hand back a writable stream on the real stdout, and divert prints.

    The model prints progress and, with ``display=True``, a full decay
    history. Those must not land between JSON lines.
    """
    fd = os.dup(1)
    os.dup2(2, 1)                     # anything printing to stdout -> stderr
    sys.stdout = sys.stderr
    return os.fdopen(fd, 'w', encoding='utf-8', newline='\n')


class Emitter:
    """Write protocol messages as one JSON object per line."""

    def __init__(self, stream):
        self._stream = stream

    def send(self, event, **fields):
        fields['event'] = event
        self._stream.write(json.dumps(fields, default=str) + '\n')
        self._stream.flush()


# --------------------------------------------------------------------------
# configuration
# --------------------------------------------------------------------------
#: Job keys that map straight onto a configuration key of the same name.
_DIRECT_KEYS = (
    'tau', 'extDT', 'measTime', 'opticalTransport',
    'micCorr', 'micCorrAuto', 'micCorrBiasThreshold', 'micCorrECThreshold',
    'diam_micelle', 'sigma_micelle',
    'sp_model', 'chou_param', 'nE_electron', 'nE_alpha',
)

#: Applied only after the cocktail, which resets them -- see _apply_config.
_DELAYED_KEYS = ('delayedFluorescence', 'tau_prompt', 'tau_delayed',
                 'eta_delayed')


def _apply_config(tl, cfg, emit):
    """Push the GUI's settings into ``config.toml``.

    Order is not cosmetic. ``modifyLScocktail`` re-derives the density, the
    effective Z and A and the whole atomic composition from the cocktail and
    the aqueous fraction, and it also resets ``tau_delayed`` and
    ``eta_delayed`` to that cocktail's fitted pair. So the cocktail goes
    first and the delayed-fluorescence parameters go last, or the window
    would show values the calculation did not use.
    """
    if not cfg:
        return

    if 'ls_cocktail' in cfg:
        tl.modifyLScocktail(cfg['ls_cocktail'],
                            float(cfg.get('fAq', 0.0)),
                            cfg.get('solvantType', 'Water'),
                            float(cfg.get('solvantConc', 0.0)))
    elif 'fAq' in cfg:
        tl.modifyfAq(float(cfg['fAq']))

    batch = {k: cfg[k] for k in _DIRECT_KEYS if k in cfg}
    if batch:
        tl.update_config_batch(batch)

    if 'effQuantum' in cfg:
        tl.modifyEffQ(cfg['effQuantum'])

    # Last, for the reason in the docstring.
    delayed = {k: cfg[k] for k in _DELAYED_KEYS if k in cfg}
    if delayed:
        tl.update_config_batch(delayed)

    emit.send('stage', text='configuration applied')


# --------------------------------------------------------------------------
# description
# --------------------------------------------------------------------------
def describe():
    """Report everything the GUI needs to build its widgets."""
    from importlib.metadata import version
    import tdcrpy.TDCR_model_lib as tl

    nuclides = sorted(
        n.split('.')[0] for n in tl.z_PenNuc.namelist()
        if n.endswith('.PenNuc.txt')
    )
    params = tl.readParameters()
    # Typical L per cocktail (2.20.37+): {name: [L, spread, measured]}.
    default_L = getattr(tl, 'defaultL', None)
    L_defaults = ({name: list(default_L(name)) for name in tl.COCKTAIL_DATA}
                  if default_L else {})
    return {
        'version': version('TDCRPy'),
        'cocktails': sorted(tl.COCKTAIL_DATA),
        'L_defaults': L_defaults,
        'nuclides': nuclides,
        # Absent before 2.20.34: the GUI then keeps a free-text field.
        'sp_models': list(getattr(tl, 'STOPPING_POWER_MODELS', ())),
        'Einterp_e': params[7],
        'config': {
            'ls_cocktail': tl.lsCocktail(),
            'fAq': params[9],
            'diam_micelle': params[8],
            'sigma_micelle': tl.sigma_micelle,
            'tau': params[10],
            'extDT': params[11],
            'measTime': params[12],
            'micCorr': bool(params[13]),
            'micCorrAuto': bool(tl.micCorrAuto),
            'micCorrBiasThreshold': tl.micCorrBiasThreshold,
            'micCorrECThreshold': tl.micCorrECThreshold,
            'effQuantum': tl.readEffQ0(),
            'opticalTransport': bool(params[15]),
            'solvantType': params[24],
            'solvantConc': params[25],
            'sp_model': params[26],
            'chou_param': params[27],
            'delayedFluorescence': bool(tl.delayedFluorescence),
            'tau_prompt': tl.tau_prompt,
            'tau_delayed': tl.tau_delayed,
            'eta_delayed': tl.eta_delayed,
            'nE_electron': params[0],
            'nE_alpha': params[1],
        },
    }


# --------------------------------------------------------------------------
# the calculations and data views
# --------------------------------------------------------------------------
def _channel_rows(out, analytical):
    """Rows for the seven channels of a 14-tuple ``(mean, u)`` result.

    The analytical model has no sampling noise and reports its uncertainties
    as 0.0; shown as such they would read as "exactly known", so they are
    dropped instead.
    """
    return [(f'Efficiency {label}', out[2 * i],
             None if analytical else out[2 * i + 1])
            for i, (_key, label) in enumerate(CHANNELS)]


def _analytical_channels(eff_S, eff_D, eff_T, symm):
    """Expand ``effA()``'s three efficiencies into the seven channels.

    ``effA()`` returns S, D and T only. For a symmetric counter the pair
    efficiencies follow from ``D = AB + BC + AC - 2T`` with AB = BC = AC, and
    the C/N doubles equal D -- the same expansion ``TDCRPy(Smodel=False)``
    makes, done here from effA's own values so they stay consistent with the
    fitted L. An asymmetric counter has no such closed form, so its pair
    channels are left empty rather than filled with the symmetric guess.
    """
    pair = (eff_D + 2 * eff_T) / 3 if symm else None
    return (eff_S, None, eff_D, None, eff_T, None,
            pair, None, pair, None, pair, None,
            eff_D if symm else None, None)


def _run_efficiency(mod, job, emit, unc_data=False):
    """``TDCRPy()``: efficiencies at a light yield the user supplies."""
    emit.send('stage', text='simulating decays')
    smodel = bool(job.get('Smodel', True))
    out = mod.TDCRPy(
        job['L'], job['rad'], job['pmf'], int(job['N']),
        float(job['kB']), float(job['V']),
        mode='eff',
        display=bool(job.get('display')),
        barp=True,
        Smodel=smodel,
        uncData=unc_data,
    )
    return [('L (fixed)', job['L'], None)] + _channel_rows(out, not smodel)


#: Columns of a result propagated over realisations of the nuclear data.
DATA_UNC_COLUMNS = (
    ('quantity', 'Quantity'), ('value', 'Value'), ('u', 'u (combined)'),
    ('u_rel', 'u / value'), ('u_MC', 'u (Monte Carlo)'),
    ('u_data', 'u (nuclear data)'),
)


def combine_realisations(runs):
    """Combine M runs, each on its own realisation of the nuclear data.

    *runs* is a list of row lists ``(quantity, value, u)``, one per run, all
    in the same order. For each quantity, with run-to-run variance ``s^2``
    and mean Monte-Carlo variance ``ū^2`` of a single run:

    * value     -- the mean over the runs;
    * u(MC)     -- ``ū / sqrt(M)``, the counting noise left in that mean;
    * u(data)   -- ``sqrt(s^2 - ū^2)``: what the data add to a single run's
      scatter (JCGM 101). An input uncertainty, so not divided by sqrt(M).

    ``s^2 - ū^2`` is a difference of two noisy numbers. Unless it exceeds
    twice the standard error of ``s^2`` it is not resolved, and the 95 %
    upper bound is reported instead, as the string ``'< x'``: claiming a
    value there would be reporting Monte-Carlo noise as a data effect.

    Returns
    -------
    list of list
        ``[quantity, value, u, u_rel, u_MC, u_data]`` per quantity.
    """
    import numpy as np
    from scipy.stats import chi2
    M = len(runs)
    out = []
    for k, first in enumerate(runs[0]):
        quantity = first[0]
        values = [r[k][1] for r in runs]
        if isinstance(values[0], (list, tuple)):          # three L values
            mean = [float(np.mean(c)) for c in zip(*values)]
            out.append([quantity, ', '.join(f'{v:.6g}' for v in mean),
                        None, None, None, None])
            continue
        if not all(isinstance(v, (int, float)) for v in values):
            out.append([quantity, values[0], None, None, None, None])
            continue
        v = np.asarray(values, dtype=float)
        # An input repeated unchanged (the fixed L) stays exactly itself.
        mean = float(v[0]) if np.all(v == v[0]) else float(v.mean())
        us = [r[k][2] for r in runs]
        if any(u is None for u in us):
            out.append([quantity, mean, None, None, None, None])
            continue
        ubar2 = float(np.mean(np.square(us)))
        s2 = float(v.var(ddof=1))
        u_mc = (ubar2 / M) ** 0.5
        excess = s2 - ubar2
        if excess > 2.0 * s2 * (2.0 / (M - 1)) ** 0.5:
            u_data = excess ** 0.5
            u = (u_mc ** 2 + excess) ** 0.5
            shown = u_data
        else:
            upper = s2 * (M - 1) / chi2.ppf(0.05, M - 1) - ubar2
            u = u_mc
            shown = f'< {max(upper, 0.0) ** 0.5:.2g}' if upper > 0 else '0'
        rel = f'{100.0 * u / mean:.2f} %' if mean else ''
        out.append([quantity, mean, u, rel, u_mc, shown])
    return out


def _run_data_realisations(mod, job, emit, once):
    """Repeat *once* over M realisations of the nuclear data and combine."""
    M = int(job.get('uncReps', 10))
    runs = []
    for i in range(M):
        emit.send('stage', text=f'nuclear-data realisation {i + 1} of {M}')
        runs.append(once())
        emit.send('progress', pct=100.0 * (i + 1) / M)
    emit.send('stage', text=f'combined {M} realisations of {job["N"]} decays '
                            f'each; u (nuclear data) shown as "< x" is not '
                            f'resolved above the Monte-Carlo noise')
    return combine_realisations(runs), DATA_UNC_COLUMNS


def _run_fit(mod, job, emit, analytical, unc_data=False):
    """``eff()`` / ``effA()``: fit the light yield to a measured TDCR."""
    emit.send('stage', text='fitting L to the measured TDCR')
    kwargs = dict(L=float(job.get('L_guess', 1.0)),
                  disp=bool(job.get('display')))
    if job.get('Lbounds'):
        kwargs['Lbounds'] = tuple(job['Lbounds'])
    if unc_data:
        kwargs['uncData'] = True

    if analytical:
        # effA() returns (L0, L, eff_S, eff_D, eff_T) -- five values, not
        # eff()'s seventeen -- and scripts index into it, so it is expanded
        # here rather than changed there.
        L0, L_opt, eff_S, eff_D, eff_T = mod.effA(
            job['TD'], job['rad'], job['pmf'],
            float(job['kB']), float(job['V']), **kwargs)
        symm = not isinstance(job['TD'], (list, tuple))
        out = _analytical_channels(eff_S, eff_D, eff_T, symm)
        u_L = None
        if not symm:
            emit.send('stage', text='analytical model, asymmetric counter: '
                                    'pair efficiencies not available')
    else:
        res = mod.eff(job['TD'], job['rad'], job['pmf'],
                      float(job['kB']), float(job['V']),
                      N=int(job['N']), **kwargs)
        L0, L_opt, out = res[0], res[1], res[2:16]
        # u_L was appended in 2.20.23, so trust the length rather than a
        # version test: older builds simply do not report it.
        u_L = res[16] if len(res) > 16 else None

    rows = [('L (fitted)', L0, u_L),
            ('L_A, L_B, L_C', [float(v) for v in L_opt], None)]
    return rows + _channel_rows(out, analytical)


#: Columns of the efficiency-curve table, as (key, heading).
CURVE_COLUMNS = (
    ('L', 'L / keV⁻¹'), ('TDCR', 'TDCR'),
    ('eff_D', 'ε_D'), ('u_eff_D', 'u(ε_D)'),
    ('eff_T', 'ε_T'), ('u_eff_T', 'u(ε_T)'),
    ('eff_S', 'ε_S'), ('u_eff_S', 'u(ε_S)'),
)


def _run_curve(mod, job, emit):
    """``TDCRPy()`` over a grid of L: the efficiency curve ε_D vs TDCR.

    The stochastic model simulates the decays once and replays the recorded
    histories at every L, exactly as :func:`tdcrpy.TDCRPy.eff` does while
    fitting, so a curve costs about one simulation rather than one per point
    -- and every point sees the same decays, so the curve is smooth rather
    than carrying independent Monte-Carlo scatter from point to point.
    """
    grid = [float(v) for v in job['L_grid']]
    smodel = bool(job.get('Smodel', True))
    args = (job['rad'], job['pmf'], int(job['N']),
            float(job['kB']), float(job['V']))
    if smodel:
        emit.send('stage', text='simulating decays')
        # No uncData: one realisation of the data would shift the whole
        # curve and look like a result; its spread needs repeated curves.
        mod.TDCRPy(grid[0], *args, record=True, barp=True,
                   display=bool(job.get('display')))
        emit.send('stage', text=f'replaying them at {len(grid)} values of L')
    rows = []
    for i, L in enumerate(grid):
        if smodel:
            out = mod.TDCRPy(L, *args, readRecHist=True)
        else:
            out = mod.TDCRPy(L, *args, Smodel=False)
        eff_S, u_S, eff_D, u_D, eff_T, u_T = out[:6]
        if not smodel:
            u_S = u_D = u_T = None
        tdcr = eff_T / eff_D if eff_D > 1e-12 else None
        rows.append([L, tdcr, eff_D, u_D, eff_T, u_T, eff_S, u_S])
        emit.send('progress', pct=100.0 * (i + 1) / len(grid))
    return rows


#: Levels of the BetaShape archive Co60-style names can take, besides an
#: integer transition number: a pure beta emitter has only 'tot'.
_BETA_LITERAL_LEVELS = ('tot', 'tot_myEstep')


def _beta_level(text):
    """Turn the GUI's free-text level field into what readBetaShape() wants:
    an int transition number, or one of the literal level names."""
    text = (text or 'tot').strip()
    if text in _BETA_LITERAL_LEVELS:
        return text
    try:
        return int(text)
    except ValueError:
        raise ValueError(f"'{text}' is not a transition number nor one of "
                         f"{_BETA_LITERAL_LEVELS}")


def _available_beta_levels(tl, rad, mode):
    """List the levels the archive actually has for *rad* and *mode*, for
    the error message when the requested one is missing."""
    prefix = f"{rad.replace('-', '')}/{mode}_{rad.replace('-', '')}_"
    out = []
    for name in tl.z_betashape.namelist():
        if name.startswith(prefix) and name.endswith('.bs'):
            tail = name[len(prefix):-3]
            if tail.startswith('trans'):
                tail = tail[len('trans'):]
            out.append(tail)
    return sorted(set(out))


#: Columns of the beta-spectrum table, as (key, heading).
BETA_SPECTRUM_COLUMNS = (('E_keV', 'E / keV'), ('p', 'p(E) / keV⁻¹'))


def _run_beta_spectrum(mod, job, emit):
    """``readBetaShape()``: the emitted-energy distribution of one beta
    transition, read straight from the shipped BetaShape archive -- no
    Monte Carlo, so this returns at once.
    """
    import tdcrpy.TDCR_model_lib as tl
    rad = job['rad']
    mode = job.get('mode', 'beta-')
    level = _beta_level(job.get('level'))
    emit.send('stage', text=f'reading the {mode} spectrum of {rad}')
    try:
        e, p = tl.readBetaShape(rad, mode, level)
    except KeyError:
        avail = _available_beta_levels(tl, rad, mode)
        raise ValueError(
            f"{rad} has no {mode} spectrum at level {level!r}. "
            f"Levels available for {rad} {mode}: {avail or '(none)'}.")
    rows = [[ei, pi] for ei, pi in zip(e, p)]
    plot = {'x': 'E_keV', 'y': ['p'], 'yerr': [None],
            'xlabel': 'E / keV', 'ylabel': 'p(E) / keV⁻¹'}
    return (rows, BETA_SPECTRUM_COLUMNS, plot,
            _decay_info(tl, rad, emit, only=(mode, str(level))))


def _run_stopping_power(mod, job, emit):
    """Stopping power over an energy grid, for the configured scintillator.

    Electrons use ``stoppingpower()`` -- with ``compare_models``, every model
    the library knows, at the same density/Z/A so the comparison is apples to
    apples. Alphas use ``stoppingpowerA()`` (ASTAR), converted from keV cm⁻¹
    to the same MeV cm⁻¹ so the two plots read on one scale.
    """
    import tdcrpy.TDCR_model_lib as tl
    grid = [float(v) for v in job['E_grid']]            # keV
    emit.send('stage', text='evaluating the stopping power')
    if job.get('particle', 'electron') == 'alpha':
        rows = [[E, tl.stoppingpowerA(E) * 1e-3] for E in grid]
        labels = ['alpha (ASTAR)']
    else:
        compare = bool(job.get('compare_models'))
        models = list(tl.STOPPING_POWER_MODELS) if compare else [None]
        rows = [[E] + [tl.stoppingpower(E * 1e3, spmodel=m) for m in models]
                for E in grid]
        labels = models if compare else [tl._stopping_globals()[3]]
    columns = (('E_keV', 'E / keV'),) + tuple((lab, lab) for lab in labels)
    plot = {'x': 'E_keV', 'y': list(labels), 'yerr': [None] * len(labels),
            'xlabel': 'E / keV', 'ylabel': 'dE/dx / MeV cm⁻¹',
            'xlog': True, 'ylog': True}
    return rows, columns, plot


def _run_quenching(mod, job, emit):
    """Birks-quenched energy Q(E) and the quenching factor Q/E.

    Computed with ``Em_e()`` / ``Em_a()``, the very functions the model
    calls -- exact integral below the interpolation threshold, the shipped
    tables above it -- so the curve is what a run actually uses. One series
    per Birks constant, to compare kB values on one plot.
    """
    import tdcrpy.TDCR_model_lib as tl
    grid = [float(v) for v in job['E_grid']]            # keV
    kBs = [float(k) for k in job['kB_list']]            # cm/keV
    alpha = job.get('particle', 'electron') == 'alpha'
    params = tl.readParameters()
    nE_e, nE_a, E_int_a, E_int_e = params[0], params[1], params[6], params[7]

    def quenched(E, kB):
        if alpha:
            return tl.Em_a(E, kB, nE_a, Et=E_int_a)
        # Em_e works in eV and cm/MeV, as TDCRPy() calls it.
        return tl.Em_e(E * 1e3, E * 1e3, kB * 1e3, nE_e,
                       Et=E_int_e * 1e3) * 1e-3

    emit.send('stage', text='integrating the Birks quenching function')
    rows = []
    for i, E in enumerate(grid):
        row = [E]
        for kB in kBs:
            q = quenched(E, kB)
            row += [q, q / E]
        rows.append(row)
        emit.send('progress', pct=100.0 * (i + 1) / len(grid))
    columns = [('E_keV', 'E / keV')]
    for kB in kBs:
        columns += [(f'Q_{kB:g}', f'Q / keV  (kB={kB:g})'),
                    (f'QE_{kB:g}', f'Q/E  (kB={kB:g})')]
    keys = [f'QE_{kB:g}' for kB in kBs]
    plot = {'x': 'E_keV', 'y': keys, 'yerr': [None] * len(keys),
            'labels': [f'kB = {kB:g} cm/keV' for kB in kBs],
            'xlabel': 'E / keV', 'ylabel': 'Q(E) / E',
            'xlog': True}
    return rows, tuple(columns), plot


#: How PenNuc's record codes read in the table.
_BRANCH_LABEL = {'beta': 'β⁻', 'beta+': 'β⁺', 'alpha': 'α',
                 'Atom_K': 'EC (K)', 'Atom_L': 'EC (L)', 'Atom_L1': 'EC (L1)',
                 'Atom_L2': 'EC (L2)', 'Atom_L3': 'EC (L3)',
                 'Atom_M': 'EC (M)', 'Atom_N': 'EC (N)', 'Atom_O': 'EC (O)'}
_TRANSITION_LABEL = {'GA': 'γ', 'EK': 'CE K', 'EL': 'CE L', 'EL1': 'CE L1',
                     'EL2': 'CE L2', 'EL3': 'CE L3', 'EM': 'CE M',
                     'EN': 'CE N'}

#: Columns of the decay-data table, as (key, heading).
DECAY_DATA_COLUMNS = (
    ('daughter', 'Daughter'), ('record', 'Record'), ('type', 'Type'),
    ('E_keV', 'E / keV'), ('P', 'P per decay'), ('u_P', 'u(P)'),
    ('from_keV', 'From level / keV'), ('to_keV', 'To level / keV'),
    ('line', 'Discrete line'),
)


def decay_data_rows(out):
    """Flatten ``readPenNuc2()`` into one row per decay branch and per
    de-excitation transition.

    Every probability in PenNuc is per decay of the parent: the branch
    probabilities of a daughter sum to its own probability (K-40, Bi-212),
    and the transition probabilities are absolute emission intensities
    (the 59.54 keV γ of Am-241 is 0.3592). The 'line' column names the
    discrete emissions -- α, γ and conversion electrons -- which are the
    ones a line plot can show; a β entry is the end point of a continuum.
    """
    (daughters, _p_daug, _Q, b_type, b_energy, b_prob, b_level, _pb,
     t_type, t_energy, t_prob, t_level, t_level_end, level_energy,
     _pt, _hl, t_unc) = out[:17]
    rows = []
    for d, name in enumerate(daughters):
        # Level index -> energy, from the de-excitation blocks; 0 is ground.
        energy_of = {0: 0.0}
        for lev, en in zip(t_level[d], level_energy[d]):
            if lev and en:
                energy_of[lev[0]] = en[0]
        for types, energies, probs, levels in zip(b_type[d], b_energy[d],
                                                 b_prob[d], b_level[d]):
            for t, e, p, lev in zip(types, energies, probs, levels):
                label = _BRANCH_LABEL.get(t, t)
                rows.append([name, 'decay', label, e, p, None, None,
                             energy_of.get(lev), 'α' if t == 'alpha' else ''])
        for lev, types, energies, probs, uncs, ends in zip(
                t_level[d], t_type[d], t_energy[d], t_prob[d], t_unc[d],
                t_level_end[d]):
            if not lev:
                continue
            for t, e, p, u, end in zip(types, energies, probs, uncs, ends):
                label = _TRANSITION_LABEL.get(t, t)
                line = 'γ' if t == 'GA' else 'conversion e⁻'
                rows.append([name, 'de-excitation', label, e, p, u,
                             energy_of.get(lev[0]), energy_of.get(end),
                             line])
    return rows


#: How BetaShape's decay-mode prefixes read.
_BETA_MODE_LABEL = {'beta-': 'β⁻', 'beta+': 'β⁺', 'capt': 'EC'}


def _beta_level_label(level):
    fine = level.endswith('_myEstep')
    base = level[:-len('_myEstep')] if fine else level
    text = 'total spectrum' if base == 'tot' else f'transition {base}'
    return text + (', fine energy step' if fine else '')


def provenance(info, only=None):
    """Turn :func:`decayDataInfo` into sections of ``[label, text, url]``.

    *only* ``(mode, level)`` keeps just that BetaShape file, for the
    spectrum view.
    """
    nuc = info['nuclear']
    items = [['Source', nuc['source'], None],
             ['DDEP evaluation date', nuc['evaluation_date'] or
              'not stated in the file', None]]
    for name, p, q, uq in nuc['daughters']:
        text = f'{name}: P = {p:g}' if p is not None else name
        if q is not None:
            text += f', Q = {q:g}' + (f' ± {uq:g}' if uq else '') + ' keV'
        items.append(['Daughter', text, None])
    for note in nuc['notes']:
        items.append(['Note', note, None])
    if nuc['tables_url']:
        items.append(['DDEP tables (LNHB)', nuc['tables_url'],
                      nuc['tables_url']])
        items.append(['DDEP evaluation comments', nuc['comments_url'],
                      nuc['comments_url']])
    sections = [{'title': 'Nuclear decay data — PenNuc', 'items': items}]

    atomic = info['atomic']
    if atomic:
        items = [['Decay', d, None] for d in atomic['decays']]
        if atomic['evaluator']:
            kind = {'FUL': 'full evaluation', 'UPD': 'update'}.get(
                atomic['type'], atomic['type'])
            items.append(['Evaluator', f'{atomic["evaluator"]} ({kind}, '
                                       f'literature cut-off '
                                       f'{atomic["cutoff"]})', None])
        sections.append({'title': 'X-rays and Auger electrons — ENSDF '
                                  'format', 'items': items})

    items = []
    for b in info['beta']:
        if only and (b['mode'], b['level']) != only:
            continue
        head = (f'{_BETA_MODE_LABEL.get(b["mode"], b["mode"])}, '
                f'{_beta_level_label(b["level"])}')
        items.append([head, b.get('program', '?'), None])
        for key, label in (('transition', 'calculated as'),
                           ('end_point', 'end point'),
                           ('shape_factor', 'experimental shape factor'),
                           ('cite', 'cite')):
            if b.get(key):
                items.append(['', f'{label}: {b[key]}', None])
    if info['beta_alternatives']:
        items.append(['Also archived', 'earlier spectrum set(s) '
                      + ', '.join(info['beta_alternatives'])
                      + ' — not used by default', None])
    if items:
        sections.append({'title': 'β spectra — BetaShape', 'items': items})
    return sections


def _decay_info(tl, rad, emit, only=None):
    """Provenance sections, or ``None`` on a build without decayDataInfo."""
    if not hasattr(tl, 'decayDataInfo'):
        return None
    try:
        return provenance(tl.decayDataInfo(rad), only)
    except Exception as exc:                  # never lose the data over it
        emit.send('stage', text=f'provenance unavailable: {exc}')
        return None


def _run_decay_data(mod, job, emit):
    """``readPenNuc2()``: the decay scheme as the model reads it."""
    import tdcrpy.TDCR_model_lib as tl
    rad = job['rad']
    emit.send('stage', text=f'reading the PenNuc decay data of {rad}')
    try:
        out = tl.readPenNuc2(rad)
    except KeyError:
        raise ValueError(f'{rad} is not in the PenNuc decay-data archive.')
    rows = decay_data_rows(out)
    emit.send('stage', text='β energies are end points; EC energies are the '
                            'transition energy. X-rays and Auger electrons '
                            'come from atomic relaxation and are not listed.')
    plot = {'kind': 'stem', 'x': 'E_keV', 'y': ['P'], 'group': 'line',
            'xlabel': 'E / keV', 'ylabel': 'emission probability per decay',
            'ylog': True}
    return rows, DECAY_DATA_COLUMNS, plot, _decay_info(tl, rad, emit)


#: Tasks whose result is a curve to plot, rather than a single set of
#: quantities: the result carries 'columns' and 'plot' alongside 'rows'.
PLOT_TASKS = {'curve', 'beta_spectrum', 'stopping_power', 'quenching',
              'decay_data'}

#: Calculations that can be repeated over realisations of the nuclear data.
DATA_UNC_TASKS = {
    'efficiency': lambda m, j, e: _run_efficiency(m, j, e, unc_data=True),
    'fit': lambda m, j, e: _run_fit(m, j, e, analytical=False,
                                    unc_data=True),
}

_RUNNERS = {
    'efficiency': lambda m, j, e: _run_efficiency(m, j, e),
    'fit': lambda m, j, e: _run_fit(m, j, e, analytical=False),
    'fit_analytical': lambda m, j, e: _run_fit(m, j, e, analytical=True),
    'curve': lambda m, j, e: _run_curve(m, j, e),
    'beta_spectrum': lambda m, j, e: _run_beta_spectrum(m, j, e),
    'stopping_power': lambda m, j, e: _run_stopping_power(m, j, e),
    'quenching': lambda m, j, e: _run_quenching(m, j, e),
    'decay_data': lambda m, j, e: _run_decay_data(m, j, e),
}


def run(job, emit):
    """Apply the configuration, then perform the requested calculation."""
    # Warnings are the model's channel for the micelle advisories and the
    # unfitted-cocktail notice, so they belong in the log rather than being
    # swallowed.
    warnings.simplefilter('always')

    emit.send('stage', text='loading TDCRPy')
    import tdcrpy.TDCRPy as mod
    import tdcrpy.TDCR_model_lib as tl

    _apply_config(tl, job.get('config'), emit)

    try:
        runner = _RUNNERS[job['task']]
    except KeyError:
        raise ValueError(f'unknown task {job.get("task")!r}')

    if job.get('uncData') and job['task'] in DATA_UNC_TASKS \
            and job.get('Smodel', True):
        once = DATA_UNC_TASKS[job['task']]
        rows, columns = _run_data_realisations(
            mod, job, emit, lambda: plain_rows(once(mod, job, emit)))
        emit.send('result', rows=plain_rows(rows), task=job['task'],
                  columns=[list(c) for c in columns], plot=None)
    elif job['task'] == 'curve':
        rows = runner(mod, job, emit)
        emit.send('result', rows=plain_rows(rows), task='curve',
                  columns=[list(c) for c in CURVE_COLUMNS],
                  plot={'x': 'TDCR', 'y': ['eff_D'], 'yerr': ['u_eff_D'],
                        'xlabel': 'TDCR  (T/D)', 'ylabel': 'ε_D',
                        'annotate': 'L'})
    elif job['task'] in PLOT_TASKS:
        rows, columns, plot, *extra = runner(mod, job, emit)
        emit.send('result', rows=plain_rows(rows), task=job['task'],
                  columns=[list(c) for c in columns], plot=plot,
                  info=extra[0] if extra else None)
    else:
        rows = runner(mod, job, emit)
        emit.send('result', rows=plain_rows(rows), task=job['task'])


def _plain(v):
    """A numpy scalar or 0-d array as a Python float; anything else as is.

    ``json.dumps(default=str)`` would otherwise send ``Em_a()``'s 0-d arrays
    as strings, which the window can neither format nor plot.
    """
    if isinstance(v, (str, bool)) or v is None:
        return v
    try:
        return float(v)
    except (TypeError, ValueError):
        return v


def plain_rows(rows):
    """Apply :func:`_plain` to every cell of every row."""
    return [[_plain(v) for v in row] for row in rows]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--describe', action='store_true',
                        help='print capabilities and current config as JSON')
    parser.add_argument('--run', action='store_true',
                        help='read one JSON job from stdin and run it')
    args = parser.parse_args(argv)

    out = _protocol_stream()
    emit = Emitter(out)

    try:
        if args.describe:
            emit.send('describe', **describe())
            return 0
        if args.run:
            job = json.loads(sys.stdin.read())
            run(job, emit)
            return 0
        parser.error('one of --describe or --run is required')
    except Exception as exc:                     # reported, never printed raw
        emit.send('error', message=f'{type(exc).__name__}: {exc}',
                  traceback=traceback.format_exc())
        return 1


if __name__ == '__main__':
    sys.exit(main())
