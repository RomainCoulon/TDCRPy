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
    return {
        'version': version('TDCRPy'),
        'cocktails': sorted(tl.COCKTAIL_DATA),
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
# the six calculations
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


def _run_efficiency(mod, job, emit):
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
        uncData=bool(job.get('uncData')),
    )
    return [('L (fixed)', job['L'], None)] + _channel_rows(out, not smodel)


def _run_fit(mod, job, emit, analytical):
    """``eff()`` / ``effA()``: fit the light yield to a measured TDCR."""
    emit.send('stage', text='fitting L to the measured TDCR')
    kwargs = dict(L=float(job.get('L_guess', 1.0)),
                  disp=bool(job.get('display')))
    if job.get('Lbounds'):
        kwargs['Lbounds'] = tuple(job['Lbounds'])

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
            ('L_A, L_B, L_C', ', '.join(f'{v:.6g}' for v in L_opt), None)]
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
        mod.TDCRPy(grid[0], *args, record=True, barp=True,
                   display=bool(job.get('display')),
                   uncData=bool(job.get('uncData')))
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
    return rows, BETA_SPECTRUM_COLUMNS, plot


def _run_stopping_power(mod, job, emit):
    """``stoppingpower()``: electron dE/dx over an energy grid, for the
    scintillator currently configured -- or, with ``compare_models``, for
    every model the library knows, at that same density/Z/A so the
    comparison is apples to apples.
    """
    import tdcrpy.TDCR_model_lib as tl
    grid = [float(v) for v in job['E_grid']]            # keV
    compare = bool(job.get('compare_models'))
    models = list(tl.STOPPING_POWER_MODELS) if compare else [None]
    emit.send('stage', text='evaluating the stopping power')
    rows = [[E] + [tl.stoppingpower(E * 1e3, spmodel=m) for m in models]
            for E in grid]
    labels = models if compare else [tl._stopping_globals()[3]]
    columns = (('E_keV', 'E / keV'),) + tuple((lab, lab) for lab in labels)
    plot = {'x': 'E_keV', 'y': list(labels), 'yerr': [None] * len(labels),
            'xlabel': 'E / keV', 'ylabel': 'dE/dx / MeV cm⁻¹',
            'xlog': True, 'ylog': True}
    return rows, columns, plot


#: Tasks whose result is a curve to plot, rather than a single set of
#: quantities: the result carries 'columns' and 'plot' alongside 'rows'.
PLOT_TASKS = {'curve', 'beta_spectrum', 'stopping_power'}

_RUNNERS = {
    'efficiency': lambda m, j, e: _run_efficiency(m, j, e),
    'fit': lambda m, j, e: _run_fit(m, j, e, analytical=False),
    'fit_analytical': lambda m, j, e: _run_fit(m, j, e, analytical=True),
    'curve': lambda m, j, e: _run_curve(m, j, e),
    'beta_spectrum': lambda m, j, e: _run_beta_spectrum(m, j, e),
    'stopping_power': lambda m, j, e: _run_stopping_power(m, j, e),
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

    if job['task'] == 'curve':
        rows = runner(mod, job, emit)
        emit.send('result', rows=rows, task='curve',
                  columns=[list(c) for c in CURVE_COLUMNS],
                  plot={'x': 'TDCR', 'y': ['eff_D'], 'yerr': ['u_eff_D'],
                        'xlabel': 'TDCR  (T/D)', 'ylabel': 'ε_D',
                        'annotate': 'L'})
    elif job['task'] in PLOT_TASKS:
        rows, columns, plot = runner(mod, job, emit)
        emit.send('result', rows=rows, task=job['task'],
                  columns=[list(c) for c in columns], plot=plot)
    else:
        rows = runner(mod, job, emit)
        emit.send('result', rows=rows, task=job['task'])


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
