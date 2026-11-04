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
# the three calculations
# --------------------------------------------------------------------------
def _run_efficiency(mod, job, emit):
    """``TDCRPy()``: efficiencies at a light yield the user supplies."""
    emit.send('stage', text='simulating decays')
    out = mod.TDCRPy(
        job['L'], job['rad'], job['pmf'], int(job['N']),
        float(job['kB']), float(job['V']),
        mode='eff',
        display=bool(job.get('display')),
        barp=True,
        Smodel=bool(job.get('Smodel', True)),
        uncData=bool(job.get('uncData')),
    )
    rows = [('L (fixed)', job['L'], None)]
    for i, (_key, label) in enumerate(CHANNELS):
        rows.append((f'Efficiency {label}', out[2 * i], out[2 * i + 1]))
    return rows


def _run_fit(mod, job, emit, analytical):
    """``eff()`` / ``effA()``: fit the light yield to a measured TDCR."""
    emit.send('stage', text='fitting L to the measured TDCR')
    kwargs = dict(L=float(job.get('L_guess', 1.0)),
                  disp=bool(job.get('display')))
    if job.get('Lbounds'):
        kwargs['Lbounds'] = tuple(job['Lbounds'])

    if analytical:
        out = mod.effA(job['TD'], job['rad'], job['pmf'],
                       float(job['kB']), float(job['V']), **kwargs)
        u_L = None
    else:
        out = mod.eff(job['TD'], job['rad'], job['pmf'],
                      float(job['kB']), float(job['V']),
                      N=int(job['N']), **kwargs)
        # u_L was appended in 2.20.23, so trust the length rather than a
        # version test: older builds simply do not report it.
        u_L = out[16] if len(out) > 16 else None

    rows = [('L (fitted)', out[0], u_L),
            ('L_A, L_B, L_C', ', '.join(f'{v:.6g}' for v in out[1]), None)]
    for i, (_key, label) in enumerate(CHANNELS):
        rows.append((f'Efficiency {label}', out[2 + 2 * i], out[3 + 2 * i]))
    return rows


_RUNNERS = {
    'efficiency': lambda m, j, e: _run_efficiency(m, j, e),
    'fit': lambda m, j, e: _run_fit(m, j, e, analytical=False),
    'fit_analytical': lambda m, j, e: _run_fit(m, j, e, analytical=True),
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
