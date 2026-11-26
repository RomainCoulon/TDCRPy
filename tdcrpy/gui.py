# -*- coding: utf-8 -*-
"""TDCRPy — a window for driving the TDCR model without writing a script.

Start it with ``tdcrpy-gui``, or ``python -m tdcrpy.gui``.

Calculations
    efficiencies at a light yield you supply (:func:`tdcrpy.TDCRPy.TDCRPy`);
    the light yield fitted to a measured TDCR ratio by the stochastic
    (:func:`~tdcrpy.TDCRPy.eff`) or the analytical
    (:func:`~tdcrpy.TDCRPy.effA`) model; and the efficiency curve, ε_D
    against TDCR, over a range of light yields.
Data & plots
    the decay data the model reads (PenNuc); the emitted beta spectrum of a
    transition (BetaShape); the stopping power of electrons (every model) or
    alpha particles (ASTAR); and the Birks ionisation-quenching function.

Every result can be exported as CSV, and every curve-shaped one is plotted
in the window and can be saved as PNG, PDF or SVG.

With the verbose switch on, a calculation also returns its first decays as
structured records (``TDCRPy(trace=)``), drawn in the Decay view tab as the
chain the model follows: decay, emission, deposition, quenching, signal.

Design notes
------------
This module deliberately does not import ``tdcrpy``. Every calculation, and
the initial enquiry that fills the widgets, runs in a separate process
(:mod:`tdcrpy._gui_worker`), so the window opens at once instead of waiting
for numba, stays responsive during a run that may take minutes, can cancel a
run by killing it, and survives a native crash in the kernel.

Nothing here keeps its own copy of the cocktail list, the nuclide list, the
configuration defaults, the typical light yields, the version or the
licence: all are read from the installed package. A hardcoded copy would be
one more thing to go stale, which this package has been bitten by
repeatedly.

Updating is handed to :mod:`tdcrpy._gui_updater`, run from a copy after the
window has closed: on Windows pip cannot replace files a running process
holds, and a failed attempt leaves a half-uninstalled package behind.
"""
import csv
import datetime
import json
import math
import os
import queue
import re
import shutil
import subprocess
import sys
import threading
import traceback
import webbrowser
import tkinter as tk
from importlib import metadata
from importlib.util import find_spec
from tkinter import filedialog, font as tkfont, messagebox, ttk

#: Matches the percentage tqdm prints, so the progress bar can follow the
#: Monte-Carlo loop without the model needing to report anything.
_TQDM_PCT = re.compile(r'(\d{1,3})%')

#: The model runs proper: each needs a source and a counting geometry.
CALC_TASKS = (
    ('efficiency', 'Efficiency at a given L'),
    ('fit', 'Fit L to a TDCR — stochastic'),
    ('fit_analytical', 'Fit L to a TDCR — analytical'),
    ('curve', 'Efficiency curve ε_D vs TDCR'),
)
#: Views of the data and sub-models the calculations rest on. None runs the
#: Monte Carlo, so they return in seconds.
DATA_TASKS = (
    ('decay_data', 'Decay data'),
    ('beta_spectrum', 'β spectrum'),
    ('stopping_power', 'Stopping power'),
    ('quenching', 'Ionisation quenching'),
)
TASKS = CALC_TASKS + DATA_TASKS

#: One sentence under the heading, saying what the selected item does.
TASK_HELP = {
    'efficiency': 'Detection efficiencies of the seven coincidence channels '
                  'for a light yield L you supply.',
    'fit': 'Fits the light yield L to a measured TDCR with the stochastic '
           'model, then reports the efficiencies at that L.',
    'fit_analytical': 'Fits L to a measured TDCR with the analytical model '
                      '— fast, valid for pure β emitters only.',
    'curve': 'ε_D against TDCR over a range of L: the decays are simulated '
             'once and replayed at every L, so the curve is smooth.',
    'decay_data': 'The decay scheme as the model reads it (PenNuc): decay '
                  'branches, γ rays and conversion electrons, per decay. '
                  'The Sources tab says where the data come from.',
    'beta_spectrum': 'The emitted energy spectrum of one β transition, from '
                     'the BetaShape archive the model samples. The Sources '
                     'tab gives its BetaShape version and citation.',
    'stopping_power': 'Stopping power dE/dx of electrons (every model side '
                      'by side, if wished) or α particles in the scintillator.',
    'quenching': 'Birks-quenched energy Q(E) and the quenching factor Q/E, '
                 'for one or several kB values.',
}

#: Most points any plotted grid (L, or energy) may have; on the efficiency
#: curve each one costs one replay.
CURVE_MAX_POINTS = 500

#: Default energy range per particle: β energies of interest run from eV to
#: MeV, α energies from 0.1 to 8 MeV (the end of the ASTAR table).
E_DEFAULTS = {'electron': '0.1, 2000, 100', 'alpha': '100, 8000, 80'}

#: Most kB values the quenching view overlays on one plot.
MAX_KB_SERIES = 8

#: Bounds on the number of nuclear-data realisations. Two is the least a
#: spread can come from; below about ten its upper bound is loose.
MIN_REPS, MAX_REPS = 2, 200

# Web resources. Every one also appears in README.md or setup.py.
REPO_URL = 'https://github.com/RomainCoulon/TDCRPy'
DOCS_URL = REPO_URL + '#readme'
NOTEBOOKS_URL = REPO_URL + '/tree/main/notebooks'
ISSUES_URL = REPO_URL + '/issues'
LICENCE_URL = REPO_URL + '/blob/main/LICENCE.md'
PYPI_URL = 'https://pypi.org/project/TDCRPy/'
PYPI_JSON_URL = 'https://pypi.org/pypi/TDCRPy/json'
HISTORY_URL = REPO_URL + '/commits/main'
PAPER_URL = 'https://doi.org/10.1016/j.apradiso.2024.111518'
REPORT_URL = 'http://dx.doi.org/10.13140/RG.2.2.15682.80321'
BIPM_URL = 'https://www.bipm.org'

# The palette of the TDCRPy logo.
NAVY = '#0B2A55'
CYAN = '#1EA7E1'
INK = '#1F2A37'
MUTED = '#5B6B80'
PANEL = '#F3F6FA'
RULE = '#D5DDE8'
SELECTED = '#DCEAF8'
WHITE = '#FFFFFF'
STRIPE = '#F6F8FB'
#: Series colours for plots, navy and cyan first.
PLOT_COLORS = ('#0B2A55', '#1EA7E1', '#E07A1F', '#2E9E5B', '#B03A48',
               '#7A5CC2', '#8C8C8C', '#C9A227')


class InputError(ValueError):
    """A problem with what the user typed, phrased for them to read."""


# --------------------------------------------------------------------------
# input handling, kept free of Tk so it can be tested headless
# --------------------------------------------------------------------------
def parse_numbers(text, field, count=None):
    """Parse a comma-separated list of numbers.

    Parameters
    ----------
    text : str
        What the user typed.
    field : str
        Field name, used in the error message.
    count : int or tuple of int, optional
        Permitted number of values. ``None`` accepts any non-empty list.

    Returns
    -------
    list of float
    """
    raw = [s.strip() for s in str(text).split(',') if s.strip()]
    if not raw:
        raise InputError(f'{field}: no value given.')
    try:
        vals = [float(s) for s in raw]
    except ValueError:
        bad = next(s for s in raw if not _is_number(s))
        raise InputError(f'{field}: {bad!r} is not a number.')
    if count is not None:
        allowed = (count,) if isinstance(count, int) else tuple(count)
        if len(vals) not in allowed:
            want = ' or '.join(str(a) for a in allowed)
            raise InputError(f'{field}: expected {want} value(s), '
                             f'got {len(vals)}.')
    return vals


def _is_number(s):
    try:
        float(s)
    except ValueError:
        return False
    return True


def parse_nuclides(text):
    """Split the nuclide field into labels, preserving the user's order."""
    rads = [s.strip() for s in str(text).split(',') if s.strip()]
    if not rads:
        raise InputError('Source: add at least one nuclide, e.g. Fe-55.')
    return rads


def parse_fractions(text, n_rad):
    """Parse activity fractions and check they match the nuclide count.

    A single nuclide needs no fraction, so an empty field becomes ``[1.0]``.
    Fractions are normalised by the library, but a set summing to something
    far from 1 is much more likely a typo than an intention, so it is
    refused here rather than silently rescaled.
    """
    text = str(text).strip()
    if not text:
        if n_rad == 1:
            return [1.0]
        raise InputError(f'Fractions: {n_rad} nuclides need {n_rad} '
                         f'fractions.')
    vals = parse_numbers(text, 'Fractions', count=n_rad)
    if any(v < 0 for v in vals):
        raise InputError('Fractions: a fraction cannot be negative.')
    total = sum(vals)
    if total <= 0:
        raise InputError('Fractions: they cannot all be zero.')
    if abs(total - 1.0) > 0.02:
        raise InputError(f'Fractions: they sum to {total:g}, not 1. '
                         f'Fix them, or give them as proportions of 1.')
    return vals


def add_to_source(source, rad, frac_text, known=()):
    """Add *rad* to the source mixture, or update its fraction.

    Parameters
    ----------
    source : list of (str, float)
        The mixture so far, in the order the user built it.
    rad : str
        The nuclide picked. Matched against *known* without regard to
        case, and stored under the archive's spelling.
    frac_text : str
        Its activity fraction, as typed.
    known : sequence of str
        The nuclides the archive has; empty skips the check (the list has
        not arrived yet).

    Returns
    -------
    list of (str, float)
        A new list; *source* is not modified.
    """
    rad = str(rad).strip()
    if not rad:
        raise InputError('Nuclide: pick one from the list.')
    if known:
        match = {k.lower(): k for k in known}.get(rad.lower())
        if match is None:
            raise InputError(f'Nuclide: {rad!r} is not in the decay-data '
                             f'archive. Pick one from the list.')
        rad = match
    frac = parse_numbers(frac_text or '1', 'Fraction', count=1)[0]
    if frac <= 0:
        raise InputError('Fraction: must be greater than zero.')
    out = list(source)
    for i, (r, _f) in enumerate(out):
        if r == rad:
            out[i] = (rad, frac)
            return out
    return out + [(rad, frac)]


def L_default_for(cocktail, defaults):
    """The default light yield for *cocktail*, and a note saying where it
    comes from.

    *defaults* is ``describe()``'s ``L_defaults``: ``{cocktail: [L, spread,
    measured]}``. Returns ``(None, '')`` when the library has no table.
    """
    entry = (defaults or {}).get(cocktail)
    if not entry:
        return None, ''
    L, spread, measured = entry
    if measured:
        return f'{L:g}', (f'typical for {cocktail}: {L:g} ± {spread:g} '
                          f'keV⁻¹ (measured)')
    return f'{L:g}', (f'no measured value for {cocktail}: a stand-in '
                      f'from the Ultima Gold family')


def normalise_source(source):
    """Scale the fractions so that they sum to 1."""
    total = sum(f for _r, f in source)
    if total <= 0:
        return list(source)
    return [(r, f / total) for r, f in source]


def source_fields(source):
    """The mixture as the comma-separated 'rad' and 'pmf' a job carries.

    A single nuclide's fraction is 1 whatever was entered: there is
    nothing for it to be a fraction of.
    """
    if not source:
        return '', ''
    rad = ', '.join(r for r, _f in source)
    if len(source) == 1:
        return rad, '1'
    return rad, ', '.join(f'{f:.6g}' for _r, f in source)


def _build_grid(text, field, log_spacing):
    """Turn "from, to, points" into an evenly- or log-spaced grid.

    Log spacing is the default because the functions plotted this way --
    efficiency vs L, stopping power vs E -- change fastest at the low end:
    an even grid would spend most of its points on the flat part.
    """
    lo, hi, n = parse_numbers(text, field, count=3)
    if n != int(n) or not 2 <= n <= CURVE_MAX_POINTS:
        raise InputError(f'{field}: the number of points must be a whole '
                         f'number from 2 to {CURVE_MAX_POINTS}.')
    if not 0 < lo < hi:
        raise InputError(f'{field}: need 0 < from < to.')
    n = int(n)
    if log_spacing:
        r = (hi / lo) ** (1.0 / (n - 1))
        grid = [lo * r ** i for i in range(n)]
    else:
        grid = [lo + (hi - lo) * i / (n - 1) for i in range(n)]
    grid[-1] = hi                    # no rounding drift on the end point
    return grid


def build_L_grid(text, log_spacing=True):
    """Light yields (keV⁻¹) to evaluate the efficiency curve at."""
    return _build_grid(text, 'L range', log_spacing)


def build_E_grid(text, log_spacing=True):
    """Energies (keV) to evaluate the stopping power or quenching at."""
    return _build_grid(text, 'Energy range', log_spacing)


def parse_particle(text):
    """The particle a stopping-power or quenching view is for."""
    particle = str(text or 'electron').strip().lower()
    if particle not in E_DEFAULTS:
        raise InputError(f'Particle: {text!r} is neither an electron nor an '
                         f'alpha particle.')
    return particle


def build_job(form):
    """Turn the form's string values into a job for the worker.

    Raises
    ------
    InputError
        With a message meant for the user, if anything is missing or
        inconsistent. Checking here rather than in the worker means a typo
        costs no process start-up and no numba import.
    """
    task = form['task']
    config = form.get('config') or {}

    # ---- data views: no source mixture, no counting geometry ------------
    if task == 'stopping_power':
        # dE/dx is a property of the scintillator (set via 'config') and the
        # model alone.
        particle = parse_particle(form.get('particle'))
        return {
            'task': task,
            'particle': particle,
            'E_grid': build_E_grid(form.get('Erange', ''),
                                   bool(form.get('Elog', True))),
            'compare_models': (particle == 'electron'
                               and bool(form.get('spCompare', False))),
            'config': config,
        }

    if task == 'quenching':
        kBs = parse_numbers(form.get('kBq') or form.get('kB', ''),
                            'Birks constant(s) kB')
        if any(k < 0 for k in kBs):
            raise InputError('Birks constant(s) kB: cannot be negative.')
        if len(kBs) > MAX_KB_SERIES:
            raise InputError(f'Birks constant(s) kB: at most {MAX_KB_SERIES} '
                             f'values on one plot.')
        return {
            'task': task,
            'particle': parse_particle(form.get('particle')),
            'E_grid': build_E_grid(form.get('Erange', ''),
                                   bool(form.get('Elog', True))),
            'kB_list': kBs,
            'config': config,
        }

    if task == 'beta_spectrum':
        rad = parse_nuclides(form.get('data_rad') or form['rad'])[0]
        return {
            'task': task,
            'rad': rad,
            'mode': form.get('betaMode') or 'beta-',
            'level': (form.get('betaLevel') or 'tot').strip(),
            'config': config,
        }

    if task == 'decay_data':
        return {'task': task,
                'rad': parse_nuclides(form.get('data_rad') or form['rad'])[0],
                'config': config}

    # ---- the model runs --------------------------------------------------
    rads = parse_nuclides(form['rad'])
    fracs = parse_fractions(form.get('pmf', ''), len(rads))

    kB = parse_numbers(form['kB'], 'Birks constant kB', count=1)[0]
    if kB < 0:
        raise InputError('Birks constant kB: cannot be negative.')
    V = parse_numbers(form['V'], 'Volume', count=1)[0]
    if V <= 0:
        raise InputError('Volume: must be greater than zero.')

    job = {
        'task': task,
        'rad': ', '.join(rads),
        'pmf': ', '.join(f'{v:g}' for v in fracs),
        'kB': kB,
        'V': V,
        'display': bool(form.get('display')),
        'config': config,
    }

    # Nuclear-data uncertainty: only where repeated runs can measure it.
    stochastic = task == 'fit' or (task == 'efficiency'
                                   and bool(form.get('Smodel', True)))
    if bool(form.get('uncData')) and stochastic:
        M = parse_numbers(form.get('uncReps') or '10', 'Realisations M',
                          count=1)[0]
        if M != int(M) or not MIN_REPS <= M <= MAX_REPS:
            raise InputError(f'Realisations M: a whole number from '
                             f'{MIN_REPS} to {MAX_REPS}.')
        job['uncData'] = True
        job['uncReps'] = int(M)

    if task == 'efficiency':
        # One value is a symmetric counter; three set the PMTs separately.
        L = parse_numbers(form['L'], 'Light yield L', count=(1, 3))
        job['L'] = L[0] if len(L) == 1 else L
        job['Smodel'] = bool(form.get('Smodel', True))
    elif task == 'curve':
        job['L_grid'] = build_L_grid(form.get('Lrange', ''),
                                     bool(form.get('Llog', True)))
        job['Smodel'] = bool(form.get('Smodel', True))
    else:
        # One value is T/D; four add T/AB, T/BC and T/AC for an asymmetric
        # counter, which is the only way to fit three separate yields.
        TD = parse_numbers(form['TD'], 'Measured TDCR', count=(1, 4))
        for v in TD:
            if not 0 < v <= 1:
                raise InputError(f'Measured TDCR: {v:g} is not a ratio '
                                 f'between 0 and 1.')
        job['TD'] = TD[0] if len(TD) == 1 else TD
        job['L_guess'] = parse_numbers(
            form.get('L_guess') or '1', 'Initial L guess', count=1)[0]
        bounds = str(form.get('Lbounds', '')).strip()
        if bounds:
            lo, hi = parse_numbers(bounds, 'L bounds', count=2)
            if lo >= hi:
                raise InputError('L bounds: the lower bound must be below '
                                 'the upper one.')
            job['Lbounds'] = [lo, hi]

    if task != 'fit_analytical':
        # The analytical curve ignores N, but TDCRPy() still takes one.
        N = parse_numbers(form['N'], 'Number of decays N', count=1)[0]
        if N < 1:
            raise InputError('Number of decays N: must be at least 1.')
        job['N'] = int(N)

    return job


def format_value(value, uncertainty):
    """Render a value and its uncertainty for the results table."""
    if isinstance(value, str):
        return value, '', ''
    if value is None:
        return '', '', ''
    if isinstance(value, (list, tuple)):
        # Three light yields, one per PMT, for an asymmetric counter.
        return ', '.join(f'{v:.6g}' for v in value), '', ''
    val = f'{value:.6g}'
    if uncertainty in (None, ''):
        return val, '', ''
    unc = f'{uncertainty:.3g}'
    rel = f'{100.0 * uncertainty / value:.2f} %' if value else ''
    return val, unc, rel


def format_cell(value, digits=5):
    """One cell of a curve or data table, for reading."""
    if value is None:
        return ''
    if isinstance(value, str):
        return value
    return f'{value:.{digits}g}'


def is_uncertainty_column(key):
    """Columns holding an uncertainty: shown to 3 significant digits."""
    return key == 'u' or key.startswith('u_')


def numberish(text):
    """Whether a string cell reads as a number ('< 0.03', '0.52 %'), and so
    belongs right-aligned with the numbers."""
    s = str(text).strip().lstrip('<').rstrip('%').strip()
    return _is_number(s) if s else False


def write_csv(path, rows, meta=()):
    """Write a results table, its header first, with ``# ``-prefixed
    provenance lines above it.

    UTF-8 with a byte-order mark, so a spreadsheet opens β, γ and ε_D
    correctly; ``pandas.read_csv(path, comment='#')`` skips the
    provenance.
    """
    with open(path, 'w', encoding='utf-8-sig', newline='') as fh:
        for line in meta:
            fh.write(f'# {line}\n')
        csv.writer(fh).writerows(rows)


def provenance_lines(info):
    """A result's provenance sections as plain lines, for a CSV header."""
    lines = []
    for section in info or ():
        lines.append(section['title'])
        for label, text, _url in section['items']:
            lines.append(f'  {label}: {text}' if label else f'    {text}')
    return lines


# --------------------------------------------------------------------------
# the decay view: wording for the records of TDCRPy(trace=)
# --------------------------------------------------------------------------
_X_RAY_LABEL = {'XKA': 'X-ray Kα', 'XKB': 'X-ray Kβ', 'XL': 'X-ray L'}
_EMITTER_LABEL = {'beta': 'β⁻', 'beta+': 'β⁺', 'alpha': 'α', 'gamma': 'γ',
                  'electron': 'conversion e⁻',
                  'annihilation': 'annihilation γ'}
#: What the interaction step did to a particle, as the view says it.
FATE_LABEL = {'escaped': 'escaped the vial', 'compton': 'Compton scattered',
              'photoelectric': 'photoabsorbed', 'pair': 'pair production',
              'partial': 'partly absorbed (wall)', 'absorbed': ''}


def particle_label(origin, kind=None):
    """How one particle of a trace reads: its origin as emitted (β⁻, Auger
    e⁻ (KLL), X-ray Kα...), or what a secondary turned into."""
    o = origin or ''
    if o in _EMITTER_LABEL:
        return _EMITTER_LABEL[o]
    if o == 'secondary':
        return 'pair e⁺' if kind == 'positron' else 'secondary e⁻'
    if o.startswith('Auger'):
        rest = o[len('Auger'):].strip()
        return f'Auger e⁻ ({rest})' if rest else 'Auger e⁻'
    if o in _X_RAY_LABEL:
        return _X_RAY_LABEL[o]
    if o.startswith('X'):
        return f'X-ray {o[1:]}'
    return o


def nuclide_name(code):
    """PenNuc's upper-case daughter code as a nuclide name: NI60 -> Ni-60."""
    m = re.fullmatch(r'([A-Z]{1,2})(\d+)([A-Z]*)', code or '')
    if not m:
        return code or ''
    sym, mass, tail = m.groups()
    return f'{sym.capitalize()}-{mass}{tail.lower()}'


def format_time(ns):
    """A delay in ns, in the unit that reads best."""
    for scale, unit in ((1e9, 's'), (1e6, 'ms'), (1e3, 'µs')):
        if ns >= scale:
            return f'{ns / scale:.3g} {unit}'
    return f'{ns:.3g} ns'


def _keV(value):
    return '?' if value is None else f'{value:.6g} keV'


def step_label(step):
    """One sampled step of the decay (branch or transition) as a sentence."""
    if step.get('kind') == 'branch':
        p, level = step.get('particle'), step.get('level')
        to = f'level {level:g}' if level else 'ground state'
        if p is None:
            return ('isomeric transition: no particle emitted' if level
                    else 'no particle and no level: nothing is deposited')
        if p.startswith('Atom_'):
            return (f'electron capture ({p[5:]} shell), '
                    f'{_keV(step.get("energy"))} → {to}')
        name = _EMITTER_LABEL.get(p, p)
        what = 'end point ' if p in ('beta', 'beta+') else ''
        return f'{name} emission, {what}{_keV(step.get("energy"))} → {to}'
    if step.get('kind') == 'transition':
        t = step.get('type', '')
        name = 'γ' if t == 'GA' else f'conversion e⁻ ({t[1:]} shell)'
        text = f'{name} {_keV(step.get("energy"))}'
        if step.get('from_keV') is not None:
            text += f', from the level at {step["from_keV"]:.6g} keV'
        if step.get('window') == 'delayed':
            text += f'  —  delayed, t = {format_time(step.get("t_ns", 0.0))}'
        return text
    return str(step)


def photoelectrons(rec, Q):
    """Mean photoelectron number per PMT (A, B, C) for light *Q* keV, as the
    semi-analytical detection model computes it: ``L · μ · Q / 3``. ``None``
    for the optical model, which samples photons instead."""
    if rec.get('optical'):
        return None
    mu = rec.get('mu') or []
    L = rec.get('L')
    if isinstance(L, (list, tuple)):
        return [l * m * Q / 3.0 for l, m in zip(L, mu)]
    m = sum(mu) / len(mu) if mu else 1.0
    return [L * m * Q / 3.0] * 3


def counted_windows(rec):
    """The windows the detection model counts: the prompt one, plus the
    delayed one when it is a pulse of its own."""
    windows = rec.get('windows') or []
    return windows[:max(1, int(rec.get('pulses', 1)))]


def blend(c0, c1, f):
    """Mix two #rrggbb colours: f = 0 gives c0, 1 gives c1."""
    f = min(1.0, max(0.0, f))
    a = [int(c0[i:i + 2], 16) for i in (1, 3, 5)]
    b = [int(c1[i:i + 2], 16) for i in (1, 3, 5)]
    return '#' + ''.join(f'{round(x + (y - x) * f):02x}' for x, y in zip(a, b))


def package_dir():
    """Locate the installed package without importing it.

    ``find_spec`` reads the metadata but does not execute the module, so the
    window can find ``config.toml`` without paying for the numba import.
    """
    spec = find_spec('tdcrpy')
    if spec is None or not spec.submodule_search_locations:
        return None
    return spec.submodule_search_locations[0]


def installed_version():
    """The installed TDCRPy version, from its metadata."""
    try:
        return metadata.version('TDCRPy')
    except metadata.PackageNotFoundError:
        return '?'


def licence_text():
    """The licence shipped with the installed package, or ``None``.

    Read from the distribution's files; an editable install keeps it beside
    the package instead.
    """
    try:
        for f in metadata.distribution('TDCRPy').files or ():
            if f.name.upper().startswith('LICEN'):
                return f.read_text(encoding='utf-8')
    except (metadata.PackageNotFoundError, OSError):
        pass
    d = package_dir()
    if d:
        path = os.path.join(os.path.dirname(d), 'LICENCE.md')
        if os.path.exists(path):
            with open(path, encoding='utf-8') as fh:
                return fh.read()
    return None


def task_label(key):
    return dict(TASKS).get(key, key)


# --------------------------------------------------------------------------
# updates
# --------------------------------------------------------------------------
def version_key(version):
    """The numeric release part of a version, for comparison.

    ``2.20.36`` -> ``(2, 20, 36)``. A suffix (``rc1``, ``.dev0``) is ignored,
    which is all this needs: PyPI's ``info.version`` is the latest stable
    release.
    """
    m = re.match(r'\s*v?(\d+(?:\.\d+)*)', str(version))
    return tuple(int(p) for p in m.group(1).split('.')) if m else ()


def is_newer(latest, current):
    """Whether *latest* is a later release than *current*."""
    a, b = version_key(latest), version_key(current)
    return bool(a) and bool(b) and a > b


def fetch_latest_release(timeout=10.0):
    """Ask PyPI for the latest TDCRPy release: ``(version, upload date)``.

    Only the package name goes out; the user-agent names the client and
    nothing about the user.
    """
    import urllib.request
    req = urllib.request.Request(PYPI_JSON_URL, headers={
        'User-Agent': 'TDCRPy-gui', 'Accept': 'application/json'})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        data = json.loads(resp.read().decode('utf-8'))
    latest = data['info']['version']
    files = data.get('releases', {}).get(latest) or [{}]
    date = (files[0].get('upload_time') or '')[:10]
    return latest, date


def install_kind():
    """How this TDCRPy was installed: ``'pip'``, ``'editable'`` (a source
    checkout, updated with git), another installer's name (``'conda'``...),
    or ``'unknown'``."""
    try:
        dist = metadata.distribution('TDCRPy')
    except metadata.PackageNotFoundError:
        return 'unknown'
    try:
        direct = json.loads(dist.read_text('direct_url.json') or '{}')
    except ValueError:
        direct = {}
    if direct.get('dir_info', {}).get('editable'):
        return 'editable'
    installer = (dist.read_text('INSTALLER') or '').strip().lower()
    if installer in ('', 'pip', 'uv'):
        return 'pip'
    return installer


def console_python():
    """The interpreter running this window, as python.exe: pythonw.exe (what
    the gui launcher uses) has no console for pip to report to."""
    exe = sys.executable
    if os.name == 'nt' and os.path.basename(exe).lower() == 'pythonw.exe':
        alt = os.path.join(os.path.dirname(exe), 'python.exe')
        if os.path.exists(alt):
            return alt
    return exe


def launch_updater(version, current):
    """Start the update helper on a copy of itself, then let the window
    close: see :mod:`tdcrpy._gui_updater` for why it must run after.

    Returns the path of the log the helper writes.
    """
    import sysconfig
    import tempfile
    pkg = package_dir()
    tmp = tempfile.mkdtemp(prefix='tdcrpy-update-')
    helper = os.path.join(tmp, 'tdcrpy_update.py')
    shutil.copyfile(os.path.join(pkg, '_gui_updater.py'), helper)
    log = os.path.join(tmp, 'update.log')
    python = console_python()
    runner = python
    if os.name == 'nt':
        w = os.path.join(os.path.dirname(python), 'pythonw.exe')
        runner = w if os.path.exists(w) else python
    launcher = os.path.join(sysconfig.get_path('scripts') or '',
                            'tdcrpy-gui.exe' if os.name == 'nt'
                            else 'tdcrpy-gui')
    cmd = [runner, '-I', helper, '--pid', str(os.getpid()),
           '--python', python, '--spec', f'TDCRPy=={version}',
           '--from-version', current, '--package-dir', pkg,
           '--log', log, '--relaunch']
    if os.path.exists(launcher):
        cmd += ['--launcher', launcher]
    flags = 0
    if os.name == 'nt':
        flags = (getattr(subprocess, 'DETACHED_PROCESS', 0)
                 | getattr(subprocess, 'CREATE_NEW_PROCESS_GROUP', 0))
    subprocess.Popen(cmd, close_fds=True, creationflags=flags, cwd=tmp)
    return log


# --------------------------------------------------------------------------
# the worker, seen from the window
# --------------------------------------------------------------------------
class Worker:
    """A single worker process, with its output pumped into a queue.

    Both pipes are drained by threads: a full pipe buffer would otherwise
    block the model mid-run, and the log can be verbose with ``display``.
    """

    def __init__(self, args, stdin_payload, out_queue):
        self._queue = out_queue
        creation = 0
        if os.name == 'nt':
            # Keeps a console window from flashing up for each run.
            creation = getattr(subprocess, 'CREATE_NO_WINDOW', 0)
        self.proc = subprocess.Popen(
            [sys.executable, '-W', 'ignore', '-u',
             '-m', 'tdcrpy._gui_worker'] + list(args),
            stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, text=True, encoding='utf-8',
            errors='replace', bufsize=1, creationflags=creation,
        )
        threading.Thread(target=self._pump_stdout, daemon=True).start()
        threading.Thread(target=self._pump_stderr, daemon=True).start()
        threading.Thread(target=self._wait, args=(stdin_payload,),
                         daemon=True).start()

    def _wait(self, payload):
        try:
            if payload is not None:
                self.proc.stdin.write(payload)
            self.proc.stdin.close()
        except OSError:
            pass
        rc = self.proc.wait()
        self._queue.put(('exit', rc))

    def _pump_stdout(self):
        for line in self.proc.stdout:
            line = line.strip()
            if not line:
                continue
            try:
                self._queue.put(('msg', json.loads(line)))
            except ValueError:
                # Not protocol: show it rather than drop it.
                self._queue.put(('log', line))

    def _pump_stderr(self):
        for line in self.proc.stderr:
            self._queue.put(('log', line.rstrip()))

    def cancel(self):
        if self.proc.poll() is None:
            self.proc.kill()


# --------------------------------------------------------------------------
# the window
# --------------------------------------------------------------------------
class App(ttk.Frame):

    def __init__(self, master):
        super().__init__(master)
        self.version = installed_version()
        # Pixels per 96-dpi pixel: images, column widths and the minimum
        # size are in pixels and must follow the screen; fonts already do.
        self.s = max(1.0, master.winfo_fpixels('1i') / 96.0)
        master.title(f'TDCRPy {self.version} — TDCR model for liquid '
                     f'scintillation counting')
        self.grid(sticky='nsew')
        master.columnconfigure(0, weight=1)
        master.rowconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(1, weight=1)

        self._queue = queue.Queue()
        self._worker = None
        self._describe = None
        self._rows = []
        self._config_snapshot = None
        self._cancelled = False
        self._last_job = None
        self._result_job = None
        self._figure = None
        self._canvas = None
        self._L_shown = {}           # the default last written per L field
        self._result_info = None
        self._update_check = None
        self._trace = []             # decay records for the Decay view
        self._trace_i = 0
        self._trace_job = None

        self.vars = {}
        self._setup_style()
        self._load_images()
        self._build_menu()
        self._build_header()
        self._build_sidebar()
        self._build_body()
        self._build_actions()
        self._build_statusbar()
        self._bind_keys()

        master.protocol('WM_DELETE_WINDOW', self._quit)
        master.report_callback_exception = self._on_tk_error

        self._snapshot_config()
        self._on_task_change()
        self.after(80, self._drain)
        self._enquire()

    # ---------------------------------------------------------------- setup
    def _px(self, n):
        return int(round(n * self.s))

    def _var(self, name, value=''):
        v = tk.StringVar(value=str(value))
        self.vars[name] = v
        return v

    def _flag(self, name, value=False):
        v = tk.BooleanVar(value=bool(value))
        self.vars[name] = v
        return v

    def _setup_style(self):
        base = tkfont.nametofont('TkDefaultFont')
        family = base.actual('family')
        size = base.actual('size')
        self.f_title = tkfont.Font(family=family, size=15, weight='bold')
        self.f_sub = tkfont.Font(family=family, size=9)
        self.f_heading = tkfont.Font(family=family, size=13, weight='bold')
        self.f_section = tkfont.Font(family=family, size=8, weight='bold')
        self.f_nav = tkfont.Font(family=family, size=max(10, size))
        self.f_bold = tkfont.Font(family=family, size=max(9, size),
                                  weight='bold')

        style = ttk.Style(self.master)
        if os.name != 'nt' and 'clam' in style.theme_names():
            style.theme_use('clam')   # the X11 'default' theme looks dated
        style.configure('Treeview', rowheight=self._px(22))
        style.configure('Treeview.Heading', font=self.f_bold)
        style.configure('Hint.TLabel', foreground=MUTED)
        style.configure('Section.TLabel', foreground=NAVY, font=self.f_bold)
        style.configure('Heading.TLabel', foreground=NAVY,
                        font=self.f_heading)
        style.configure('Status.TLabel', foreground=MUTED)

    def _load_images(self):
        """The logos and window icon from ``tdcrpy/assets``; any missing
        file just leaves its place empty."""
        self.img = {}
        d = package_dir()
        assets = os.path.join(d, 'assets') if d else None
        hidpi = self.s >= 1.75

        def load(name):
            if not assets:
                return None
            for candidate in ((f'{name}_2x.png', f'{name}.png') if hidpi
                              else (f'{name}.png',)):
                path = os.path.join(assets, candidate)
                if os.path.exists(path):
                    try:
                        return tk.PhotoImage(file=path)
                    except tk.TclError:
                        return None
            return None

        for name in ('tdcrpy_logo_header', 'tdcrpy_logo_about', 'bipm_logo'):
            self.img[name] = load(name)
        icons = [load(f'tdcrpy_icon_{n}') for n in (256, 48, 32, 16)]
        icons = [i for i in icons if i is not None]
        if icons:
            try:
                self.master.iconphoto(True, *icons)
            except tk.TclError:
                pass
        self.img['icons'] = icons

    def _build_menu(self):
        bar = tk.Menu(self.master)
        m = tk.Menu(bar, tearoff=False)
        m.add_command(label='Export results as CSV…', accelerator='Ctrl+S',
                      command=self._save)
        m.add_command(label='Save plot as image…', command=self._save_plot)
        m.add_separator()
        m.add_command(label='Exit', accelerator='Ctrl+Q', command=self._quit)
        bar.add_cascade(label='File', menu=m)

        m = tk.Menu(bar, tearoff=False)
        m.add_command(label='Run', accelerator='F5', command=self._run)
        m.add_command(label='Cancel', accelerator='Esc', command=self._cancel)
        bar.add_cascade(label='Run', menu=m)

        m = tk.Menu(bar, tearoff=False)
        m.add_command(label='Undo my configuration changes',
                      command=self._restore_snapshot)
        m.add_command(label='Restore shipped defaults',
                      command=self._restore_defaults)
        m.add_separator()
        m.add_command(label='Open the configuration folder',
                      command=self._open_config_folder)
        bar.add_cascade(label='Tools', menu=m)

        m = tk.Menu(bar, tearoff=False)
        m.add_command(label='Documentation (online)',
                      command=lambda: self._open(DOCS_URL))
        m.add_command(label='Example notebooks',
                      command=lambda: self._open(NOTEBOOKS_URL))
        m.add_separator()
        m.add_command(label='Reference paper — Appl. Radiat. Isot. (2024)',
                      command=lambda: self._open(PAPER_URL))
        m.add_command(label='BIPM technical report',
                      command=lambda: self._open(REPORT_URL))
        m.add_separator()
        m.add_command(label='Report an issue…',
                      command=lambda: self._open(ISSUES_URL))
        m.add_command(label='TDCRPy on PyPI',
                      command=lambda: self._open(PYPI_URL))
        m.add_command(label='Check for updates…',
                      command=self._check_updates)
        m.add_separator()
        m.add_command(label='Licence', command=self._show_licence)
        m.add_command(label='About TDCRPy', command=self._show_about)
        bar.add_cascade(label='Help', menu=m)
        self.master.configure(menu=bar)

    def _bind_keys(self):
        root = self.master
        root.bind('<F5>', lambda e: self._run())
        root.bind('<Escape>', lambda e: self._cancel())
        root.bind('<Control-s>', lambda e: self._save())
        root.bind('<Control-q>', lambda e: self._quit())

    # ------------------------------------------------------------- header
    def _build_header(self):
        head = tk.Frame(self, bg=WHITE)
        head.grid(row=0, column=0, columnspan=2, sticky='ew')
        head.columnconfigure(2, weight=1)

        logo = self.img.get('tdcrpy_logo_header')
        if logo is not None:
            mark = tk.Label(head, image=logo, bg=WHITE, cursor='hand2')
            mark.grid(row=0, column=0, rowspan=2,
                      padx=(self._px(14), self._px(12)), pady=self._px(8))
            mark.bind('<Button-1>', lambda e: self._show_about())
        else:
            tk.Label(head, text='TDCRPy', font=self.f_title, fg=NAVY,
                     bg=WHITE).grid(row=0, column=0, rowspan=2,
                                    padx=self._px(14), pady=self._px(8))

        tk.Label(head, text='Liquid scintillation counting — TDCR model',
                 font=self.f_title, fg=NAVY, bg=WHITE, anchor='w').grid(
            row=0, column=1, sticky='sw')
        tk.Label(head, text=f'Version {self.version}   ·   '
                            f'Photo-physical stochastic model   ·   '
                            f'MIT licence',
                 font=self.f_sub, fg=MUTED, bg=WHITE, anchor='w').grid(
            row=1, column=1, sticky='nw')

        bipm = self.img.get('bipm_logo')
        if bipm is not None:
            mark = tk.Label(head, image=bipm, bg=WHITE, cursor='hand2')
        else:
            mark = tk.Label(head, text='BIPM', font=self.f_title, fg=NAVY,
                            bg=WHITE, cursor='hand2')
        mark.grid(row=0, column=3, rowspan=2, padx=self._px(16))
        mark.bind('<Button-1>', lambda e: self._open(BIPM_URL))
        _Tooltip(mark, 'Developed and maintained by the Bureau International '
                       'des Poids et Mesures — open bipm.org')

        tk.Frame(self, bg=CYAN, height=self._px(3)).grid(
            row=0, column=0, columnspan=2, sticky='sew')

    # ------------------------------------------------------------ sidebar
    def _build_sidebar(self):
        side = tk.Frame(self, bg=PANEL, width=self._px(220))
        side.grid(row=1, column=0, rowspan=2, sticky='nsw')
        side.grid_propagate(False)
        side.columnconfigure(0, weight=1)
        self._var('task', TASKS[0][0])

        r = 0
        for title, tasks in (('CALCULATIONS', CALC_TASKS),
                             ('DATA & PLOTS', DATA_TASKS)):
            tk.Label(side, text=title, font=self.f_section, fg=MUTED,
                     bg=PANEL, anchor='w').grid(
                row=r, column=0, sticky='ew', padx=self._px(16),
                pady=(self._px(16), self._px(4)))
            r += 1
            for key, label in tasks:
                tk.Radiobutton(
                    side, text=label, value=key, variable=self.vars['task'],
                    command=self._on_task_change, indicatoron=False,
                    anchor='w', justify='left', font=self.f_nav,
                    fg=INK, bg=PANEL, activebackground=SELECTED,
                    selectcolor=SELECTED, relief='flat', offrelief='flat',
                    bd=0, highlightthickness=0, cursor='hand2',
                    padx=self._px(18), pady=self._px(5),
                ).grid(row=r, column=0, sticky='ew')
                r += 1

        side.rowconfigure(r, weight=1)
        r += 1
        self.side_help = tk.Label(
            side, font=self.f_sub, fg=MUTED, bg=PANEL, anchor='nw',
            justify='left', wraplength=self._px(190),
            text='F5 runs the selected item.\nEvery result can be exported '
                 'as CSV, and every plot saved as PNG, PDF or SVG.')
        self.side_help.grid(row=r, column=0, sticky='ew', padx=self._px(16),
                            pady=self._px(12))

    # --------------------------------------------------------------- body
    def _build_body(self):
        pane = ttk.PanedWindow(self, orient='horizontal')
        pane.grid(row=1, column=1, sticky='nsew', padx=self._px(10),
                  pady=(self._px(10), 0))

        left = ttk.Frame(pane)
        left.columnconfigure(0, weight=1)
        left.rowconfigure(2, weight=1)
        self.heading = ttk.Label(left, text='', style='Heading.TLabel')
        self.heading.grid(row=0, column=0, sticky='w')
        self.subheading = ttk.Label(left, text='', style='Hint.TLabel',
                                    wraplength=self._px(470), justify='left')
        self.subheading.grid(row=1, column=0, sticky='w', pady=(0, 6))
        self.inputs = ttk.Notebook(left)
        self.inputs.grid(row=2, column=0, sticky='nsew')
        self._tab_source(self.inputs)
        self._tab_data(self.inputs)
        self._tab_scintillator(self.inputs)
        self._tab_micelles(self.inputs)
        self._tab_delayed(self.inputs)
        self._tab_counter(self.inputs)
        pane.add(left, weight=1)

        right = ttk.Frame(pane)
        right.columnconfigure(0, weight=1)
        right.rowconfigure(0, weight=1)
        self.results = ttk.Notebook(right)
        self.results.grid(row=0, column=0, sticky='nsew')
        self._tab_table(self.results)
        self._tab_plot(self.results)
        self._tab_decay(self.results)
        self._tab_sources(self.results)
        self._tab_log(self.results)
        pane.add(right, weight=2)

    def _row(self, parent, r, label, name, default='', width=16, hint=None):
        ttk.Label(parent, text=label).grid(row=r, column=0, sticky='w',
                                           pady=2)
        e = ttk.Entry(parent, textvariable=self._var(name, default),
                      width=width)
        e.grid(row=r, column=1, sticky='w', pady=2)
        if hint:
            ttk.Label(parent, text=hint, style='Hint.TLabel').grid(
                row=r, column=2, sticky='w', padx=(8, 0))
        return e

    def _section(self, parent, r, text):
        ttk.Label(parent, text=text, style='Section.TLabel').grid(
            row=r, column=0, columnspan=3, sticky='w', pady=(10, 2))

    def _tab_source(self, nb):
        f = ttk.Frame(nb, padding=12)
        nb.add(f, text='Source')
        self.tab_source = f
        # The job still reads 'rad' and 'pmf' as comma-separated text; the
        # list below is the only thing that writes them.
        self._var('rad', 'Fe-55')
        self._var('pmf', '1')
        self.source = [('Fe-55', 1.0)]

        self._section(f, 0, 'Radionuclide source')
        pick = ttk.Frame(f)
        pick.grid(row=1, column=0, columnspan=3, sticky='w')
        ttk.Label(pick, text='Nuclide').pack(side='left')
        self.nuclide_box = ttk.Combobox(
            pick, textvariable=self._var('pick_rad', 'Fe-55'), width=12,
            values=())
        self.nuclide_box.pack(side='left', padx=(6, 12))
        self.nuclide_box.bind('<Return>', lambda e: self._add_nuclide())
        self.nuclide_box.bind('<<ComboboxSelected>>',
                              lambda e: self.frac_entry.focus_set())
        ttk.Label(pick, text='Fraction').pack(side='left')
        self.frac_entry = ttk.Entry(
            pick, textvariable=self._var('pick_frac', '1'), width=8)
        self.frac_entry.pack(side='left', padx=(6, 8))
        self.frac_entry.bind('<Return>', lambda e: self._add_nuclide())
        self.add_button = ttk.Button(pick, text='Add', width=6,
                                     command=self._add_nuclide)
        self.add_button.pack(side='left')
        _Tooltip(self.add_button, 'Add the nuclide to the source, or update '
                                  'its fraction if it is already there')

        lst = ttk.Frame(f)
        lst.grid(row=2, column=0, columnspan=3, sticky='ew', pady=(6, 0))
        self.source_tree = ttk.Treeview(
            lst, columns=('rad', 'frac'), show='headings', height=3,
            selectmode='browse')
        self.source_tree.heading('rad', text='Nuclide')
        self.source_tree.heading('frac', text='Activity fraction')
        self.source_tree.column('rad', width=self._px(110), anchor='w')
        self.source_tree.column('frac', width=self._px(120), anchor='e')
        self.source_tree.pack(side='left')
        self.source_tree.bind('<<TreeviewSelect>>', self._on_source_select)
        self.source_tree.bind('<Delete>', lambda e: self._remove_nuclide())
        side = ttk.Frame(lst)
        side.pack(side='left', fill='y', padx=(8, 0))
        self.source_buttons = []
        for text, cmd, tip in (
                ('Remove', self._remove_nuclide, 'Remove the selected nuclide'),
                ('Normalise', self._normalise_source,
                 'Scale the fractions so that they sum to 1 — enter '
                 'activities or ratios, then normalise'),
                ('Clear', self._clear_source, 'Empty the source')):
            b = ttk.Button(side, text=text, width=10, command=cmd)
            b.pack(anchor='w', pady=(0, 3))
            _Tooltip(b, tip)
            self.source_buttons.append(b)
        self.source_sum = ttk.Label(f, style='Hint.TLabel')
        self.source_sum.grid(row=3, column=0, columnspan=3, sticky='w')
        self._refresh_source()

        self._section(f, 4, 'Light yield')
        self.l_entry = self._row(f, 5, 'Light yield L', 'L', '',
                                 hint='keV⁻¹; 3 values for asymmetric PMTs')
        self.L_hint = ttk.Label(f, style='Hint.TLabel')
        self.L_hint.grid(row=6, column=1, columnspan=2, sticky='w')
        self.td_entry = self._row(f, 7, 'Measured TDCR', 'TD', '0.9',
                                  hint='T/D; 4 values for an asymmetric '
                                       'counter')
        self.guess_entry = self._row(f, 8, 'Initial L guess', 'L_guess', '')
        self.bounds_entry = self._row(
            f, 9, 'L bounds', 'Lbounds', '',
            hint='blank lets the fit bracket the physical root itself')
        self.range_entry = self._row(
            f, 10, 'L range', 'Lrange', '0.1, 10, 25',
            hint='from, to, number of points (keV⁻¹)')
        self.log_check = ttk.Checkbutton(
            f, text='log spacing', variable=self._flag('Llog', True))
        self.log_check.grid(row=11, column=1, sticky='w')

        self._section(f, 12, 'Sample and statistics')
        self.kb_entry = self._row(f, 13, 'Birks constant kB', 'kB', '1.0e-5',
                                  hint='cm keV⁻¹')
        self.v_entry = self._row(f, 14, 'Volume', 'V', '10', hint='mL')
        self.n_entry = self._row(f, 15, 'Number of decays N', 'N', '10000',
                                 hint='≥ 10 000 recommended (JCGM 101)')

        self.smodel_check = ttk.Checkbutton(
            f, text='Stochastic model (uncheck for the analytical one, '
                    'valid for pure β emitters only)',
            variable=self._flag('Smodel', True),
            command=self._on_task_change)
        self.smodel_check.grid(row=16, column=0, columnspan=3, sticky='w',
                               pady=(8, 0))
        self.unc_check = ttk.Checkbutton(
            f, text='Propagate nuclear-data uncertainty',
            variable=self._flag('uncData'), command=self._on_task_change)
        self.unc_check.grid(row=17, column=0, columnspan=3, sticky='w')
        _Tooltip(self.unc_check,
                 'Repeats the run M times, each on one random realisation of '
                 'the γ, conversion-electron, X-ray and Auger emission '
                 'probabilities (drawn from their stated uncertainties). The '
                 'spread between runs, less the Monte-Carlo noise, is the '
                 'nuclear-data uncertainty. Takes M times as long.')
        self.reps_entry = self._row(
            f, 18, '      Realisations M', 'uncReps', '10', width=8,
            hint='runs of N decays each; total time × M')
        self.display_check = ttk.Checkbutton(
            f, text='Verbose log and decay view (every decay)',
            variable=self._flag('display'))
        self.display_check.grid(row=19, column=0, columnspan=3, sticky='w')
        _Tooltip(self.display_check,
                 'Prints every decay to the Log tab, and draws the first '
                 'ones in the Decay view tab: the sampled decay, the '
                 'particles emitted, the energy they deposit, the light '
                 'left after quenching and the PMT signal. Slows a run '
                 'down; keep N small.')

    def _tab_data(self, nb):
        f = ttk.Frame(nb, padding=12)
        nb.add(f, text='Data & plots')
        self.tab_data = f

        ttk.Label(f, text='Nuclide').grid(row=0, column=0, sticky='w')
        self.data_nuclide = ttk.Combobox(
            f, textvariable=self._var('data_rad', 'Fe-55'), width=16,
            values=())
        self.data_nuclide.grid(row=0, column=1, sticky='w', pady=2)
        ttk.Label(f, text='decay data and β spectrum',
                  style='Hint.TLabel').grid(row=0, column=2, sticky='w',
                                            padx=(8, 0))

        self._section(f, 1, 'β spectrum')
        ttk.Label(f, text='Decay mode').grid(row=2, column=0, sticky='w',
                                             pady=2)
        self.beta_mode_box = ttk.Combobox(
            f, textvariable=self._var('betaMode', 'beta-'), width=10,
            values=('beta-', 'beta+'), state='readonly')
        self.beta_mode_box.grid(row=2, column=1, sticky='w', pady=2)
        self.beta_level_entry = self._row(
            f, 3, 'Transition level', 'betaLevel', 'tot',
            hint="daughter level, or 'tot' for a pure β emitter")

        self._section(f, 4, 'Stopping power and quenching')
        ttk.Label(f, text='Particle').grid(row=5, column=0, sticky='w',
                                           pady=2)
        prow = ttk.Frame(f)
        prow.grid(row=5, column=1, columnspan=2, sticky='w')
        self._var('particle', 'electron')
        self.particle_radios = [
            ttk.Radiobutton(prow, text=text, value=value,
                            variable=self.vars['particle'],
                            command=self._on_particle_change)
            for value, text in (('electron', 'electron'),
                                ('alpha', 'α particle'))]
        for rb in self.particle_radios:
            rb.pack(side='left', padx=(0, 12))
        self._particle_shown = 'electron'

        self.erange_entry = self._row(
            f, 6, 'Energy range', 'Erange', E_DEFAULTS['electron'],
            hint='from, to, number of points (keV)')
        self.elog_check = ttk.Checkbutton(
            f, text='log spacing', variable=self._flag('Elog', True))
        self.elog_check.grid(row=7, column=1, sticky='w')
        self.spcompare_check = ttk.Checkbutton(
            f, text='overlay every electron stopping-power model',
            variable=self._flag('spCompare', False))
        self.spcompare_check.grid(row=8, column=0, columnspan=3, sticky='w')
        self.kbq_entry = self._row(
            f, 9, 'Birks constant(s) kB', 'kBq', '1.0e-5', width=24,
            hint='cm keV⁻¹; several, comma-separated, to compare')

        ttk.Label(
            f, style='Hint.TLabel', wraplength=self._px(430),
            justify='left',
            text='Stopping power and quenching use the cocktail and '
                 'stopping-power model set on the Scintillator tab.',
        ).grid(row=10, column=0, columnspan=3, sticky='w', pady=(12, 0))

    #: Where the stopping-power choice acts. Said in the window because the
    #: answer -- not everywhere -- is not what the field suggests.
    SP_NOTE = ('Electrons, low energy. It enters wherever the quenching '
               'integral is evaluated directly: below {where}, and for any '
               'partial deposit. Above it the shipped quenching tables are '
               'interpolated and are not recomputed.')

    def _tab_scintillator(self, nb):
        f = ttk.Frame(nb, padding=12)
        nb.add(f, text='Scintillator')
        ttk.Label(f, text='Cocktail').grid(row=0, column=0, sticky='w')
        self.cocktail_box = ttk.Combobox(
            f, textvariable=self._var('ls_cocktail'), width=26,
            state='disabled', values=())
        self.cocktail_box.grid(row=0, column=1, sticky='w', pady=2)
        self.cocktail_box.bind('<<ComboboxSelected>>',
                               lambda e: self._apply_L_default())
        self._row(f, 1, 'Aqueous fraction fAq', 'fAq', '0.0',
                  hint='0 for an organic source')
        ttk.Label(f, text='Aqueous solvent').grid(row=2, column=0, sticky='w')
        ttk.Combobox(f, textvariable=self._var('solvantType', 'Water'),
                     values=('Water', 'HCl', 'HNO3', 'NaOH'),
                     width=14).grid(row=2, column=1, sticky='w', pady=2)
        self._row(f, 3, 'Solvent concentration', 'solvantConc', '0.0',
                  hint='mol L⁻¹')
        self._section(f, 4, 'Stopping power and quenching')
        ttk.Label(f, text='Stopping-power model').grid(row=5, column=0,
                                                       sticky='w')
        # Editable until the library says which models it has; read-only
        # afterwards, so a typo cannot reach config.toml.
        self.sp_box = ttk.Combobox(
            f, textvariable=self._var('sp_model', 'tan_xia'), width=20,
            values=())
        self.sp_box.grid(row=5, column=1, sticky='w', pady=2)
        self.sp_note = ttk.Label(f, wraplength=self._px(430),
                                 style='Hint.TLabel', justify='left',
                                 text=self.SP_NOTE.format(where='Einterp_e'))
        self.sp_note.grid(row=6, column=0, columnspan=3, sticky='w',
                          pady=(0, 6))
        self._row(f, 7, 'Chou parameter', 'chou_param', '0',
                  hint='cm² MeV⁻²')
        self._row(f, 8, 'Electron quench bins', 'nE_electron', '1000')
        self._row(f, 9, 'Alpha quench bins', 'nE_alpha', '1000')

    def _tab_micelles(self, nb):
        f = ttk.Frame(nb, padding=12)
        nb.add(f, text='Micelles')
        ttk.Label(
            f, wraplength=self._px(430), style='Hint.TLabel', justify='left',
            text='With an aqueous fraction the source sits in water droplets, '
                 'and an electron slow enough to stop inside one never '
                 'reaches the solvent. Left on, the correction switches '
                 'itself on for the nuclides measured to need it.',
        ).grid(row=0, column=0, columnspan=3, sticky='w', pady=(0, 8))
        ttk.Checkbutton(f, text='Self-activate for nuclides that need it '
                                '(micCorrAuto)',
                        variable=self._flag('micCorrAuto', True)).grid(
            row=1, column=0, columnspan=3, sticky='w')
        ttk.Checkbutton(f, text='Force the correction on for every nuclide '
                                '(micCorr)',
                        variable=self._flag('micCorr')).grid(
            row=2, column=0, columnspan=3, sticky='w')
        self._row(f, 3, 'Bias threshold', 'micCorrBiasThreshold', '0.01',
                  hint='% of the light; the measured-bias gate')
        self._row(f, 4, 'P(EC) threshold', 'micCorrECThreshold', '0.01',
                  hint='fallback, for an unscreened nuclide')
        self._row(f, 5, 'Micelle diameter', 'diam_micelle', '4',
                  hint='nm — a diameter, not a radius')
        self._row(f, 6, 'Size dispersion σ', 'sigma_micelle', '0',
                  hint='nm; 0 is monodisperse')

    def _tab_delayed(self, nb):
        f = ttk.Frame(nb, padding=12)
        nb.add(f, text='Delayed fluorescence')
        ttk.Label(
            f, wraplength=self._px(430), style='Hint.TLabel', justify='left',
            text='Prompt plus delayed (triplet-triplet annihilation) light, '
                 'which makes the yield depend on the coincidence resolving '
                 'time. This is a different model, not a refinement: turning '
                 'it on changes every efficiency reported. Off by default.',
        ).grid(row=0, column=0, columnspan=3, sticky='w', pady=(0, 8))
        ttk.Checkbutton(f, text='Enable delayed fluorescence',
                        variable=self._flag('delayedFluorescence')).grid(
            row=1, column=0, columnspan=3, sticky='w')
        self._row(f, 2, 'Prompt τₚ', 'tau_prompt', '3.0', hint='ns')
        self._row(f, 3, 'Delayed τ_d', 'tau_delayed', '85.0',
                  hint='ns — reset by the cocktail choice')
        self._row(f, 4, 'Yield η', 'eta_delayed', '0.03',
                  hint='reset by the cocktail choice')

    def _tab_counter(self, nb):
        f = ttk.Frame(nb, padding=12)
        nb.add(f, text='Counter')
        self._row(f, 0, 'Coincidence resolving time', 'tau', '50', hint='ns')
        self._row(f, 1, 'Extended dead time', 'extDT', '30',
                  hint='µs')
        self._row(f, 2, 'Measurement time', 'measTime', '60', hint='min')
        self._row(f, 3, 'PMT quantum efficiencies', 'effQuantum',
                  '0.25, 0.25, 0.25', width=20, hint='A, B, C')
        ttk.Checkbutton(f, text='Full optical Monte-Carlo transport '
                                '(slower)',
                        variable=self._flag('opticalTransport')).grid(
            row=4, column=0, columnspan=3, sticky='w', pady=(8, 0))

    # ------------------------------------------------------------ results
    def _tab_table(self, nb):
        f = ttk.Frame(nb, padding=6)
        nb.add(f, text='Table')
        f.columnconfigure(0, weight=1)
        f.rowconfigure(1, weight=1)
        self.caption = ttk.Label(f, text='No results yet.',
                                 style='Hint.TLabel')
        self.caption.grid(row=0, column=0, sticky='w', pady=(0, 4))
        self.tree = ttk.Treeview(f, show='headings', height=12)
        self.tree.tag_configure('odd', background=STRIPE)
        self._set_columns(self.QUANTITY_COLUMNS)
        self.tree.grid(row=1, column=0, sticky='nsew')
        ttk.Scrollbar(f, orient='vertical', command=self.tree.yview).grid(
            row=1, column=1, sticky='ns')
        xs = ttk.Scrollbar(f, orient='horizontal', command=self.tree.xview)
        xs.grid(row=2, column=0, sticky='ew')
        self.tree.configure(xscrollcommand=xs.set)
        self.tab_table = f

    def _tab_plot(self, nb):
        f = ttk.Frame(nb)
        nb.add(f, text='Plot')
        self.plot_frame = f
        self.plot_placeholder = ttk.Label(
            f, style='Hint.TLabel', justify='center',
            text='Plots appear here.\n\nRun the efficiency curve, or any '
                 'item under Data & plots.')
        self.plot_placeholder.place(relx=0.5, rely=0.5, anchor='center')
        self.tab_plot = f

    # ------------------------------------------------------- decay view
    #: Colours of the energy stages, shared by the stage strip and the bars.
    C_EMIT, C_DEP, C_LIGHT = '#9AA9BD', CYAN, '#E07A1F'

    def _tab_decay(self, nb):
        """One simulated decay at a time, drawn as the chain the model
        follows: decay -> emission -> deposition -> quenching -> signal."""
        f = ttk.Frame(nb, padding=6)
        nb.add(f, text='Decay view')
        f.columnconfigure(0, weight=1)
        f.rowconfigure(1, weight=1)

        bar = ttk.Frame(f)
        bar.grid(row=0, column=0, columnspan=2, sticky='ew', pady=(0, 4))
        self.trace_prev = ttk.Button(bar, text='◀', width=3,
                                     command=lambda: self._step_trace(-1))
        self.trace_prev.pack(side='left')
        ttk.Label(bar, text=' Decay ').pack(side='left')
        self.trace_pos = tk.StringVar(value='1')
        self.trace_spin = ttk.Spinbox(bar, from_=1, to=1, width=6,
                                      textvariable=self.trace_pos,
                                      command=self._on_trace_pos)
        self.trace_spin.pack(side='left')
        self.trace_spin.bind('<Return>', self._on_trace_pos)
        self.trace_spin.bind('<FocusOut>', self._on_trace_pos)
        self.trace_count = ttk.Label(bar, text=' of 0 ')
        self.trace_count.pack(side='left')
        self.trace_next = ttk.Button(bar, text='▶', width=3,
                                     command=lambda: self._step_trace(1))
        self.trace_next.pack(side='left')
        self.trace_info = ttk.Label(bar, text='', style='Hint.TLabel')
        self.trace_info.pack(side='left', padx=(12, 0))
        self.trace_save = ttk.Button(bar, text='Save decays…',
                                     command=self._save_trace)
        self.trace_save.pack(side='right')
        _Tooltip(self.trace_save, 'Every recorded decay as JSON: steps, '
                                  'particles with emitted, deposited and '
                                  'quenched energies, and the detection '
                                  'probabilities.')

        c = tk.Canvas(f, background=WHITE, highlightthickness=0)
        c.grid(row=1, column=0, sticky='nsew')
        ys = ttk.Scrollbar(f, orient='vertical', command=c.yview)
        ys.grid(row=1, column=1, sticky='ns')
        xs = ttk.Scrollbar(f, orient='horizontal', command=c.xview)
        xs.grid(row=2, column=0, sticky='ew')
        c.configure(yscrollcommand=ys.set, xscrollcommand=xs.set)
        c.bind('<MouseWheel>', lambda e: c.yview_scroll(
            int(-e.delta / 120) or (-1 if e.delta > 0 else 1), 'units'))
        c.bind('<Button-4>', lambda e: c.yview_scroll(-1, 'units'))
        c.bind('<Button-5>', lambda e: c.yview_scroll(1, 'units'))
        c.bind('<Button-1>', lambda e: c.focus_set())
        c.bind('<Left>', lambda e: self._step_trace(-1))
        c.bind('<Right>', lambda e: self._step_trace(1))
        c.bind('<Configure>', self._on_decay_resize)
        self.trace_canvas = c
        self.tab_decay = f
        self._show_trace([], None)

    def _show_trace(self, trace, job):
        """Take a new set of decay records and draw the first."""
        self._trace = list(trace or [])
        self._trace_job = job
        self._trace_i = 0
        n = len(self._trace)
        self.trace_spin.configure(to=max(1, n))
        self.trace_pos.set('1')
        self.trace_count.configure(text=f' of {n} ')
        state = 'normal' if n else 'disabled'
        for w in (self.trace_prev, self.trace_next, self.trace_spin,
                  self.trace_save):
            w.configure(state=state)
        self._draw_decay()
        self.trace_canvas.yview_moveto(0)
        self.trace_canvas.xview_moveto(0)
        if n:
            self.results.select(self.tab_decay)

    def _step_trace(self, delta):
        if self._trace:
            self._goto_trace(self._trace_i + delta)

    def _on_trace_pos(self, _event=None):
        try:
            self._goto_trace(int(float(self.trace_pos.get())) - 1)
        except ValueError:
            self.trace_pos.set(str(self._trace_i + 1))

    def _goto_trace(self, i):
        if not self._trace:
            return
        self._trace_i = min(len(self._trace) - 1, max(0, i))
        self.trace_pos.set(str(self._trace_i + 1))
        self._draw_decay()

    def _draw_decay(self):
        """Draw the current decay record, laid out for the canvas width."""
        c = self.trace_canvas
        c.delete('all')
        px = self._px
        font = tkfont.nametofont('TkDefaultFont')
        width = max(c.winfo_width(), px(600))
        self._decay_drawn_width = c.winfo_width()
        x0, right = px(14), width - px(14)
        if not self._trace:
            c.create_text(
                px(20), px(20), anchor='nw', fill=MUTED, font=font,
                width=min(px(620), right - px(20)),
                text='Each simulated decay is drawn here, from the decay '
                     'itself through energy deposition and quenching to '
                     'the PMT signal.\n\nTick "Verbose log and decay view" '
                     'on the Source tab, keep N small, and run Efficiency, '
                     'Fit L (stochastic) or the Efficiency curve. The '
                     'analytical model does not simulate decays one by '
                     'one, so it has nothing to show here.')
            self.trace_info.configure(text='')
            c.configure(scrollregion=c.bbox('all'))
            return

        rec = self._trace[self._trace_i]
        job = self._trace_job or {}
        L = rec.get('L')
        L_text = (', '.join(f'{v:.4g}' for v in L)
                  if isinstance(L, (list, tuple)) else f'{L:.4g}')
        where = {'fit': ' (fitted)', 'curve': ' (first point of the curve)'}
        self.trace_info.configure(
            text=f'{rec.get("nuclide")}  ·  L = {L_text} keV⁻¹'
                 f'{where.get(job.get("task"), "")}  ·  '
                 f'{len(self._trace)} decay(s) recorded'
                 + (f' of N = {job["N"]}' if job.get('N')
                    and int(job['N']) > len(self._trace) else ''))

        def text(x, y, s, wrap=None, **kw):
            """Draw text; return the y just below it."""
            kw.setdefault('font', font)
            kw.setdefault('fill', INK)
            item = c.create_text(x, y, anchor='nw', text=s,
                                 width=wrap or 0, **kw)
            return c.bbox(item)[3]

        y = px(10)
        # ---- the five stages, as a strip of boxes with arrows between
        kB, V = job.get('kB'), job.get('V')
        stages = (
            ('1  Decay', 'PenNuc data', NAVY),
            ('2  Emission', 'relaxation, β spectra', self.C_EMIT),
            ('3  Deposition', f'in {V:g} mL' if V else 'in the vial',
             self.C_DEP),
            ('4  Quenching', f'Birks, kB = {kB:g}' if kB else 'Birks',
             self.C_LIGHT),
            ('5  Signal', 'PMTs A, B, C', NAVY),
        )
        gap = px(18)
        w = (right - x0 - 4 * gap) / 5
        # Text first, so the boxes can be as tall as the longest wraps.
        bottom = y + px(40)
        for i, (title, sub, colour) in enumerate(stages):
            xa = x0 + i * (w + gap)
            c.create_text(xa + px(11), y + px(5), anchor='nw', text=title,
                          font=self.f_bold, fill=NAVY)
            item = c.create_text(xa + px(11), y + px(25), anchor='nw',
                                 text=sub, font=self.f_sub, fill=MUTED,
                                 width=w - px(14))
            bottom = max(bottom, c.bbox(item)[3])
        h = bottom - y + px(6)
        for i, (_title, _sub, colour) in enumerate(stages):
            xa = x0 + i * (w + gap)
            box = c.create_rectangle(xa, y, xa + w, y + h, fill=PANEL,
                                     outline=RULE)
            strip = c.create_rectangle(xa, y, xa + px(5), y + h,
                                       fill=colour, outline=colour)
            c.tag_lower(strip)
            c.tag_lower(box)
            if i:
                c.create_line(xa - gap + px(3), y + h / 2, xa - px(3),
                              y + h / 2, arrow='last', fill=MUTED,
                              width=px(2))
        y += h + px(18)
        wrap = right - x0 - px(28)
        xt = x0 + px(28)

        # ---- 1 decay
        y = self._decay_heading(
            c, x0, y, '1', NAVY,
            f'Decay {rec.get("decay")}:  {rec.get("nuclide")}  →  '
            f'{nuclide_name(rec.get("daughter"))}')
        for step in rec.get('steps') or ():
            y = text(xt, y, '→  ' + step_label(step), wrap) + px(3)
        y += px(10)

        # ---- 2, 3, 4 per coincidence window
        emit_r = xt + px(268)               # right edge of the Emitted column
        dep_x = emit_r + px(16)
        group = max(px(120), (right - dep_x) / 2)
        light_x = dep_x + group
        # The bar gives way to its label, "0.1936 keV  (34 %)".
        bw = max(px(40), min(px(150), group - px(150)))
        for win in rec.get('windows') or ():
            delayed = win.get('label') == 'delayed'
            title = ('Prompt' if not delayed else
                     f'Delayed, t = {format_time(rec.get("t1", 0) * 1e9)} '
                     f'after the decay')
            y = self._decay_heading(c, x0, y, '2', self.C_DEP,
                                    f'{title}:  emitted → deposited → light')
            for x, head, anchor in ((xt, 'Particle', 'nw'),
                                    (emit_r, 'Emitted', 'ne'),
                                    (dep_x, 'Deposited', 'nw'),
                                    (light_x, 'Light after quenching',
                                     'nw')):
                c.create_text(x, y, anchor=anchor, text=head,
                              font=self.f_section, fill=MUTED)
            y += px(18)
            parts = win.get('particles') or ()
            if not parts:
                y = text(xt, y, 'nothing emitted in this window',
                         fill=MUTED) + px(6)
            for p in parts:
                e0, ed, eq = p['E0'], p['E_dep'], p['E_q']
                name = particle_label(p['origin'], p['type'])
                c.create_text(xt, y, anchor='nw', text=name, font=font,
                              fill=INK)
                note = FATE_LABEL.get(p.get('fate'), '')
                if note:
                    c.create_text(xt + font.measure(name) + px(8),
                                  y + px(2), anchor='nw', text=note,
                                  font=self.f_sub, fill=MUTED)
                c.create_text(emit_r, y, anchor='ne', font=font, fill=INK,
                              text=f'{e0:.4g} keV')
                self._energy_bar(c, dep_x, y, bw, ed, e0, self.C_DEP,
                                 f'{ed:.4g} keV')
                share = f'  ({100 * eq / ed:.0f} %)' if ed > 0 else ''
                self._energy_bar(c, light_x, y, bw, eq, e0, self.C_LIGHT,
                                 f'{eq:.4g} keV{share}' if ed > 0 else '—')
                y += px(22)
            c.create_line(xt, y, right, y, fill=RULE)
            y = text(light_x, y + px(4),
                     f'Q = {win.get("Q", 0.0):.4g} keV of light',
                     font=self.f_bold, fill=self.C_LIGHT) + px(12)
        y = text(xt, y, 'Bars: share of the energy each particle was '
                        'emitted with. Deposited light is quenched by the '
                        'Birks law; the percentage is what survives.',
                 wrap, font=self.f_sub, fill=MUTED) + px(14)

        # ---- 5 signal
        y = self._decay_heading(c, x0, y, '5', NAVY,
                                'Signal:  photoelectrons and coincidences')
        windows = rec.get('windows') or []
        counted = counted_windows(rec)
        if len(windows) > len(counted):
            y = text(xt, y, 'The delayed light is not counted: it arrives '
                            'within the extended dead time after the prompt '
                            'pulse, or after the end of the count.',
                     wrap, fill=MUTED) + px(8)
        elif len(counted) > 1:
            y = text(xt, y, 'Two pulses: the delayed light comes after the '
                            'dead time, so it can be counted again. The '
                            'values below add both up and can exceed 1.',
                     wrap, fill=MUTED) + px(8)

        mu = rec.get('mu') or []
        spacing = min(px(170), (right - xt) / 3)
        for win in counted:
            lam = photoelectrons(rec, win.get('Q', 0.0))
            label = ('Delayed pulse' if win.get('label') == 'delayed'
                     else 'Prompt pulse')
            if lam is None:
                y = text(xt, y, f'{label}: optical transport model. The '
                                f'photons are sampled to each PMT, so there '
                                f'is no mean photoelectron number to show.',
                         wrap) + px(8)
                continue
            if isinstance(L, (list, tuple)):
                formula = (f'{label}:  λ_i = L_i · μ_i · Q / 3,  with '
                           f'Q = {win["Q"]:.4g} keV')
            else:
                m = sum(mu) / len(mu) if mu else 1.0
                formula = (f'{label}:  λ = L · μ · Q / 3 = {L:.4g} × '
                           f'{m:.3g} × {win["Q"]:.4g} / 3 = {lam[0]:.3g} '
                           f'photoelectrons per PMT, on average')
            y = text(xt, y, formula, wrap) + px(8)
            r = px(20)
            for k, (name, lm) in enumerate(zip('ABC', lam)):
                p_fire = 1.0 - math.exp(-lm)
                cx = xt + r + k * spacing
                c.create_oval(cx - r, y, cx + r, y + 2 * r,
                              fill=blend(WHITE, CYAN, p_fire), outline=NAVY,
                              width=px(2))
                c.create_text(cx, y + r, text=name, font=self.f_bold,
                              fill=NAVY)
                c.create_text(cx + r + px(6), y + px(4), anchor='nw',
                              font=self.f_sub, fill=INK,
                              text=f'λ = {lm:.3g}\nP(fires) = '
                                   f'{100 * p_fire:.1f} %')
            y += 2 * r + px(14)

        sig = rec.get('signal') or {}
        scale = max(1.0, float(len(counted)))
        rows = (('S', 'Single: at least one PMT'),
                ('D', 'Double (logical sum): two or more'),
                ('T', 'Triple: all three PMTs'),
                ('D2', 'Double, C/N counter'))
        bar_x = xt + max(font.measure(t) for _k, t in rows) + px(14)
        bar_w = max(px(60), min(px(280), right - bar_x - px(70)))
        for key, label in rows:
            v = float(sig.get(key, 0.0))
            c.create_text(xt, y, anchor='nw', font=font, fill=INK,
                          text=label)
            self._energy_bar(c, bar_x, y, bar_w, v, scale, NAVY, f'{v:.4f}')
            y += px(22)
        d = sig.get('D', 0.0)
        if d > 0:
            text(xt, y + px(4),
                 f'T / D for this decay = {sig.get("T", 0.0) / d:.4f}. The '
                 f'efficiencies of the run are these probabilities averaged '
                 f'over all N decays.', wrap, font=self.f_sub, fill=MUTED)
        bbox = c.bbox('all')
        c.configure(scrollregion=(0, 0, bbox[2] + px(10), bbox[3] + px(16)))

    def _on_decay_resize(self, event):
        """Lay the drawing out again when the pane changes width."""
        if abs(event.width - getattr(self, '_decay_drawn_width', 0)) > 4:
            self._draw_decay()

    def _decay_heading(self, c, x, y, number, colour, text):
        """A numbered section heading in the decay view; returns the next y."""
        px = self._px
        c.create_rectangle(x, y + px(2), x + px(18), y + px(18),
                           fill=colour, outline=colour)
        c.create_text(x + px(28), y, anchor='nw', text=text,
                      font=self.f_bold, fill=NAVY)
        c.create_text(x + px(9), y + px(10), text=number[0],
                      font=self.f_section, fill=WHITE)
        return y + px(24)

    def _energy_bar(self, c, x, y, width, value, full, colour, text):
        """A horizontal bar for value/full with its label to the right."""
        px = self._px
        top, bottom = y + px(3), y + px(15)
        c.create_rectangle(x, top, x + width, bottom, fill=STRIPE,
                           outline=RULE)
        frac = value / full if full and full > 0 else 0.0
        if frac > 0:
            c.create_rectangle(x, top, x + max(px(1), width *
                                                min(1.0, frac)), bottom,
                               fill=colour, outline=colour)
        c.create_text(x + width + px(6), y, anchor='nw', text=text,
                      font=self.f_sub, fill=INK)

    def _save_trace(self):
        if not self._trace:
            return
        job = self._trace_job or {}
        stem = '_'.join(str(p) for p in ('decays', job.get('rad', ''))
                        if p).replace(', ', '-')
        path = filedialog.asksaveasfilename(
            title='Save the recorded decays', defaultextension='.json',
            initialfile=f'{stem}.json',
            filetypes=[('JSON', '*.json'), ('All files', '*.*')])
        if not path:
            return
        payload = {'tdcrpy': self.version,
                   'run': self._run_summary(job) if job else '',
                   'saved': datetime.datetime.now().isoformat(
                       timespec='seconds'),
                   'decays': self._trace}
        try:
            with open(path, 'w', encoding='utf-8') as fh:
                json.dump(payload, fh, indent=1, ensure_ascii=False)
        except OSError as exc:
            messagebox.showerror('TDCRPy', f'Could not save:\n{exc}')
            return
        self._append_log(f'{len(self._trace)} decay(s) written to {path}')

    def _tab_sources(self, nb):
        """Where the data behind a decay-data or spectrum view come from."""
        f = ttk.Frame(nb, padding=6)
        nb.add(f, text='Sources')
        f.columnconfigure(0, weight=1)
        f.rowconfigure(0, weight=1)
        self.sources_text = tk.Text(
            f, wrap='word', relief='flat', state='disabled', background=WHITE,
            padx=10, pady=8, font=tkfont.nametofont('TkDefaultFont'))
        self.sources_text.grid(row=0, column=0, sticky='nsew')
        ttk.Scrollbar(f, orient='vertical',
                      command=self.sources_text.yview).grid(row=0, column=1,
                                                            sticky='ns')
        t = self.sources_text
        t.tag_configure('title', font=self.f_bold, foreground=NAVY,
                        spacing1=self._px(8), spacing3=self._px(3))
        t.tag_configure('label', foreground=MUTED)
        t.tag_configure('link', foreground=CYAN, underline=True)
        t.tag_bind('link', '<Enter>', lambda e: t.configure(cursor='hand2'))
        t.tag_bind('link', '<Leave>', lambda e: t.configure(cursor=''))
        self._show_info(None)
        self.tab_sources = f

    def _show_info(self, info):
        """Fill the Sources tab from a result's provenance sections."""
        t = self.sources_text
        t.configure(state='normal')
        t.delete('1.0', 'end')
        if not info:
            t.insert('end', 'Decay data and β spectra list here where their '
                            'data come from: the DDEP evaluation, the '
                            'atomic-data evaluator and the BetaShape version '
                            'of every spectrum.', 'label')
        for i, section in enumerate(info or ()):
            t.insert('end', section['title'] + '\n', 'title')
            for label, text, url in section['items']:
                if label:
                    t.insert('end', f'{label}:  ', 'label')
                else:
                    t.insert('end', '      ')
                if url:
                    tag = f'url{i}_{t.index("end")}'
                    t.insert('end', text, ('link', tag))
                    t.tag_bind(tag, '<Button-1>',
                               lambda e, u=url: self._open(u))
                else:
                    t.insert('end', text)
                t.insert('end', '\n')
        t.configure(state='disabled')

    def _tab_log(self, nb):
        f = ttk.Frame(nb, padding=6)
        nb.add(f, text='Log')
        f.columnconfigure(0, weight=1)
        f.rowconfigure(0, weight=1)
        self.log = tk.Text(f, height=10, wrap='none', relief='flat',
                           font=('Consolas', 9), state='disabled',
                           background=STRIPE, padx=6, pady=4)
        self.log.grid(row=0, column=0, sticky='nsew')
        ttk.Scrollbar(f, orient='vertical', command=self.log.yview).grid(
            row=0, column=1, sticky='ns')
        self.tab_log = f

    #: Results table for a single calculation: one quantity per row.
    QUANTITY_COLUMNS = (
        ('quantity', 'Quantity', 200, 'w'),
        ('value', 'Value', 110, 'e'),
        ('u', 'u', 90, 'e'),
        ('rel', 'u / value', 80, 'e'),
    )

    def _column_spec(self, columns, rows):
        """Columns sized to their heading and content; text columns read
        left-aligned, numbers right-aligned."""
        body = tkfont.nametofont('TkDefaultFont')
        sample = rows[:300]
        spec = []
        for i, (key, title) in enumerate(columns):
            digits = 3 if is_uncertainty_column(key) else 5
            cells = [format_cell(r[i], digits) for r in sample if i < len(r)]
            text = any(isinstance(r[i], str) and r[i] and not numberish(r[i])
                       for r in sample if i < len(r))
            width = max([self.f_bold.measure(title)] +
                        [body.measure(c) for c in cells]) + 18
            # Already in screen pixels: undo the scaling _set_columns applies.
            spec.append((key, title, width / self.s, 'w' if text else 'e'))
        return spec

    def _set_columns(self, spec, stretch=True):
        """Rebuild the results table's columns from (key, title, width,
        anchor) tuples. Without *stretch* the columns keep their width and
        a wide table scrolls sideways instead of squeezing its last ones."""
        self.tree.configure(columns=[c[0] for c in spec])
        for key, title, width, anchor in spec:
            self.tree.heading(key, text=title)
            self.tree.column(key, width=self._px(width), anchor=anchor,
                             stretch=stretch, minwidth=20)

    # ------------------------------------------------- action and status bar
    def _build_actions(self):
        bar = ttk.Frame(self, padding=(self._px(10), self._px(8)))
        bar.grid(row=2, column=1, sticky='ew')
        bar.columnconfigure(5, weight=1)

        self.run_button = tk.Button(
            bar, text='▶  Run', command=self._run, font=self.f_bold,
            fg=WHITE, bg=NAVY, activeforeground=WHITE, activebackground=CYAN,
            disabledforeground='#9FB0C6', relief='flat', bd=0,
            padx=self._px(18), pady=self._px(4), cursor='hand2')
        self.run_button.grid(row=0, column=0)
        _Tooltip(self.run_button, 'Run the selected item (F5)')
        self.cancel_button = ttk.Button(bar, text='Cancel', state='disabled',
                                        command=self._cancel)
        self.cancel_button.grid(row=0, column=1, padx=(6, 0))
        ttk.Separator(bar, orient='vertical').grid(row=0, column=2,
                                                   sticky='ns', padx=10)
        self.export_button = ttk.Button(bar, text='Export CSV…',
                                        command=self._save, state='disabled')
        self.export_button.grid(row=0, column=3)
        self.saveplot_button = ttk.Button(bar, text='Save plot…',
                                          command=self._save_plot,
                                          state='disabled')
        self.saveplot_button.grid(row=0, column=4, padx=(6, 0))

        self.progress = ttk.Progressbar(bar, mode='determinate',
                                        length=self._px(180))
        self.progress.grid(row=0, column=6, sticky='e', padx=8)
        self.status = ttk.Label(bar, text='ready', width=28, anchor='e')
        self.status.grid(row=0, column=7, sticky='e')

    def _build_statusbar(self):
        bar = tk.Frame(self, bg=PANEL)
        bar.grid(row=3, column=0, columnspan=2, sticky='ew')
        tk.Frame(bar, bg=RULE, height=1).pack(side='top', fill='x')
        tk.Label(bar, bg=PANEL, fg=MUTED, font=self.f_sub, anchor='w',
                 text='Running writes these settings to the package '
                      'config.toml, shared with scripts and notebooks — '
                      'Tools ▸ Undo my configuration changes puts it back.'
                 ).pack(side='left', padx=self._px(12), pady=self._px(4))
        for text, action in (('MIT licence', self._show_licence),
                             ('BIPM', lambda: self._open(BIPM_URL)),
                             ('Documentation',
                              lambda: self._open(DOCS_URL))):
            link = tk.Label(bar, text=text, bg=PANEL, fg=NAVY,
                            font=self.f_sub, cursor='hand2')
            link.pack(side='right', padx=(0, self._px(12)))
            link.bind('<Button-1>', lambda e, a=action: a())

    # -------------------------------------------------------------- config
    def _config_paths(self):
        d = package_dir()
        if d is None:
            return None, None
        return (os.path.join(d, 'config.toml'),
                os.path.join(d, 'configDefault.toml'))

    def _snapshot_config(self):
        """Keep the configuration as it was before this window touched it."""
        live, _ = self._config_paths()
        if not live or not os.path.exists(live):
            return
        try:
            with open(live, encoding='utf-8') as fh:
                self._config_snapshot = fh.read()
        except OSError:
            self._config_snapshot = None

    def _restore_snapshot(self):
        live, _ = self._config_paths()
        if not live or self._config_snapshot is None:
            messagebox.showwarning(
                'TDCRPy', 'No snapshot was taken, so there is nothing to '
                          'put back.')
            return
        if not messagebox.askyesno(
                'TDCRPy', 'Put config.toml back as it was when this window '
                          'opened?'):
            return
        try:
            with open(live, 'w', encoding='utf-8') as fh:
                fh.write(self._config_snapshot)
        except OSError as exc:
            messagebox.showerror('TDCRPy', f'Could not write config.toml:'
                                           f'\n{exc}')
            return
        self._append_log('config.toml restored to the session snapshot')
        self._enquire()

    def _restore_defaults(self):
        live, default = self._config_paths()
        if not live or not default or not os.path.exists(default):
            messagebox.showwarning('TDCRPy',
                                   'configDefault.toml was not found.')
            return
        if not messagebox.askyesno(
                'TDCRPy', 'Replace config.toml with the shipped defaults? '
                          'Any settings you had saved there will be lost.'):
            return
        try:
            shutil.copyfile(default, live)
        except OSError as exc:
            messagebox.showerror('TDCRPy', f'Could not copy the defaults:'
                                           f'\n{exc}')
            return
        self._append_log('config.toml replaced with the shipped defaults')
        self._enquire()

    def _open_config_folder(self):
        d = package_dir()
        if not d:
            return
        try:
            if os.name == 'nt':
                os.startfile(d)
            else:
                webbrowser.open('file://' + d)
        except OSError as exc:
            messagebox.showerror('TDCRPy', f'Could not open {d}:\n{exc}')

    def collect_config(self):
        """Gather the configuration fields as the worker expects them."""
        name = ('ls_cocktail', 'fAq', 'solvantType', 'solvantConc',
                'sp_model', 'chou_param', 'nE_electron', 'nE_alpha',
                'micCorr', 'micCorrAuto', 'micCorrBiasThreshold',
                'micCorrECThreshold', 'diam_micelle', 'sigma_micelle',
                'delayedFluorescence', 'tau_prompt', 'tau_delayed',
                'eta_delayed', 'tau', 'extDT', 'measTime', 'effQuantum',
                'opticalTransport')
        cfg = {}
        for key in name:
            var = self.vars.get(key)
            if var is None:
                continue
            value = var.get()
            if isinstance(value, str) and not value.strip():
                continue          # leave whatever is in config.toml alone
            cfg[key] = value
        return cfg

    # ----------------------------------------------------------- behaviour
    def _on_task_change(self):
        """Show only the fields the chosen item uses, and bring its inputs
        forward."""
        task = self.vars['task'].get()
        self.heading.configure(text=task_label(task))
        self.subheading.configure(text=TASK_HELP.get(task, ''))

        fitting = task in ('fit', 'fit_analytical')
        curve = task == 'curve'
        beta = task == 'beta_spectrum'
        decay = task == 'decay_data'
        stopping = task == 'stopping_power'
        quench = task == 'quenching'
        data = task in dict(DATA_TASKS)
        electron = self.vars['particle'].get() == 'electron'
        # Repeated runs can only measure the data uncertainty where there
        # is sampling to repeat: not the analytical model, not the curve.
        unc_ok = task == 'fit' or (task == 'efficiency'
                                   and bool(self.vars['Smodel'].get()))
        widgets = [
            (self.l_entry, task == 'efficiency'),
            (self.td_entry, fitting),
            (self.guess_entry, fitting),
            (self.bounds_entry, fitting),
            (self.range_entry, curve),
            (self.log_check, curve),
            (self.nuclide_box, not data),
            (self.frac_entry, not data),
            (self.add_button, not data),
            (self.kb_entry, not data),
            (self.v_entry, not data),
            (self.n_entry, task != 'fit_analytical' and not data),
            (self.smodel_check, not fitting and not data),
            (self.unc_check, unc_ok),
            (self.reps_entry, unc_ok and bool(self.vars['uncData'].get())),
            (self.display_check, not data),
        ] + [(b, not data) for b in self.source_buttons]
        for widget, shown in widgets + [
            (self.data_nuclide, beta or decay),
            (self.beta_level_entry, beta),
            (self.erange_entry, stopping or quench),
            (self.elog_check, stopping or quench),
            (self.spcompare_check, stopping and electron),
            (self.kbq_entry, quench),
        ]:
            widget.configure(state='normal' if shown else 'disabled')
        for rb in self.particle_radios:
            rb.configure(state='normal' if stopping or quench else 'disabled')
        # A combobox, not an entry: 'normal' would let free text through.
        self.beta_mode_box.configure(state='readonly' if beta else 'disabled')
        self.inputs.select(self.tab_data if data else self.tab_source)

    def _on_particle_change(self):
        """Swap the energy range to the new particle's default -- unless the
        user typed one of their own."""
        new = self.vars['particle'].get()
        old = self._particle_shown
        if self.vars['Erange'].get().strip() == E_DEFAULTS.get(old):
            self.vars['Erange'].set(E_DEFAULTS[new])
        self._particle_shown = new
        self._on_task_change()

    # ------------------------------------------------------------- source
    def _known_nuclides(self):
        return list(self.nuclide_box.cget('values') or ())

    def _add_nuclide(self):
        """Add the picked nuclide with its fraction, or update its fraction
        if it is already in the source."""
        try:
            self.source = add_to_source(
                self.source, self.vars['pick_rad'].get(),
                self.vars['pick_frac'].get(), self._known_nuclides())
        except InputError as exc:
            messagebox.showerror('Check the source', str(exc))
            return
        self._refresh_source()
        self.nuclide_box.focus_set()
        self.nuclide_box.selection_range(0, 'end')

    def _remove_nuclide(self):
        sel = self.source_tree.selection()
        if not sel:
            return
        rad = self.source_tree.item(sel[0], 'values')[0]
        self.source = [(r, f) for r, f in self.source if r != rad]
        self._refresh_source()

    def _clear_source(self):
        self.source = []
        self._refresh_source()

    def _normalise_source(self):
        self.source = normalise_source(self.source)
        self._refresh_source()

    def _on_source_select(self, _event=None):
        """Load the selected row into the picker, ready to edit."""
        sel = self.source_tree.selection()
        if sel:
            rad, frac = self.source_tree.item(sel[0], 'values')
            self.vars['pick_rad'].set(rad)
            self.vars['pick_frac'].set(frac)

    def _refresh_source(self):
        """Redraw the source list and write it into 'rad' and 'pmf'."""
        tree = self.source_tree
        tree.delete(*tree.get_children())
        for rad, frac in self.source:
            tree.insert('', 'end', values=(rad, f'{frac:g}'))
        rad, pmf = source_fields(self.source)
        self.vars['rad'].set(rad)
        self.vars['pmf'].set(pmf)
        total = sum(f for _r, f in self.source)
        if not self.source:
            text = 'Add at least one nuclide.'
        elif len(self.source) == 1:
            text = 'Single nuclide: its fraction is 1 whatever is entered.'
        elif abs(total - 1.0) > 0.02:
            text = (f'Fractions sum to {total:g}, not 1 — press Normalise '
                    f'if they are activities or ratios.')
        else:
            text = f'Fractions sum to {total:g}.'
        self.source_sum.configure(text=text)
        if self.source and 'data_rad' in self.vars \
                and not self.vars['data_rad'].get().strip():
            self.vars['data_rad'].set(self.source[0][0])

    def _open(self, url):
        try:
            webbrowser.open(url)
        except webbrowser.Error as exc:
            messagebox.showerror('TDCRPy', f'Could not open {url}:\n{exc}')

    #: Lines kept in the log pane. "Verbose" prints every decay, so a run of
    #: 10 000 can emit hundreds of thousands of lines; an unbounded Text
    #: widget would consume the memory and then the frame rate.
    LOG_LIMIT = 4000
    LOG_TRIM = 1000

    def _append_log(self, text):
        self.log.configure(state='normal')
        self.log.insert('end', text + '\n')
        lines = int(self.log.index('end-1c').split('.')[0])
        if lines > self.LOG_LIMIT:
            self.log.delete('1.0', f'{self.LOG_TRIM}.0')
        self.log.see('end')
        self.log.configure(state='disabled')

    def _set_busy(self, busy, status=None):
        self.run_button.configure(state='disabled' if busy else 'normal')
        self.cancel_button.configure(state='normal' if busy else 'disabled')
        if status:
            self.status.configure(text=status)
        if not busy:
            self.progress.configure(value=0)

    # ------------------------------------------------------------ the runs
    def _enquire(self):
        """Ask the library what it offers, and fill the widgets with it."""
        if self._worker is not None:
            return
        self._set_busy(True, 'reading configuration…')
        self._worker = Worker(['--describe'], None, self._queue)

    def _run_summary(self, job):
        """One log line describing the job -- built field by field, since
        the data views carry none of 'rad', 'kB', 'V' or 'N'."""
        parts = [f'run: {job["task"]}']
        if 'rad' in job:
            parts.append(job['rad'])
        if 'particle' in job:
            parts.append(job['particle'])
        if 'kB' in job:
            parts.append(f'kB={job["kB"]:g}')
        if 'kB_list' in job:
            parts.append('kB=' + ', '.join(f'{k:g}' for k in job['kB_list']))
        if 'V' in job:
            parts.append(f'V={job["V"]:g}')
        if 'N' in job:
            parts.append(f'N={job["N"]}')
        if 'L_grid' in job:
            g = job['L_grid']
            parts.append(f'L={g[0]:g}…{g[-1]:g} ({len(g)} points)')
        if 'E_grid' in job:
            g = job['E_grid']
            parts.append(f'E={g[0]:g}…{g[-1]:g} keV ({len(g)} points)')
        if 'mode' in job:
            parts.append(f'{job["mode"]}  level={job.get("level")}')
        if job.get('uncData'):
            parts.append(f'nuclear data: {job.get("uncReps")} realisations')
        return '  '.join(parts)

    def _run(self):
        if self._worker is not None:
            if self._last_job is None and self._describe is None:
                self.status.configure(text='still loading — one moment')
            else:
                messagebox.showinfo(
                    'TDCRPy', 'A calculation is already running. Wait for '
                              'it, or press Cancel.')
            return
        form = {k: v.get() for k, v in self.vars.items()}
        form['config'] = self.collect_config()
        try:
            job = build_job(form)
        except InputError as exc:
            messagebox.showerror('Check the inputs', str(exc))
            return

        self._append_log('-' * 68)
        self._append_log(self._run_summary(job))
        self._last_job = job
        self.progress.configure(mode='indeterminate')
        self.progress.start(40)
        self._set_busy(True, 'running…')
        self._worker = Worker(['--run'], json.dumps(job), self._queue)

    def _cancel(self):
        if self._worker is not None:
            # Remember the intent. The exit code cannot be trusted to say
            # this: killing a process yields 1 on Windows and -9 on POSIX,
            # and 1 is also what an ordinary failure returns, so a cancelled
            # run would otherwise be reported to the user as "failed".
            self._cancelled = True
            self._worker.cancel()
            self._append_log('cancelled')

    def _drain(self):
        """Move worker messages onto the widgets, in the main thread.

        Capped per tick: a verbose run produces lines faster than they can be
        rendered, and draining the whole queue in one go would starve the
        event loop and make Cancel unresponsive -- exactly when it is wanted.
        """
        try:
            for _ in range(200):
                kind, payload = self._queue.get_nowait()
                if kind == 'log':
                    self._handle_log(payload)
                elif kind == 'msg':
                    self._handle_message(payload)
                elif kind == 'exit':
                    self._handle_exit(payload)
                elif kind == 'update':
                    self._show_update(payload)
        except queue.Empty:
            pass
        self.after(80, self._drain)

    def _handle_log(self, text):
        pct = _TQDM_PCT.search(text)
        if pct and '|' in text:
            if str(self.progress.cget('mode')) == 'indeterminate':
                self.progress.stop()
                self.progress.configure(mode='determinate')
            self.progress.configure(value=min(100, int(pct.group(1))))
            return                      # the bar says it; the log needn't
        if text.strip():
            self._append_log(text)

    def _handle_message(self, msg):
        event = msg.get('event')
        if event == 'stage':
            self.status.configure(text=msg.get('text', '')[:60])
            self._append_log(f'• {msg.get("text", "")}')
        elif event == 'progress':
            if str(self.progress.cget('mode')) == 'indeterminate':
                self.progress.stop()
                self.progress.configure(mode='determinate')
            self.progress.configure(value=min(100.0, float(msg['pct'])))
        elif event == 'describe':
            self._apply_describe(msg)
        elif event == 'result':
            self._result_job = self._last_job
            self._result_info = msg.get('info')
            self._show_info(self._result_info)
            if msg.get('columns'):
                self._show_curve(msg['columns'], msg.get('rows') or [],
                                 msg.get('plot') or {})
            else:
                self._show_result(msg.get('rows') or [])
            self.export_button.configure(state='normal')
            task = (self._result_job or {}).get('task')
            if 'trace' in msg:
                self._show_trace(msg['trace'] or [], self._result_job)
            elif task in dict(CALC_TASKS):
                self._show_trace([], None)      # never a stale decay
        elif event == 'error':
            self._append_log(msg.get('traceback') or msg.get('message', ''))
            self.results.select(self.tab_log)
            messagebox.showerror('TDCRPy failed', msg.get('message',
                                                          'unknown error'))

    def _apply_describe(self, msg):
        self._describe = msg
        cocktails = msg.get('cocktails') or []
        if cocktails:
            self.cocktail_box.configure(values=cocktails, state='readonly')
        nuclides = msg.get('nuclides') or []
        if nuclides:
            # Editable, so a name can be typed; Add checks it is in the list.
            self.nuclide_box.configure(values=nuclides)
            self.data_nuclide.configure(values=nuclides)
        sp_models = msg.get('sp_models') or []
        if sp_models:
            self.sp_box.configure(values=sp_models, state='readonly')
        if msg.get('Einterp_e') is not None:
            self.sp_note.configure(text=self.SP_NOTE.format(
                where=f'Einterp_e = {msg["Einterp_e"]:g} keV'))
        for key, value in (msg.get('config') or {}).items():
            var = self.vars.get(key)
            if var is None:
                continue
            if isinstance(var, tk.BooleanVar):
                var.set(bool(value))
            else:
                var.set('' if value is None else str(value))
        self._apply_L_default()
        self._append_log(f'TDCRPy {msg.get("version", self.version)}: '
                         f'configuration loaded from the installed package '
                         f'({len(nuclides)} nuclides, '
                         f'{len(cocktails)} cocktails)')

    def _apply_L_default(self):
        """Put the cocktail's typical L into the L fields -- unless the
        user has typed a value of their own there."""
        defaults = (self._describe or {}).get('L_defaults') or {}
        value, hint = L_default_for(self.vars['ls_cocktail'].get(), defaults)
        if value is None:                     # a library without the table
            value, hint = '1.0', ''
        for name in ('L', 'L_guess'):
            current = self.vars[name].get().strip()
            if not current or current == self._L_shown.get(name):
                self.vars[name].set(value)
                self._L_shown[name] = value
        self.L_hint.configure(text=hint)

    def _caption(self, nrows):
        job = self._result_job or {}
        what = task_label(job.get('task', ''))
        subject = job.get('rad') or job.get('particle', '')
        return f'{what} — {subject}   ·   {nrows} row(s)'

    def _show_result(self, rows):
        self.tree.delete(*self.tree.get_children())
        self._set_columns(self.QUANTITY_COLUMNS)
        self._rows = [['quantity', 'value', 'u', 'u_rel']]
        for i, (quantity, value, unc) in enumerate(rows):
            cells = (quantity,) + format_value(value, unc)
            self.tree.insert('', 'end', values=cells,
                             tags=('odd',) if i % 2 else ())
            self._rows.append(list(cells))
        self.caption.configure(text=self._caption(len(rows)))
        self.status.configure(text='done')
        self.results.select(self.tab_table)

    def _show_curve(self, columns, rows, plot):
        """Fill the table with the curve or data, and draw it."""
        self.tree.delete(*self.tree.get_children())
        self._set_columns(self._column_spec(columns, rows), stretch=False)
        # Full precision in the saved file; the table is for reading.
        self._rows = [[key for key, _title in columns]]
        digits = [3 if is_uncertainty_column(k) else 5 for k, _t in columns]
        for i, row in enumerate(rows):
            self.tree.insert('', 'end',
                             values=[format_cell(v, d)
                                     for v, d in zip(row, digits)],
                             tags=('odd',) if i % 2 else ())
            self._rows.append(['' if v is None else
                               (v if isinstance(v, str) else repr(v))
                               for v in row])
        self.caption.configure(text=self._caption(len(rows)))
        self.status.configure(text='done')
        if not plot:                        # a table, nothing to draw
            self.results.select(self.tab_table)
            return
        try:
            self._plot([k for k, _t in columns], rows, plot)
        except ImportError as exc:
            self._append_log(f'no plot: matplotlib is not importable ({exc})')
            self.results.select(self.tab_table)
            return
        self.results.select(self.tab_plot)

    # --------------------------------------------------------------- plots
    def _ensure_canvas(self):
        """Create the plot area the first time it is needed.

        matplotlib is imported here, not at the top: the window opens
        without paying for it, and many sessions never draw.
        """
        if self._canvas is not None:
            return
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                       NavigationToolbar2Tk)
        self.plot_placeholder.place_forget()
        self._figure = Figure(figsize=(6.4, 4.6), dpi=100)
        self._canvas = FigureCanvasTkAgg(self._figure, master=self.plot_frame)
        try:
            toolbar = NavigationToolbar2Tk(self._canvas, self.plot_frame,
                                           pack_toolbar=False)
            toolbar.pack(side='bottom', fill='x')
        except TypeError:            # matplotlib < 3.3
            toolbar = NavigationToolbar2Tk(self._canvas, self.plot_frame)
        toolbar.update()
        self._canvas.get_tk_widget().pack(side='top', fill='both',
                                          expand=True)

    def _plot_title(self, job):
        """A one-line caption naming what is actually plotted."""
        task = job.get('task')
        cocktail = (job.get('config') or {}).get('ls_cocktail', '')
        if task == 'curve':
            model = 'stochastic' if job.get('Smodel', True) else 'analytical'
            return (f'{job.get("rad", "")}   kB = {job.get("kB", 0):g} '
                    f'cm/keV   {model}')
        if task == 'beta_spectrum':
            return (f'{job.get("rad", "")}   {job.get("mode", "beta-")}   '
                    f'level {job.get("level", "tot")}   (BetaShape)')
        if task == 'decay_data':
            return f'{job.get("rad", "")} — discrete emissions (PenNuc)'
        if task == 'stopping_power':
            return f'{job.get("particle", "electron")} stopping power   ' \
                   f'{cocktail}'
        if task == 'quenching':
            return f'Birks quenching, {job.get("particle", "electron")}s   ' \
                   f'{cocktail}'
        return ''

    def _plot(self, keys, rows, plot):
        """Draw *rows* as *plot* says.

        *plot* names the columns and the style: ``x``, ``y`` (a list; more
        than one overlays several series with a legend), optional ``labels``
        (legend text per ``y``), ``yerr`` (same length as ``y``, ``None``
        where there is no uncertainty), ``xlabel``, ``ylabel``,
        ``xlog``/``ylog``, ``annotate`` (a column to label at the first and
        last point, used for L on the efficiency curve), and ``kind`` --
        ``'stem'`` draws vertical lines per discrete emission, one colour per
        value of the ``group`` column, skipping rows where it is empty.
        """
        self._ensure_canvas()
        col = {k: i for i, k in enumerate(keys)}
        xk = plot['x']
        fig = self._figure
        fig.clear()
        ax = fig.add_subplot(111)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)

        if plot.get('kind') == 'stem':
            self._plot_stems(ax, col, rows, plot)
        else:
            self._plot_lines(ax, col, rows, plot)

        if plot.get('xlog'):
            ax.set_xscale('log')
        if plot.get('ylog'):
            ax.set_yscale('log')
        ax.set_xlabel(plot.get('xlabel', xk))
        ax.set_ylabel(plot.get('ylabel', ''))
        ax.grid(alpha=0.25, which='both')
        ax.set_title(self._plot_title(self._result_job or {}), fontsize=10,
                     color=NAVY, loc='left')
        fig.text(0.995, 0.005, f'TDCRPy {self.version}', ha='right',
                 va='bottom', fontsize=7, color='#9AA5B1')
        fig.tight_layout()
        self._canvas.draw()
        self.saveplot_button.configure(state='normal')

    def _plot_lines(self, ax, col, rows, plot):
        xk = plot['x']
        pts = [r for r in rows if r[col[xk]] is not None]
        x = [r[col[xk]] for r in pts]
        y_keys = plot['y']
        yerr_keys = plot.get('yerr') or [None] * len(y_keys)
        labels = plot.get('labels') or y_keys
        multi = len(y_keys) > 1
        for i, (yk, ek, lab) in enumerate(zip(y_keys, yerr_keys, labels)):
            color = PLOT_COLORS[i % len(PLOT_COLORS)]
            y = [r[col[yk]] for r in pts]
            if ek and any(r[col[ek]] is not None for r in pts):
                yerr = [r[col[ek]] or 0.0 for r in pts]
                ax.errorbar(x, y, yerr=yerr, fmt='o-', ms=3, lw=1.2,
                            capsize=2, color=color,
                            label=lab if multi else None)
            else:
                style = 'o-' if len(x) <= 60 else '-'
                ax.plot(x, y, style, ms=3, lw=1.4, color=color,
                        label=lab if multi else None)
        if multi:
            ax.legend(fontsize=8, framealpha=0.92, edgecolor=RULE)
        akey = plot.get('annotate')
        if akey and pts:
            y0 = [r[col[y_keys[0]]] for r in pts]
            # The last label goes left of its point, or it falls off the axes.
            for r, yv, dx, ha in ((pts[0], y0[0], 6, 'left'),
                                 (pts[-1], y0[-1], -6, 'right')):
                ax.annotate(f'{akey} = {r[col[akey]]:.3g}', (r[col[xk]], yv),
                            textcoords='offset points', xytext=(dx, -12),
                            ha=ha, fontsize=8, color=MUTED)

    def _plot_stems(self, ax, col, rows, plot):
        xk, yk, gk = plot['x'], plot['y'][0], plot['group']
        groups = []
        for r in rows:
            g = r[col[gk]]
            if g and g not in groups:
                groups.append(g)
        positive = [r[col[yk]] for r in rows
                    if r[col[gk]] and r[col[yk]] and r[col[yk]] > 0]
        if not positive:
            ax.text(0.5, 0.5, 'No discrete emissions in the decay data.',
                    ha='center', va='center', transform=ax.transAxes,
                    color=MUTED)
            return
        base = min(positive) / 3.0 if plot.get('ylog') else 0.0
        for i, g in enumerate(groups):
            pts = [r for r in rows if r[col[gk]] == g
                   and r[col[yk]] and r[col[yk]] > 0]
            x = [r[col[xk]] for r in pts]
            y = [r[col[yk]] for r in pts]
            color = PLOT_COLORS[i % len(PLOT_COLORS)]
            ax.vlines(x, base, y, colors=color, lw=1.2, label=g)
            ax.plot(x, y, 'o', ms=3, color=color)
        ax.set_ylim(bottom=base)
        ax.legend(fontsize=8, framealpha=0.92, edgecolor=RULE,
                  loc='upper center')

    def _save_plot(self):
        if self._figure is None or not self._figure.axes:
            messagebox.showinfo('TDCRPy', 'There is no plot to save yet.')
            return
        job = self._result_job or {}
        stem = '_'.join(str(p) for p in (job.get('task', 'plot'),
                                         job.get('rad') or
                                         job.get('particle', ''))
                        if p).replace(', ', '-')
        path = filedialog.asksaveasfilename(
            title='Save plot', defaultextension='.png',
            initialfile=f'{stem}.png',
            filetypes=[('PNG image', '*.png'), ('PDF document', '*.pdf'),
                       ('SVG drawing', '*.svg')])
        if not path:
            return
        try:
            self._figure.savefig(path, dpi=300, bbox_inches='tight')
        except (OSError, ValueError) as exc:
            messagebox.showerror('TDCRPy', f'Could not save the plot:\n{exc}')
            return
        self._append_log(f'plot written to {path}')

    # ------------------------------------------------------------- endings
    def _handle_exit(self, rc):
        cancelled, self._cancelled = self._cancelled, False
        self._worker = None
        self.progress.stop()
        self.progress.configure(mode='determinate', value=0)
        self._set_busy(False)
        if cancelled:
            self.status.configure(text='cancelled')
            return
        if rc == 0:
            if self.status.cget('text') in ('reading configuration…',
                                            'running…'):
                self.status.configure(text='ready')
            return
        # A negative code is a signal on POSIX; on Windows a native crash in
        # the numba kernel shows up as an NTSTATUS such as 0xC0000005.
        if rc < 0 or rc in (128 + 9, 139, 3221225477):
            self.status.configure(text=f'crashed (code {rc})')
            self._append_log(f'the worker crashed (code {rc}); this is a '
                             f'fault in the Monte-Carlo kernel, not in the '
                             f'inputs, and re-running usually succeeds')
        else:
            self.status.configure(text=f'failed (code {rc})')
        self.results.select(self.tab_log)

    def _save(self):
        if len(self._rows) < 2:
            messagebox.showinfo('TDCRPy', 'There are no results to export '
                                          'yet.')
            return
        job = self._result_job or {}
        stem = '_'.join(str(p) for p in (job.get('task', 'results'),
                                         job.get('rad') or
                                         job.get('particle', ''))
                        if p).replace(', ', '-')
        path = filedialog.asksaveasfilename(
            title='Export results as CSV', defaultextension='.csv',
            initialfile=f'{stem}.csv',
            filetypes=[('CSV', '*.csv'), ('All files', '*.*')])
        if not path:
            return
        meta = [f'TDCRPy {self.version}',
                f'{task_label(job.get("task", ""))}',
                self._run_summary(job) if job else '',
                f'exported {datetime.datetime.now().isoformat(timespec="seconds")}']
        cfg = job.get('config') or {}
        if cfg:
            meta.append('configuration: ' + '; '.join(
                f'{k}={v}' for k, v in sorted(cfg.items())))
        meta += provenance_lines(self._result_info)
        try:
            write_csv(path, self._rows, [m for m in meta if m])
        except OSError as exc:
            messagebox.showerror('TDCRPy', f'Could not save:\n{exc}')
            return
        self._append_log(f'results written to {path}')

    def _quit(self):
        if self._worker is not None and self._last_job is not None:
            if not messagebox.askyesno(
                    'TDCRPy', 'A calculation is running. Stop it and quit?'):
                return
        if self._worker is not None:
            self._worker.cancel()
        self.master.destroy()

    def _on_tk_error(self, exc_type, exc, tb):
        """An error in the window's own code: say so instead of vanishing
        (under pythonw there is no console for the traceback)."""
        text = ''.join(traceback.format_exception(exc_type, exc, tb))
        try:
            self._append_log(text)
        except tk.TclError:
            pass
        messagebox.showerror('TDCRPy — internal error',
                             f'{exc_type.__name__}: {exc}\n\nThe full '
                             f'traceback is in the Log tab. Please report it '
                             f'at {ISSUES_URL}')

    # ------------------------------------------------------------- dialogs
    def _dialog(self, title):
        win = tk.Toplevel(self.master)
        win.title(title)
        win.transient(self.master)
        win.resizable(False, False)
        win.configure(bg=WHITE)
        win.bind('<Escape>', lambda e: win.destroy())
        return win

    def _link(self, parent, text, url):
        lab = tk.Label(parent, text=text, fg=CYAN, bg=WHITE, cursor='hand2',
                       font=self.f_sub)
        lab.bind('<Button-1>', lambda e: self._open(url))
        return lab

    def _show_about(self):
        win = self._dialog('About TDCRPy')
        pad = self._px(22)
        logo = self.img.get('tdcrpy_logo_about')
        if logo is not None:
            tk.Label(win, image=logo, bg=WHITE).pack(pady=(pad, 4))
        tk.Label(win, text=f'TDCRPy {self.version}', font=self.f_title,
                 fg=NAVY, bg=WHITE).pack()
        tk.Label(win, text='A photo-physical stochastic model for liquid '
                           'scintillation counting\n(TDCR and CIEMAT/NIST '
                           'methods)',
                 font=self.f_sub, fg=INK, bg=WHITE, justify='center').pack(
            pady=(4, 12))
        try:
            author = metadata.metadata('TDCRPy').get('Author')
        except metadata.PackageNotFoundError:
            author = None
        if author:
            tk.Label(win, text=f'Author: {author}', font=self.f_sub,
                     fg=MUTED, bg=WHITE).pack()
        bipm = self.img.get('bipm_logo')
        if bipm is not None:
            mark = tk.Label(win, image=bipm, bg=WHITE, cursor='hand2')
            mark.pack(pady=(12, 2))
            mark.bind('<Button-1>', lambda e: self._open(BIPM_URL))
        tk.Label(win, text='Developed and maintained by the Bureau '
                           'International des Poids et Mesures',
                 font=self.f_sub, fg=MUTED, bg=WHITE).pack(pady=(2, 10))
        links = tk.Frame(win, bg=WHITE)
        links.pack()
        for text, url in (('Documentation', DOCS_URL),
                          ('Reference paper', PAPER_URL),
                          ('Source code', REPO_URL),
                          ('Report an issue', ISSUES_URL)):
            self._link(links, text, url).pack(side='left', padx=8)
        tk.Label(win, text='Released under the MIT licence. Provided "as is", '
                           'without warranty of any kind.',
                 font=self.f_sub, fg=MUTED, bg=WHITE).pack(pady=(12, 4))
        buttons = tk.Frame(win, bg=WHITE)
        buttons.pack(pady=(6, pad))
        ttk.Button(buttons, text='Licence…',
                   command=lambda: (win.destroy(), self._show_licence())
                   ).pack(side='left', padx=4)
        ttk.Button(buttons, text='Check for updates',
                   command=lambda: (win.destroy(), self._check_updates())
                   ).pack(side='left', padx=4)
        ttk.Button(buttons, text='Close', command=win.destroy).pack(
            side='left', padx=4)
        win.grab_set()

    # -------------------------------------------------------------- updates
    def _check_updates(self):
        """Ask PyPI, off the main thread, for the latest release."""
        if self._update_check is not None and self._update_check.is_alive():
            return
        self.status.configure(text='checking for updates…')
        q = self._queue

        def ask():
            try:
                latest, date = fetch_latest_release()
                q.put(('update', {'latest': latest, 'date': date}))
            except Exception as exc:       # network, proxy, PyPI down...
                q.put(('update', {'error': f'{type(exc).__name__}: {exc}'}))

        self._update_check = threading.Thread(target=ask, daemon=True)
        self._update_check.start()

    def _show_update(self, info):
        if self.status.cget('text') == 'checking for updates…':
            self.status.configure(text='ready')
        if info.get('error'):
            messagebox.showerror(
                'Check for updates', 'Could not reach PyPI to check for a '
                f'new version.\n\n{info["error"]}\n\nCheck the network '
                f'connection (or proxy), or look at {PYPI_URL}')
            return
        latest, current = info['latest'], self.version
        self._append_log(f'latest TDCRPy on PyPI: {latest} '
                         f'({info.get("date") or "date unknown"}); '
                         f'this one: {current}')
        if not is_newer(latest, current):
            messagebox.showinfo('Check for updates',
                                f'TDCRPy {current} is up to date — the latest '
                                f'release on PyPI is {latest}.')
            return

        kind = install_kind()
        win = self._dialog('Update available')
        pad = self._px(20)
        tk.Label(win, text=f'TDCRPy {latest} is available', font=self.f_title,
                 fg=NAVY, bg=WHITE).pack(padx=pad, pady=(pad, 4))
        when = f'Released {info["date"]}. ' if info.get('date') else ''
        tk.Label(win, text=f'{when}You have {current}.', font=self.f_sub,
                 fg=MUTED, bg=WHITE).pack(padx=pad)
        if kind == 'pip':
            note = ('Update now closes this window, installs the new version '
                    'with pip and reopens it — about a minute.\nClose any '
                    'notebook or script using TDCRPy first: Windows locks '
                    'the files it has open.')
        elif kind == 'editable':
            note = ('This TDCRPy runs from a source checkout (an editable '
                    'install), so pip must not replace it.\nUpdate it with '
                    '"git pull" in that folder.')
        else:
            note = (f'This TDCRPy was installed by {kind}; update it with '
                    f'{kind}, so that {kind} keeps track of it.')
        tk.Label(win, text=note, font=self.f_sub, fg=INK, bg=WHITE,
                 justify='left', wraplength=self._px(420)).pack(
            padx=pad, pady=(12, 8))
        self._link(win, 'What changed (commit history)', HISTORY_URL).pack()
        buttons = tk.Frame(win, bg=WHITE)
        buttons.pack(pady=(14, pad))
        if kind == 'pip':
            ttk.Button(buttons, text='Update now',
                       command=lambda: (win.destroy(),
                                        self._start_update(latest))
                       ).pack(side='left', padx=4)
        ttk.Button(buttons, text='Later', command=win.destroy).pack(
            side='left', padx=4)
        win.grab_set()

    def _start_update(self, version):
        if self._worker is not None:
            messagebox.showinfo('TDCRPy', 'Wait for the calculation to '
                                          'finish, or cancel it, then update.')
            return
        try:
            log = launch_updater(version, self.version)
        except Exception as exc:
            messagebox.showerror('TDCRPy', f'Could not start the update:'
                                           f'\n{exc}')
            return
        self._append_log(f'update to {version} handed over; log: {log}')
        self.master.destroy()

    def _show_licence(self):
        text = licence_text()
        win = self._dialog('TDCRPy — licence')
        if text is None:
            tk.Label(win, text='TDCRPy is released under the MIT licence.',
                     bg=WHITE, fg=INK).pack(padx=20, pady=(20, 6))
            self._link(win, 'Read it online', LICENCE_URL).pack(pady=(0, 20))
            return
        box = tk.Text(win, width=78, height=24, wrap='word', relief='flat',
                      font=('Consolas', 9), padx=12, pady=10, bg=STRIPE)
        box.insert('1.0', text)
        box.configure(state='disabled')
        box.pack(padx=12, pady=(12, 6))
        ttk.Button(win, text='Close', command=win.destroy).pack(pady=(0, 12))


class _Tooltip:
    """A small hint shown while the pointer rests on a widget."""

    def __init__(self, widget, text):
        self.widget, self.text, self.tip = widget, text, None
        widget.bind('<Enter>', self._show, add='+')
        widget.bind('<Leave>', self._hide, add='+')

    def _show(self, _event=None):
        if self.tip is not None:
            return
        x = self.widget.winfo_rootx() + 10
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 4
        self.tip = tk.Toplevel(self.widget)
        self.tip.wm_overrideredirect(True)
        self.tip.wm_geometry(f'+{x}+{y}')
        tk.Label(self.tip, text=self.text, bg='#FFFFE8', fg=INK,
                 relief='solid', bd=1, padx=6, pady=3,
                 wraplength=320, justify='left').pack()

    def _hide(self, _event=None):
        if self.tip is not None:
            self.tip.destroy()
            self.tip = None


def _report(text, title='TDCRPy'):
    """Say something without a console to say it on.

    ``tdcrpy-gui`` is installed as a *gui_script* so that Windows does not
    put a terminal behind the window. Such a launcher runs under
    ``pythonw.exe``, where ``sys.stdout`` and ``sys.stderr`` are ``None`` --
    so a plain ``print`` here would raise and the launcher would fail with
    nothing shown at all. Fall back to a dialog, which is the only channel
    that always exists once Tk is importable.
    """
    stream = sys.stderr or sys.stdout
    if stream is not None:
        try:
            stream.write(text + '\n')
            return
        except (OSError, ValueError):
            pass
    try:
        root = tk.Tk()
        root.withdraw()
        messagebox.showinfo(title, text)
        root.destroy()
    except tk.TclError:
        pass


def _dpi_aware():
    """Ask Windows for real pixels, so text is sharp on a scaled display
    instead of bitmap-stretched. Harmless elsewhere."""
    if os.name != 'nt':
        return
    try:
        import ctypes
        ctypes.windll.shcore.SetProcessDpiAwareness(1)
    except (AttributeError, OSError):
        pass


def main(argv=None):
    """Open the window. Returns the process exit status."""
    argv = list(sys.argv[1:] if argv is None else argv)
    # A launcher that opens a window regardless of its arguments is a
    # nuisance: --help should answer, not launch.
    if any(a in ('-h', '--help') for a in argv):
        _report(__doc__.strip(), 'TDCRPy — help')
        return 0
    if any(a in ('-V', '--version') for a in argv):
        _report(f'TDCRPy {installed_version()}', 'TDCRPy — version')
        return 0
    if argv:
        _report(f'tdcrpy-gui takes no arguments (got {" ".join(argv)}). '
                f'Try --help.')
        return 2

    if package_dir() is None:
        _report('TDCRPy is not importable from this interpreter.')
        return 1
    _dpi_aware()
    try:
        root = tk.Tk(className='TDCRPy')
    except tk.TclError as exc:
        # Nearly always a Python built without tcl/tk, or a headless session.
        _report(f'Could not open a window: {exc}\n\n'
                f'On Windows this usually means Python was installed without '
                f'the "tcl/tk and IDLE" option. Check with:\n'
                f'    python -c "import tkinter"')
        return 1
    app = App(root)
    sw, sh = root.winfo_screenwidth(), root.winfo_screenheight()
    w, h = min(app._px(1360), sw - 40), min(app._px(860), sh - 90)
    root.minsize(min(app._px(1100), w), min(app._px(680), h))
    root.geometry(f'{w}x{h}+{max(0, (sw - w) // 2)}+{max(0, (sh - h - 40) // 2)}')
    root.mainloop()
    return 0


if __name__ == '__main__':
    sys.exit(main())
