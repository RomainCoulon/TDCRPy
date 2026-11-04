# -*- coding: utf-8 -*-
"""A window for driving TDCRPy, for people who would rather not write a script.

Start it with ``tdcrpy-gui``, or ``python -m tdcrpy.gui``.

Six calculations are offered: efficiencies at a light yield you supply
(:func:`tdcrpy.TDCRPy.TDCRPy`); the light yield fitted to a measured TDCR
ratio by either the stochastic (:func:`~tdcrpy.TDCRPy.eff`) or the analytical
(:func:`~tdcrpy.TDCRPy.effA`) model, both reported over the same seven
channels; the efficiency curve, ε_D against TDCR, over a range of light
yields; the emitted beta spectrum of one nuclide and transition
(:func:`~tdcrpy.TDCR_model_lib.readBetaShape`); and the electron stopping
power dE/dx against energy (:func:`~tdcrpy.TDCR_model_lib.stoppingpower`),
for the configured scintillator or for every model side by side.

Design notes
------------
This module deliberately does not import ``tdcrpy``. Every calculation, and
the initial enquiry that fills the widgets, runs in a separate process
(:mod:`tdcrpy._gui_worker`), so the window opens at once instead of waiting
for numba, stays responsive during a run that may take minutes, can cancel a
run by killing it, and survives a native crash in the kernel.

Nothing here keeps its own copy of the cocktail list, the nuclide list or the
configuration defaults: all three are read from the library at start-up. A
hardcoded copy would be one more thing to go stale, which this package has
been bitten by repeatedly.
"""
import json
import os
import queue
import re
import shutil
import subprocess
import sys
import threading
import tkinter as tk
from importlib.util import find_spec
from tkinter import filedialog, messagebox, ttk

#: Matches the percentage tqdm prints, so the progress bar can follow the
#: Monte-Carlo loop without the model needing to report anything.
_TQDM_PCT = re.compile(r'(\d{1,3})%')

TASKS = (
    ('efficiency', 'Efficiency at L'),
    ('fit', 'Fit L to a TDCR — stochastic'),
    ('fit_analytical', 'Fit L to a TDCR — analytical'),
    ('curve', 'Curve ε_D vs TDCR over an L range'),
    ('beta_spectrum', 'Beta spectrum'),
    ('stopping_power', 'Stopping power'),
)

#: Most points any plotted grid (L, or energy for the stopping power) may
#: have; on the efficiency curve each one costs one replay.
CURVE_MAX_POINTS = 500


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
        raise InputError('Nuclide: give at least one, e.g. "Fe-55".')
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
    """Energies (keV) to evaluate the stopping power at."""
    return _build_grid(text, 'Energy range', log_spacing)


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

    if task == 'stopping_power':
        # No nuclide, no Monte Carlo: dE/dx is a property of the
        # scintillator (set via 'config') and the model alone.
        return {
            'task': task,
            'E_grid': build_E_grid(form.get('Erange', ''),
                                   bool(form.get('Elog', True))),
            'compare_models': bool(form.get('spCompare', False)),
            'config': form.get('config') or {},
        }

    if task == 'beta_spectrum':
        rad = parse_nuclides(form['rad'])[0]
        return {
            'task': task,
            'rad': rad,
            'mode': form.get('betaMode') or 'beta-',
            'level': (form.get('betaLevel') or 'tot').strip(),
            'config': form.get('config') or {},
        }

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
        'uncData': bool(form.get('uncData')),
        'config': form.get('config') or {},
    }

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
    val = f'{value:.6g}'
    if uncertainty in (None, ''):
        return val, '', ''
    unc = f'{uncertainty:.3g}'
    rel = f'{100.0 * uncertainty / value:.2f} %' if value else ''
    return val, unc, rel


def package_dir():
    """Locate the installed package without importing it.

    ``find_spec`` reads the metadata but does not execute the module, so the
    window can find ``config.toml`` without paying for the numba import.
    """
    spec = find_spec('tdcrpy')
    if spec is None or not spec.submodule_search_locations:
        return None
    return spec.submodule_search_locations[0]


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
        super().__init__(master, padding=8)
        self.master.title('TDCRPy')
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
        self._plot_win = None

        self.vars = {}
        self._build_header()
        self._build_inputs()
        self._build_output()
        self._build_footer()

        self._snapshot_config()
        self._on_task_change()
        self.after(80, self._drain)
        self._enquire()

    # ---------------------------------------------------------------- setup
    def _var(self, name, value=''):
        v = tk.StringVar(value=str(value))
        self.vars[name] = v
        return v

    def _flag(self, name, value=False):
        v = tk.BooleanVar(value=bool(value))
        self.vars[name] = v
        return v

    def _build_header(self):
        bar = ttk.Frame(self)
        bar.grid(row=0, column=0, columnspan=2, sticky='ew', pady=(0, 8))
        ttk.Label(bar, text='Calculation', font=('', 10, 'bold')).pack(
            side='left', padx=(0, 10))
        self._var('task', TASKS[0][0])
        for key, label in TASKS:
            ttk.Radiobutton(bar, text=label, value=key,
                            variable=self.vars['task'],
                            command=self._on_task_change).pack(side='left',
                                                               padx=6)
        self.version_label = ttk.Label(bar, text='starting…',
                                       foreground='#666')
        self.version_label.pack(side='right')

    def _build_inputs(self):
        nb = ttk.Notebook(self)
        nb.grid(row=1, column=0, sticky='nsew', padx=(0, 8))
        self._tab_source(nb)
        self._tab_scintillator(nb)
        self._tab_micelles(nb)
        self._tab_delayed(nb)
        self._tab_counter(nb)

    def _row(self, parent, r, label, name, default='', width=16, hint=None):
        ttk.Label(parent, text=label).grid(row=r, column=0, sticky='w',
                                           pady=2)
        e = ttk.Entry(parent, textvariable=self._var(name, default),
                      width=width)
        e.grid(row=r, column=1, sticky='w', pady=2)
        if hint:
            ttk.Label(parent, text=hint, foreground='#666').grid(
                row=r, column=2, sticky='w', padx=(8, 0))
        return e

    def _tab_source(self, nb):
        f = ttk.Frame(nb, padding=10)
        nb.add(f, text='Source')

        ttk.Label(f, text='Nuclide(s)').grid(row=0, column=0, sticky='w')
        self.nuclide_entry = ttk.Entry(
            f, textvariable=self._var('rad', 'Fe-55'), width=28)
        self.nuclide_entry.grid(row=0, column=1, sticky='w', pady=2)

        picker = ttk.Frame(f)
        picker.grid(row=1, column=1, sticky='w', pady=(0, 6))
        self.nuclide_box = ttk.Combobox(picker, width=14, state='disabled',
                                        values=())
        self.nuclide_box.pack(side='left')
        ttk.Button(picker, text='Add', width=5,
                   command=self._add_nuclide).pack(side='left', padx=4)
        ttk.Label(f, text='comma-separated for a mixture',
                  foreground='#666').grid(row=0, column=2, sticky='w',
                                          padx=(8, 0))

        self.pmf_entry = self._row(
            f, 2, 'Activity fractions', 'pmf', '',
            hint='blank for a single nuclide, else e.g. 0.8, 0.2')
        ttk.Separator(f, orient='horizontal').grid(
            row=3, column=0, columnspan=3, sticky='ew', pady=8)

        self.l_entry = self._row(f, 4, 'Light yield L', 'L', '1.0',
                                 hint='keV⁻¹; 3 values for '
                                      'asymmetric PMTs')
        self.td_entry = self._row(f, 5, 'Measured TDCR', 'TD', '0.9',
                                  hint='T/D; 4 values for an asymmetric '
                                       'counter')
        self.guess_entry = self._row(f, 6, 'Initial L guess', 'L_guess', '1.0')
        self.bounds_entry = self._row(
            f, 7, 'L bounds', 'Lbounds', '',
            hint='blank lets the fit bracket the physical root itself')

        self.range_entry = self._row(
            f, 8, 'L range', 'Lrange', '0.1, 10, 25',
            hint='from, to, number of points (keV⁻¹)')
        self.log_check = ttk.Checkbutton(
            f, text='log spacing', variable=self._flag('Llog', True))
        self.log_check.grid(row=9, column=1, sticky='w')

        self.kb_entry = self._row(f, 10, 'Birks constant kB', 'kB', '1.0e-5',
                                  hint='cm keV⁻¹')
        self.v_entry = self._row(f, 11, 'Volume', 'V', '10', hint='mL')
        self.n_entry = self._row(f, 12, 'Number of decays N', 'N', '10000',
                                 hint='≥ 10 000 recommended (JCGM 101)')

        self.smodel_check = ttk.Checkbutton(
            f, text='Stochastic model (uncheck for the analytical one, '
                    'valid for pure β emitters only)',
            variable=self._flag('Smodel', True))
        self.smodel_check.grid(row=13, column=0, columnspan=3, sticky='w',
                               pady=(8, 0))
        ttk.Checkbutton(f, text='Propagate nuclear-data uncertainty',
                        variable=self._flag('uncData')).grid(
            row=14, column=0, columnspan=3, sticky='w')
        ttk.Checkbutton(f, text='Verbose log (prints every decay)',
                        variable=self._flag('display')).grid(
            row=15, column=0, columnspan=3, sticky='w')

        ttk.Separator(f, orient='horizontal').grid(
            row=16, column=0, columnspan=3, sticky='ew', pady=8)
        ttk.Label(f, text='Beta decay mode').grid(row=17, column=0,
                                                   sticky='w', pady=2)
        self.beta_mode_box = ttk.Combobox(
            f, textvariable=self._var('betaMode', 'beta-'), width=10,
            values=('beta-', 'beta+'), state='readonly')
        self.beta_mode_box.grid(row=17, column=1, sticky='w', pady=2)
        self.beta_level_entry = self._row(
            f, 18, 'Transition level', 'betaLevel', 'tot',
            hint="the daughter level before the transition, or 'tot' "
                 "for a pure beta emitter")

    #: Where the stopping-power choice acts. Said in the window because the
    #: answer -- not everywhere -- is not what the field suggests.
    SP_NOTE = ('Electrons, low energy. It enters wherever the quenching '
               'integral is evaluated directly: below {where}, and for any '
               'partial deposit. Above it the shipped quenching tables are '
               'interpolated and are not recomputed.')

    def _tab_scintillator(self, nb):
        f = ttk.Frame(nb, padding=10)
        nb.add(f, text='Scintillator')
        ttk.Label(f, text='Cocktail').grid(row=0, column=0, sticky='w')
        self.cocktail_box = ttk.Combobox(
            f, textvariable=self._var('ls_cocktail'), width=26,
            state='disabled', values=())
        self.cocktail_box.grid(row=0, column=1, sticky='w', pady=2)
        self._row(f, 1, 'Aqueous fraction fAq', 'fAq', '0.0',
                  hint='0 for an organic source')
        ttk.Label(f, text='Aqueous solvent').grid(row=2, column=0, sticky='w')
        ttk.Combobox(f, textvariable=self._var('solvantType', 'Water'),
                     values=('Water', 'HCl', 'HNO3', 'NaOH'),
                     width=14).grid(row=2, column=1, sticky='w', pady=2)
        self._row(f, 3, 'Solvent concentration', 'solvantConc', '0.0',
                  hint='mol L⁻¹')
        ttk.Separator(f, orient='horizontal').grid(
            row=4, column=0, columnspan=3, sticky='ew', pady=8)
        ttk.Label(f, text='Stopping-power model').grid(row=5, column=0,
                                                       sticky='w')
        # Editable until the library says which models it has; read-only
        # afterwards, so a typo cannot reach config.toml.
        self.sp_box = ttk.Combobox(
            f, textvariable=self._var('sp_model', 'tan_xia'), width=20,
            values=())
        self.sp_box.grid(row=5, column=1, sticky='w', pady=2)
        self.sp_note = ttk.Label(f, wraplength=430, foreground='#666',
                                 text=self.SP_NOTE.format(where='Einterp_e'))
        self.sp_note.grid(row=6, column=0, columnspan=3, sticky='w',
                          pady=(0, 6))
        self._row(f, 7, 'Chou parameter', 'chou_param', '0',
                  hint='cm² MeV⁻²')
        self._row(f, 8, 'Electron quench bins', 'nE_electron', '1000')
        self._row(f, 9, 'Alpha quench bins', 'nE_alpha', '1000')

        ttk.Separator(f, orient='horizontal').grid(
            row=10, column=0, columnspan=3, sticky='ew', pady=8)
        self.erange_entry = self._row(
            f, 11, 'Energy range (plot)', 'Erange', '0.1, 2000, 100',
            hint='from, to, number of points (keV)')
        self.elog_check = ttk.Checkbutton(
            f, text='log spacing', variable=self._flag('Elog', True))
        self.elog_check.grid(row=12, column=1, sticky='w')
        self.spcompare_check = ttk.Checkbutton(
            f, text='overlay every stopping-power model (same ρ, Z, A)',
            variable=self._flag('spCompare', False))
        self.spcompare_check.grid(row=13, column=0, columnspan=3, sticky='w')

    def _tab_micelles(self, nb):
        f = ttk.Frame(nb, padding=10)
        nb.add(f, text='Micelles')
        ttk.Label(
            f, wraplength=430, foreground='#444',
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
        f = ttk.Frame(nb, padding=10)
        nb.add(f, text='Delayed fluorescence')
        ttk.Label(
            f, wraplength=430, foreground='#444',
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
        f = ttk.Frame(nb, padding=10)
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

    def _build_output(self):
        right = ttk.Frame(self)
        right.grid(row=1, column=1, sticky='nsew')
        right.columnconfigure(0, weight=1)
        right.rowconfigure(0, weight=3)
        right.rowconfigure(2, weight=2)

        self.tree = ttk.Treeview(right, show='headings', height=10)
        self._set_columns(self.QUANTITY_COLUMNS)
        self.tree.grid(row=0, column=0, sticky='nsew')
        ttk.Scrollbar(right, orient='vertical',
                      command=self.tree.yview).grid(row=0, column=1,
                                                    sticky='ns')

        ttk.Label(right, text='Log', foreground='#666').grid(
            row=1, column=0, sticky='w', pady=(8, 2))
        self.log = tk.Text(right, height=10, wrap='none',
                           font=('Consolas', 9), state='disabled')
        self.log.grid(row=2, column=0, sticky='nsew')
        ttk.Scrollbar(right, orient='vertical',
                      command=self.log.yview).grid(row=2, column=1,
                                                   sticky='ns')

    #: Results table for a single calculation: one quantity per row.
    QUANTITY_COLUMNS = (
        ('quantity', 'Quantity', 210, 'w'),
        ('value', 'Value', 120, 'e'),
        ('u', 'u', 100, 'e'),
        ('rel', 'u / value', 90, 'e'),
    )

    def _set_columns(self, spec):
        """Rebuild the results table's columns from (key, title, width,
        anchor) tuples."""
        self.tree.configure(columns=[c[0] for c in spec])
        for key, title, width, anchor in spec:
            self.tree.heading(key, text=title)
            self.tree.column(key, width=width, anchor=anchor, stretch=True)

    def _build_footer(self):
        bar = ttk.Frame(self)
        bar.grid(row=2, column=0, columnspan=2, sticky='ew', pady=(8, 0))
        bar.columnconfigure(3, weight=1)

        self.run_button = ttk.Button(bar, text='Run', command=self._run)
        self.run_button.grid(row=0, column=0)
        self.cancel_button = ttk.Button(bar, text='Cancel', state='disabled',
                                        command=self._cancel)
        self.cancel_button.grid(row=0, column=1, padx=4)
        ttk.Button(bar, text='Save results…',
                   command=self._save).grid(row=0, column=2, padx=4)

        self.progress = ttk.Progressbar(bar, mode='determinate', length=200)
        self.progress.grid(row=0, column=3, sticky='e', padx=8)
        self.status = ttk.Label(bar, text='ready')
        self.status.grid(row=0, column=4, sticky='e')

        tools = ttk.Frame(self)
        tools.grid(row=3, column=0, columnspan=2, sticky='ew', pady=(6, 0))
        ttk.Button(tools, text='Undo my config changes',
                   command=self._restore_snapshot).pack(side='left')
        ttk.Button(tools, text='Restore shipped defaults',
                   command=self._restore_defaults).pack(side='left', padx=6)
        ttk.Label(
            tools, foreground='#666', wraplength=700,
            text='Running writes these settings to the package config.toml, '
                 'which other scripts share. The snapshot taken when this '
                 'window opened can be put back at any time.',
        ).pack(side='left', padx=8)

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
        """Show only the fields the chosen calculation uses."""
        task = self.vars['task'].get()
        fitting = task in ('fit', 'fit_analytical')
        curve = task == 'curve'
        beta = task == 'beta_spectrum'
        stopping = task == 'stopping_power'
        # Neither reads a nuclide mixture, runs a Monte Carlo, nor needs a
        # counting geometry: both are a direct plot of a library function.
        diagnostic = beta or stopping
        for widget, shown in (
            (self.l_entry, task == 'efficiency'),
            (self.td_entry, fitting),
            (self.guess_entry, fitting),
            (self.bounds_entry, fitting),
            (self.range_entry, curve),
            (self.log_check, curve),
            (self.pmf_entry, not diagnostic),
            (self.kb_entry, not diagnostic),
            (self.v_entry, not diagnostic),
            (self.n_entry, task != 'fit_analytical' and not diagnostic),
            (self.smodel_check, not fitting and not diagnostic),
            (self.nuclide_entry, not stopping),
            (self.beta_level_entry, beta),
            (self.erange_entry, stopping),
            (self.elog_check, stopping),
            (self.spcompare_check, stopping),
        ):
            widget.configure(state='normal' if shown else 'disabled')
        # A combobox, not an entry: 'normal' would let free text through.
        self.beta_mode_box.configure(state='readonly' if beta else 'disabled')

    def _add_nuclide(self):
        pick = self.nuclide_box.get().strip()
        if not pick:
            return
        current = self.vars['rad'].get().strip()
        parts = [p.strip() for p in current.split(',') if p.strip()]
        if pick not in parts:
            parts.append(pick)
        self.vars['rad'].set(', '.join(parts))

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
        the diagnostic tasks (beta spectrum, stopping power) carry none of
        'rad', 'kB', 'V' or 'N'."""
        parts = [f'run: {job["task"]}']
        if 'rad' in job:
            parts.append(job['rad'])
        if 'kB' in job:
            parts.append(f'kB={job["kB"]:g}')
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
        return '  '.join(parts)

    def _run(self):
        if self._worker is not None:
            messagebox.showinfo(
                'TDCRPy', 'A calculation is already running. Wait for it, or '
                          'press Cancel.')
            return
        form = {k: v.get() for k, v in self.vars.items()}
        form['config'] = self.collect_config()
        try:
            job = build_job(form)
        except InputError as exc:
            messagebox.showerror('Check the inputs', str(exc))
            return

        self.tree.delete(*self.tree.get_children())
        self._rows = []
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
            self.status.configure(text=msg.get('text', ''))
            self._append_log(f'• {msg.get("text", "")}')
        elif event == 'progress':
            if str(self.progress.cget('mode')) == 'indeterminate':
                self.progress.stop()
                self.progress.configure(mode='determinate')
            self.progress.configure(value=min(100.0, float(msg['pct'])))
        elif event == 'describe':
            self._apply_describe(msg)
        elif event == 'result':
            if msg.get('columns'):
                self._show_curve(msg['columns'], msg.get('rows') or [],
                                 msg.get('plot') or {})
            else:
                self._show_result(msg.get('rows') or [])
        elif event == 'error':
            self._append_log(msg.get('traceback') or msg.get('message', ''))
            messagebox.showerror('TDCRPy failed', msg.get('message',
                                                          'unknown error'))

    def _apply_describe(self, msg):
        self._describe = msg
        self.version_label.configure(text=f'TDCRPy {msg.get("version", "?")}')
        cocktails = msg.get('cocktails') or []
        if cocktails:
            self.cocktail_box.configure(values=cocktails, state='readonly')
        nuclides = msg.get('nuclides') or []
        if nuclides:
            self.nuclide_box.configure(values=nuclides, state='readonly')
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
        self._append_log(f'configuration loaded from the installed package '
                         f'({len(nuclides)} nuclides, '
                         f'{len(cocktails)} cocktails)')

    def _show_result(self, rows):
        self._set_columns(self.QUANTITY_COLUMNS)
        self._rows = [['quantity', 'value', 'u', 'u_rel']]
        for quantity, value, unc in rows:
            cells = (quantity,) + format_value(value, unc)
            self.tree.insert('', 'end', values=cells)
            self._rows.append(list(cells))
        self.status.configure(text='done')

    def _show_curve(self, columns, rows, plot):
        """Fill the table with the curve, and draw it."""
        self._set_columns([(key, title, 82, 'e') for key, title in columns])
        # Full precision in the saved file; the table is for reading.
        self._rows = [[key for key, _title in columns]]
        for row in rows:
            self.tree.insert('', 'end', values=[
                '' if v is None else f'{v:.5g}' for v in row])
            self._rows.append(['' if v is None else repr(v) for v in row])
        self.status.configure(text='done')
        try:
            self._plot_curve([k for k, _t in columns], rows, plot)
        except ImportError as exc:
            self._append_log(f'no plot: matplotlib is not importable ({exc})')

    def _plot_title(self, job):
        """A one-line caption naming what is actually plotted."""
        task = job.get('task')
        if task == 'curve':
            model = 'stochastic' if job.get('Smodel', True) else 'analytical'
            return (f'{job.get("rad", "")}   kB = {job.get("kB", 0):g} '
                   f'cm/keV   {model}')
        if task == 'beta_spectrum':
            return (f'{job.get("rad", "")}   {job.get("mode", "beta-")}   '
                   f'level {job.get("level", "tot")}')
        if task == 'stopping_power':
            return 'electron stopping power'
        return ''

    def _plot_curve(self, keys, rows, plot):
        """One or more y-columns against one x-column, in a window of its
        own, with matplotlib's toolbar for zoom and Save.

        *plot* names which columns to draw and how: ``x``, ``y`` (a list --
        more than one overlays several series with a legend, used for
        'overlay every stopping-power model'), optionally ``yerr`` (same
        length as ``y``, entries ``None`` where there is no uncertainty),
        ``xlabel``, ``ylabel``, ``xlog``/``ylog``, and ``annotate`` (a column
        to label at the first and last point, used for L on the efficiency
        curve).

        Imported here, not at the top: the window is meant to open without
        paying for matplotlib, and most runs never draw a curve.
        """
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                       NavigationToolbar2Tk)
        col = {k: i for i, k in enumerate(keys)}
        xk = plot['x']
        pts = [r for r in rows if r[col[xk]] is not None]
        x = [r[col[xk]] for r in pts]
        y_keys = plot['y']
        yerr_keys = plot.get('yerr') or [None] * len(y_keys)
        multi = len(y_keys) > 1

        if self._plot_win is None or not self._plot_win.winfo_exists():
            win = tk.Toplevel(self.master)
            win.title('TDCRPy — plot')
            fig = Figure(figsize=(6.4, 4.6), dpi=100)
            canvas = FigureCanvasTkAgg(fig, master=win)
            NavigationToolbar2Tk(canvas, win).update()
            canvas.get_tk_widget().pack(fill='both', expand=True)
            self._plot_win, self._plot_fig, self._plot_canvas = win, fig, canvas
        fig = self._plot_fig
        fig.clear()
        ax = fig.add_subplot(111)
        for yk, ek in zip(y_keys, yerr_keys):
            y = [r[col[yk]] for r in pts]
            if ek and any(r[col[ek]] is not None for r in pts):
                yerr = [r[col[ek]] or 0.0 for r in pts]
                ax.errorbar(x, y, yerr=yerr, fmt='o-', ms=3, lw=1, capsize=2,
                            label=yk if multi else None)
            else:
                style = 'o-' if len(x) <= 60 else '-'
                ax.plot(x, y, style, ms=3, lw=1, label=yk if multi else None)
        if plot.get('xlog'):
            ax.set_xscale('log')
        if plot.get('ylog'):
            ax.set_yscale('log')
        ax.set_xlabel(plot.get('xlabel', xk))
        ax.set_ylabel(plot.get('ylabel', ''))
        ax.grid(alpha=0.3, which='both')
        if multi:
            ax.legend(fontsize=8)
        ax.set_title(self._plot_title(self._last_job or {}), fontsize=10)
        akey = plot.get('annotate')
        if akey and pts:
            y0 = [r[col[y_keys[0]]] for r in pts]
            # The last label goes left of its point, or it falls off the axes.
            for r, yv, dx, ha in ((pts[0], y0[0], 6, 'left'),
                                 (pts[-1], y0[-1], -6, 'right')):
                ax.annotate(f'{akey} = {r[col[akey]]:.3g}', (r[col[xk]], yv),
                            textcoords='offset points', xytext=(dx, -12),
                            ha=ha, fontsize=8, color='#555')
        fig.tight_layout()
        self._plot_canvas.draw()
        self._plot_win.lift()

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
            self.status.configure(text='done' if self._rows else 'ready')
            return
        # A crash inside the numba kernel is a real possibility in this
        # package, not a mishap to hide behind "failed": name it, so nobody
        # hunts for a mistake in their inputs.
        if rc < 0 or rc in (128 + 9, 139, 3221225477):
            self.status.configure(text=f'crashed (code {rc})')
            self._append_log(
                f'the worker died with code {rc}, which is a crash inside '
                f'the Monte-Carlo kernel rather than a problem with the '
                f'inputs. Re-running usually succeeds.')
        else:
            self.status.configure(text=f'failed (code {rc})')
            self._append_log(f'worker exited with code {rc}')

    def _save(self):
        if not self._rows:
            messagebox.showinfo('TDCRPy', 'There are no results to save yet.')
            return
        path = filedialog.asksaveasfilename(
            title='Save results', defaultextension='.csv',
            filetypes=[('CSV', '*.csv'), ('All files', '*.*')])
        if not path:
            return
        try:
            with open(path, 'w', encoding='utf-8', newline='') as fh:
                import csv
                # First row is the header, built with the table.
                csv.writer(fh).writerows(self._rows)
        except OSError as exc:
            messagebox.showerror('TDCRPy', f'Could not save:\n{exc}')
            return
        self._append_log(f'results written to {path}')


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


def main(argv=None):
    """Open the window. Returns the process exit status."""
    argv = list(sys.argv[1:] if argv is None else argv)
    # A launcher that opens a window regardless of its arguments is a
    # nuisance: --help should answer, not launch.
    if any(a in ('-h', '--help') for a in argv):
        _report(__doc__.strip(), 'TDCRPy — help')
        return 0
    if any(a in ('-V', '--version') for a in argv):
        from importlib.metadata import version
        _report(f'TDCRPy {version("TDCRPy")}', 'TDCRPy — version')
        return 0
    if argv:
        _report(f'tdcrpy-gui takes no arguments (got {" ".join(argv)}). '
                f'Try --help.')
        return 2

    if package_dir() is None:
        _report('TDCRPy is not importable from this interpreter.')
        return 1
    try:
        root = tk.Tk()
    except tk.TclError as exc:
        # Nearly always a Python built without tcl/tk, or a headless session.
        _report(f'Could not open a window: {exc}\n\n'
                f'On Windows this usually means Python was installed without '
                f'the "tcl/tk and IDLE" option. Check with:\n'
                f'    python -c "import tkinter"')
        return 1
    App(root)
    root.minsize(1040, 700)
    root.mainloop()
    return 0


if __name__ == '__main__':
    sys.exit(main())
