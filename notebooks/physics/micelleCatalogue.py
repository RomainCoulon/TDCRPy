# -*- coding: utf-8 -*-
"""Micelle confinement screening across TDCRPy's whole nuclide catalogue.

micelleCheck.ipynb measures six nuclides at their measured TDCR. This asks the
same physical question of all 230 nuclides in the PenNuc library, for which no
measured ratios exist: how much of the light is lost INSIDE the water droplets,
over and above the (1 - fAq) bulk dilution that a refitted L absorbs?

    F          = sum(retained energy) / sum(quenched energy)   -- energy-weighted
    dilution   = fAq                        (universal, degenerate with L)
    confinement = (1 - fAq - F) / (1 - fAq) -- the part L CANNOT absorb

Confinement is the quantity the literature calls "the micelle effect", and it
is the one that biases an efficiency. It needs no TDCR fit: one recorded run
per nuclide gives the Birks-quenched electron spectrum, and the transport
kernel then gives the retention of each electron. That is both cheaper and far
more precise than differencing two Monte-Carlo efficiencies, where a 5 % light
loss turns into a few tenths of a percent in eff_D and drowns.

Why this is worth running: v2.20.29 made the correction self-activate for
nuclides whose electron-capture probability reaches micCorrECThreshold = 0.01.
That threshold was chosen by an ARGUMENT -- P(EC) x 5 % of the light is about
0.05 %, under the Monte-Carlo noise of a normal run -- not by measurement. This
script measures the confinement loss of every nuclide against its P(EC), so
the threshold can be checked rather than assumed, and so any nuclide that
needs the correction for some other reason shows up.

Conditions follow micelleCheck.ipynb exactly: Ultima Gold + 10 % water, 4 nm
monodisperse droplets, kB = 1e-5 cm/keV, V = 10 mL.

Rows are appended to a JSON-lines file as they are produced. A run over the
whole catalogue takes hours, and a late crash that lost every earlier nuclide
would be a real cost, not a theoretical one. Re-running resumes.

Usage
-----
    python micelleCatalogue.py                 # the whole catalogue
    python micelleCatalogue.py H-3 Fe-55       # just these
    MICELLE_N_REC=20000 python micelleCatalogue.py   # tighter statistics
"""
import contextlib
import io
import json
import math
import os
import sys
import time
import warnings
from importlib.metadata import version

import numpy as np
from numba import njit

warnings.filterwarnings('ignore')

import tdcrpy.TDCR_model_lib as tl
import tdcrpy.TDCRPy as tm

# ---- conditions, identical to micelleCheck.ipynb -----------------------
KB = 1e-5                 # cm/keV
V = 10                    # mL
FAQ = 0.10                # aqueous volume fraction
DIAM = 4                  # nm, micelle DIAMETER (the kernel uses DIAM/2)
SIGMA = 0.0               # monodisperse
COCKTAIL = 'Ultima Gold'

#: decays recorded per nuclide. micelleCheck used 20000 for six nuclides; the
#: screening trades precision for coverage and reports the uncertainty it gets.
N_REC = int(os.environ.get('MICELLE_N_REC', 4000))
#: cap on electrons transported per nuclide. The metric is energy-weighted over
#: a random subsample, so this bounds the cost for a nuclide with a long
#: conversion-electron cascade without biasing F.
MAX_ELECTRONS = int(os.environ.get('MICELLE_MAX_E', 40000))

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   'micelleCatalogue.jsonl')


@njit(fastmath=True)
def _count_kernel_vec(E0_keV, r_d_nm, sigma_m, f_w):
    """Retained energy and droplet crossings, one transport per energy.

    Body identical to micelleCheck.ipynb's _count_kernel, which is itself a
    line-by-line mirror of TDCR_model_lib._mc_kernel with a crossing counter
    added -- same transport, same sampling, same order of random draws. The
    only change here is that the outer loop runs over an ARRAY of energies
    instead of n_samples at one fixed energy, which is what makes a
    catalogue-wide run affordable. Validated against the library kernel at
    import; if the library changes, the check fails rather than the numbers
    quietly drifting.
    """
    n = E0_keV.shape[0]
    e_out = np.zeros(n)
    n_out = np.zeros(n)
    r_eff = tl.mean_truncated_radius(r_d_nm, sigma_m)
    lam = (4.0 / 3.0) * r_eff * (1.0 - f_w) / f_w

    for i in range(n):
        R = tl.calc_range_numba(E0_keV[i])
        if R <= 0.0:
            continue
        # hydrophilic tracer: the decay starts inside a droplet, and escaping
        # it is the first crossing
        current_r = tl.get_micelle_radius(r_d_nm, sigma_m)
        u1 = np.random.rand()
        r_pos = current_r * (u1 ** (1.0 / 3.0))
        cos_t = np.random.uniform(-1.0, 1.0)
        b = 2.0 * r_pos * cos_t
        c = (r_pos ** 2) - (current_r ** 2)
        delta = max((b ** 2) - (4.0 * c), 0.0)
        R -= (-b + np.sqrt(delta)) / 2.0
        n_cross = 1
        active = True

        e_act = 0.0
        while R > 0.0:
            u = np.random.uniform(1e-10, 1.0)
            if active:
                d = -lam * np.log(u)
            else:
                current_r = tl.get_micelle_radius(r_d_nm, sigma_m)
                d = 2.0 * current_r * np.sqrt(u)
                n_cross += 1
            step = min(d, R)
            if active:
                e_act += (tl.calc_energy_numba(R)
                          - tl.calc_energy_numba(R - step))
            R -= step
            active = not active
        e_out[i] = e_act
        n_out[i] = n_cross
    return e_out, n_out


def set_state():
    """Configure the cocktail and reload, with the correction OFF for real.

    micCorrAuto must be pinned off. Since v2.20.29 the correction
    self-activates for an electron-capture nuclide whenever fAq > 0, so at
    FAQ = 0.10 a plain modifyMicCorr(False) would leave Fe-55, Cr-51 and every
    other EC nuclide CORRECTED -- and this script would be recording the
    already-corrected electron spectrum and reporting the confinement of a
    cocktail against itself.
    """
    global tl, tm
    with contextlib.redirect_stdout(io.StringIO()):
        tl.modifyLScocktail(COCKTAIL, FAQ, 'Water', 0.0)
    tl.update_config_batch({
        'diam_micelle': DIAM, 'sigma_micelle': SIGMA,
        'micCorr': False, 'micCorrAuto': False,
        'delayedFluorescence': False,
    })
    import importlib
    tl = importlib.reload(tl)
    tm = importlib.reload(tm)
    assert not tl.micCorr and not tl.micCorrAuto, 'correction is not off'


def validate_kernel():
    """The mirrored kernel must still reproduce the library's retention."""
    r_d = DIAM / 2.0
    print('validating the mirrored kernel against the library:')
    ok = True
    for e in (0.2, 1.0, 5.7, 50.0):
        arr = np.full(60000, e)
        mine = _count_kernel_vec(arr, r_d, SIGMA, FAQ)[0].mean() / e
        theirs = tl.pure_mc_efficient_energy_numba(
            e, r_d_nm=r_d, sigma_r=SIGMA, f_w=FAQ, num_samples=60000)[0] / e
        agree = abs(mine - theirs) < 0.01
        ok &= agree
        print(f'  E = {e:5.2f} keV   mirror {mine:.4f}   library {theirs:.4f}'
              f'   {"ok" if agree else "MISMATCH"}')
    if not ok:
        raise AssertionError('the mirror no longer reproduces the library')
    print('  -> mirror validated; its crossing counts can be used\n')


def quenched_electrons(rad, n):
    """Birks-quenched electron energies (keV) over n decays, correction off."""
    tm.TDCRPy(2.0, rad, '1', n, KB, V, record=True)
    out = []
    with open(tl.record_file_path(2)) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            p = line.split()
            if len(p) == 5 and p[0] == '1':
                e = float(p[1]) * 1e-3
                if e > 0:
                    out.append(e)
    tl.cleanup_record_files()
    return np.asarray(out, dtype=float)


def ratio_uncertainty(num, den):
    """u(sum(num)/sum(den)) for a ratio of sums, WITH the covariance.

    num and den are paired per electron and strongly correlated -- a high-
    energy electron contributes to both -- so treating them as independent,
    or either as binomial, overstates the uncertainty. Same estimator as
    eff()'s u_L.
    """
    n = num.size
    if n < 2:
        return float('nan')
    mu_n, mu_d = num.mean(), den.mean()
    if mu_d == 0:
        return float('nan')
    v_n, v_d = num.var(ddof=1), den.var(ddof=1)
    cov = np.cov(num, den, ddof=1)[0, 1]
    r = mu_n / mu_d
    var = (r ** 2 / n) * (v_n / mu_n ** 2 + v_d / mu_d ** 2
                          - 2 * cov / (mu_n * mu_d))
    return math.sqrt(max(var, 0.0))


def study(rad, rng):
    """One nuclide: confinement, crossings, and the EC probability."""
    t0 = time.time()
    e_all = quenched_electrons(rad, N_REC)
    row = dict(rad=rad, n_decays=N_REC, n_electrons_total=int(e_all.size),
               p_ec=float(tl.ec_branch_probability(rad)))
    if e_all.size == 0:
        # A pure alpha or pure gamma emitter deposits no electron energy that
        # TDCRPy applies the correction to. Worth recording, not skipping:
        # physically the alpha light IS confined too, but TDCRPy's correction
        # reaches only electrons and positrons, so for these nuclides the
        # correction is inert by construction.
        row.update(status='no-electrons', F=None, conf=None, u_conf=None,
                   E_mean=None, n_used=0, secs=round(time.time() - t0, 1))
        return row

    e = e_all
    if e.size > MAX_ELECTRONS:
        e = rng.choice(e_all, MAX_ELECTRONS, replace=False)
    ret, nx = _count_kernel_vec(e, DIAM / 2.0, SIGMA, FAQ)

    f = float(ret.sum() / e.sum())
    u_f = ratio_uncertainty(ret, e)
    q = np.percentile(nx, [10, 50, 90])
    row.update(
        status='ok', n_used=int(e.size),
        E_mean=float(e.mean()), E_median=float(np.median(e)),
        F=f, u_F=float(u_f),
        dilution_pct=100.0 * FAQ,
        conf=100.0 * (1.0 - FAQ - f) / (1.0 - FAQ),
        u_conf=100.0 * float(u_f) / (1.0 - FAQ),
        n_mean=float(nx.mean()), n_sd=float(nx.std(ddof=1)),
        n_p10=float(q[0]), n_med=float(q[1]), n_p90=float(q[2]),
        secs=round(time.time() - t0, 1))
    return row


def main():
    key = dict(tdcrpy=version('TDCRPy'), N_REC=N_REC, FAQ=FAQ, DIAM=DIAM,
               SIGMA=SIGMA, KB=KB, V=V, COCKTAIL=COCKTAIL,
               MAX_ELECTRONS=MAX_ELECTRONS)
    print(f'TDCRPy {version("TDCRPy")} -- micelle confinement across the '
          f'catalogue')
    print(f'conditions: {COCKTAIL} + {FAQ:.0%} water, {DIAM} nm monodisperse, '
          f'kB = {KB:g} cm/keV, V = {V} mL')
    print(f'statistics: {N_REC} decays/nuclide, '
          f'<= {MAX_ELECTRONS} electrons transported\n')

    set_state()
    print(f'rho={tl.RHO:.4f} g/cm3, Z={tl.Z:.4f}, A={tl.A:.4f}; '
          f'micCorr={tl.micCorr}, micCorrAuto={tl.micCorrAuto}\n')
    validate_kernel()

    wanted = sys.argv[1:]
    if not wanted:
        wanted = sorted(n.split('.')[0] for n in tl.z_PenNuc.namelist()
                        if n.endswith('.PenNuc.txt'))

    # ---- resume: keep rows that match these exact conditions
    done = {}
    if os.path.exists(OUT):
        with open(OUT, encoding='utf-8') as fh:
            for line in fh:
                try:
                    r = json.loads(line)
                except ValueError:
                    continue
                if r.get('key') == key:
                    done[r['rad']] = r
        if done:
            print(f'resuming: {len(done)} nuclides already done under these '
                  f'conditions\n')

    todo = [r for r in wanted if r not in done]
    print(f'{len(wanted)} nuclides requested, {len(todo)} to run\n')
    hdr = (f'{"nuclide":10s} {"P(EC)":>7s} {"n_e":>7s} {"<E>/keV":>9s} '
           f'{"conf %":>14s} {"crossings med":>14s} {"s":>6s}')
    print(hdr)
    print('-' * len(hdr))

    rng = np.random.default_rng(12345)
    t_start = time.time()
    for i, rad in enumerate(todo, 1):
        try:
            row = study(rad, rng)
        except Exception as exc:
            row = dict(rad=rad, status=f'error: {type(exc).__name__}: {exc}',
                       p_ec=None, conf=None)
        row['key'] = key
        with open(OUT, 'a', encoding='utf-8') as fh:
            fh.write(json.dumps(row) + '\n')
        if row['status'] == 'ok':
            print(f'{rad:10s} {row["p_ec"]:7.4f} {row["n_used"]:7d} '
                  f'{row["E_mean"]:9.2f} '
                  f'{row["conf"]:8.3f} +/-{row["u_conf"]:5.3f} '
                  f'{row["n_med"]:14.1f} {row["secs"]:6.1f}', flush=True)
        else:
            print(f'{rad:10s} {"":7s} {"":7s} {"":9s} {row["status"]:>14s} '
                  f'{"":14s} {row.get("secs", 0):6.1f}', flush=True)
        if i % 20 == 0:
            el = time.time() - t_start
            print(f'   [{i}/{len(todo)} done, {el/60:.0f} min elapsed, '
                  f'~{el/i*(len(todo)-i)/60:.0f} min left]', flush=True)

    print(f'\ntotal {(time.time() - t_start)/60:.1f} min; rows in {OUT}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
