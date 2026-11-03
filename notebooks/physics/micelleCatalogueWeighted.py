# -*- coding: utf-8 -*-
"""Second pass: light-WEIGHTED confinement, the quantity that predicts bias.

micelleCatalogue.py measures confinement against the ELECTRON light only --
its record filter keeps kpar 1 (electron/Auger) and divides by the sum of
those energies. For a beta emitter that is the whole scintillating signal and
the number is the bias. For an alpha emitter it is not: Po-209 came out at
6.12 % confinement, third worst in the catalogue, but a 4.9 MeV alpha carries
almost all of its light and TDCRPy's correction does not touch alphas at all.
Reporting 6.12 % as an efficiency bias would be wrong by orders of magnitude.

What actually predicts the bias is the share of the light that is both
(a) carried by particles the correction acts on and (b) confined:

    f_e      = electron+positron quenched light / TOTAL quenched light
    conf_e   = confinement of that component   (pass 1's number)
    bias     = f_e * conf_e                    -- light-weighted confinement

Two corrections to pass 1 besides the weighting:

  * positrons (kpar 3) are included. The correction applies to them in
    _quench, and pass 1's filter dropped them -- so for a beta-plus emitter
    pass 1 divided by too small a denominator and over-reported conf_e.
  * alphas (kpar 4) and photons (kpar 2) enter the TOTAL light. Photons
    deposit through the electrons they create, which are already recorded as
    kpar 1, so counting the photon line too would double-count; the quenched
    energy written for a photon is zero in _quench, so summing it is harmless
    and it is summed only for the record.

Run after micelleCatalogue.py. Only the nuclides whose pass-1 confinement
reaches the threshold are worth re-running, since below it the weighting can
only make the number smaller.

Usage
-----
    python micelleCatalogueWeighted.py            # everything >= 0.05 %
    python micelleCatalogueWeighted.py Po-209     # named nuclides
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

warnings.filterwarnings('ignore')

import tdcrpy.TDCR_model_lib as tl
import tdcrpy.TDCRPy as tm

from micelleCatalogue import (FAQ, DIAM, SIGMA, KB, V, COCKTAIL, N_REC,
                              _count_kernel_vec, set_state, ratio_uncertainty)

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'micelleCatalogue.jsonl')
OUT = os.path.join(HERE, 'micelleCatalogueWeighted.jsonl')

#: pass-1 confinement (%) at or above which the weighting is worth computing
FLOOR = float(os.environ.get('MICELLE_FLOOR', 0.05))


def recorded_light(rad, n):
    """Quenched energies by particle class over n decays, correction off.

    Returns
    -------
    scint : ndarray
        Electron, Auger and positron quenched energies (keV) -- the light the
        micelle correction acts on.
    total : float
        Sum of all quenched energies (keV), alphas included.
    """
    tm.TDCRPy(2.0, rad, '1', n, KB, V, record=True)
    scint, total = [], 0.0
    with open(tl.record_file_path(2)) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            p = line.split()
            if len(p) != 5:
                continue
            e = float(p[1]) * 1e-3
            if e <= 0:
                continue
            total += e
            if p[0] in ('1', '3'):        # electron/Auger, positron
                scint.append(e)
    tl.cleanup_record_files()
    return np.asarray(scint, dtype=float), total


def main():
    pass1 = {}
    if os.path.exists(SRC):
        with open(SRC, encoding='utf-8') as fh:
            for line in fh:
                try:
                    r = json.loads(line)
                except ValueError:
                    continue
                pass1[r['rad']] = r

    standalone = os.environ.get('MICELLE_ALL') == '1'
    wanted = sys.argv[1:]
    if not wanted:
        if standalone:
            # the whole library, independent of any earlier pass
            wanted = sorted(n.split('.')[0] for n in tl.z_PenNuc.namelist()
                            if n.endswith('.PenNuc.txt'))
        else:
            wanted = sorted(k for k, v in pass1.items()
                            if v.get('status') == 'ok'
                            and v.get('conf', 0) >= FLOOR)

    key = dict(tdcrpy=version('TDCRPy'), N_REC=N_REC, FAQ=FAQ, DIAM=DIAM,
               SIGMA=SIGMA, KB=KB, V=V, COCKTAIL=COCKTAIL, pass_=2)
    done = set()
    if os.path.exists(OUT):
        with open(OUT, encoding='utf-8') as fh:
            for line in fh:
                try:
                    r = json.loads(line)
                except ValueError:
                    continue
                if r.get('key') == key:
                    done.add(r['rad'])

    todo = [r for r in wanted if r not in done]
    set_state()
    print(f'TDCRPy {version("TDCRPy")} -- light-weighted confinement, pass 2')
    print(f'micCorr={tl.micCorr}, micCorrAuto={tl.micCorrAuto}; '
          f'{len(wanted)} nuclides above {FLOOR} %, {len(todo)} to run\n')
    hdr = (f'{"nuclide":10s} {"f_e":>7s} {"conf_e %":>16s} '
           f'{"bias %":>16s} {"s":>6s}')
    print(hdr)
    print('-' * len(hdr))

    t0 = time.time()
    for i, rad in enumerate(todo, 1):
        t1 = time.time()
        try:
            scint, total = recorded_light(rad, N_REC)
            if scint.size == 0 or total <= 0:
                # no electron or positron light at all: the correction reaches
                # nothing, so the bias is a measured zero rather than a gap
                row = dict(rad=rad, status='no-light', f_e=0.0, bias=0.0,
                           u_bias=0.0, conf_e=None, n_scint=0,
                           p_ec=float(tl.ec_branch_probability(rad)))
            else:
                ret, _ = _count_kernel_vec(scint, DIAM / 2.0, SIGMA, FAQ)
                f = float(ret.sum() / scint.sum())
                u_f = ratio_uncertainty(ret, scint)
                conf_e = 100.0 * (1.0 - FAQ - f) / (1.0 - FAQ)
                u_conf = 100.0 * float(u_f) / (1.0 - FAQ)
                f_e = float(scint.sum() / total)
                prev = pass1.get(rad, {})
                row = dict(rad=rad, status='ok', f_e=f_e,
                           scint_keV=float(scint.sum()), total_keV=total,
                           n_scint=int(scint.size), F=f,
                           E_mean=float(scint.mean()),
                           conf_e=conf_e, u_conf_e=u_conf,
                           bias=f_e * conf_e, u_bias=f_e * u_conf,
                           conf_pass1=prev.get('conf'),
                           p_ec=float(tl.ec_branch_probability(rad)))
        except Exception as exc:
            row = dict(rad=rad, status=f'error: {type(exc).__name__}: {exc}',
                       f_e=None, bias=None)
        row['secs'] = round(time.time() - t1, 1)
        row['key'] = key
        with open(OUT, 'a', encoding='utf-8') as fh:
            fh.write(json.dumps(row) + '\n')
        if row['status'] == 'ok':
            print(f'{rad:10s} {row["f_e"]:7.4f} '
                  f'{row["conf_e"]:10.3f} +/-{row["u_conf_e"]:5.3f} '
                  f'{row["bias"]:10.4f} +/-{row["u_bias"]:5.4f} '
                  f'{row["secs"]:6.1f}', flush=True)
        else:
            print(f'{rad:10s} {row["status"]:>40s} {row["secs"]:6.1f}',
                  flush=True)
        if i % 20 == 0:
            el = time.time() - t0
            print(f'   [{i}/{len(todo)}, {el/60:.0f} min, '
                  f'~{el/i*(len(todo)-i)/60:.0f} min left]', flush=True)

    print(f'\ntotal {(time.time() - t0)/60:.1f} min; rows in {OUT}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
