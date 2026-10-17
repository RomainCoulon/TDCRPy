# -*- coding: utf-8 -*-
"""Validation of the measured-bias gate (v2.20.31).

v2.20.29 self-activated the reverse-micelle correction for a nuclide whose
electron-capture probability reached micCorrECThreshold = 0.01. The
catalogue-wide screening (notebooks/physics/micelleCatalogue.py) measured the
quantity that proxy stood in for -- the light-weighted confinement loss -- and
showed the proxy was the weakest of three candidates:

    Spearman rho against the measured bias
      median droplet crossings   -0.821
      mean electron energy       -0.767
      P(EC)                      +0.482

v2.20.31 gates on the measured bias instead, at micCorrBiasThreshold = 0.01 %.
This script checks the swap did what it was supposed to, by replaying both
gates over the whole shipped table. It is pure arithmetic on data already
measured -- no Monte Carlo, so it runs in a second and is exact.

Usage
-----
    python validate_micelleBiasGate.py
"""
import os
import sys
import warnings

warnings.filterwarnings('ignore')

import tdcrpy.TDCR_model_lib as lib

NEW_THRESHOLD = 0.01     # per cent of light
OLD_THRESHOLD = 0.01     # probability, the v2.20.29 gate
SKIP = ('Ja-1', 'AT211PO211')


def main():
    table = lib._load_micelle_bias()
    if not table:
        sys.exit('the bias table did not load')
    rads = sorted(r for r in table if r not in SKIP)
    print(f'TDCRPy {lib.version("TDCRPy") if hasattr(lib, "version") else ""}'
          f' -- measured-bias gate validation over {len(rads)} nuclides')
    print(f'new gate: measured bias >= {NEW_THRESHOLD} % of the light')
    print(f'old gate: P(EC) >= {OLD_THRESHOLD}\n')

    rows = []
    for r in rads:
        bias = table[r][0]
        rows.append((r, bias, lib.ec_branch_probability(r)))

    new_on = {r for r, b, _ in rows if b >= NEW_THRESHOLD}
    old_on = {r for r, _, p in rows if p >= OLD_THRESHOLD}

    print(f'activates under the new gate : {len(new_on)}')
    print(f'activates under the old gate : {len(old_on)}\n')

    gained = sorted((r, b, p) for r, b, p in rows
                    if r in new_on and r not in old_on)
    dropped = sorted((r, b, p) for r, b, p in rows
                     if r in old_on and r not in new_on)

    ok = True

    # ---- 1. the new gate must be exactly "bias >= threshold" -------------
    for r, b, _ in rows:
        if (b >= NEW_THRESHOLD) != lib.is_micelle_critical(r):
            ok = False
            print(f'  FAIL {r}: bias {b} but is_micelle_critical '
                  f'{lib.is_micelle_critical(r)}')
    print(f'1. gate agrees with the table for all {len(rows)} nuclides: '
          f'{"ok" if ok else "FAILED"}')

    # ---- 2. every nuclide the old gate missed is now covered -------------
    print(f'\n2. NOW COVERED, missed before ({len(gained)}) -- '
          f'sorted by what was at stake:')
    for r, b, p in sorted(gained, key=lambda t: -t[1]):
        print(f'   {r:10s} bias {b:8.4f} %   P(EC) {p:.4g}')
    worst = max((b for _, b, _ in gained), default=0.0)
    print(f'   largest bias previously ungated: {worst:.4f} %')
    if worst < 0.5:
        ok = False
        print('   FAIL: expected Pu-241 (0.598 %) and H-3 (0.531 %) here')

    # ---- 3. what the new gate stops paying for ---------------------------
    print(f'\n3. NO LONGER ACTIVATED ({len(dropped)}) -- '
          f'EC branch but hard electrons:')
    for r, b, p in sorted(dropped, key=lambda t: -t[1])[:12]:
        print(f'   {r:10s} bias {b:8.4f} %   P(EC) {p:.4g}')
    if dropped:
        biggest = max(b for _, b, _ in dropped)
        print(f'   largest bias now skipped: {biggest:.4f} % '
              f'(threshold {NEW_THRESHOLD} %)')
        if biggest >= NEW_THRESHOLD:
            ok = False
            print('   FAIL: dropped a nuclide at or above the threshold')

    # ---- 4. no nuclide with a serious bias is left out -------------------
    serious = [(r, b) for r, b, _ in rows if b >= 0.05]
    missed = [(r, b) for r, b in serious if r not in new_on]
    print(f'\n4. nuclides with bias >= 0.05 %: {len(serious)}; '
          f'left ungated: {len(missed)}')
    if missed:
        ok = False
        for r, b in missed:
            print(f'   FAIL {r}: {b} % and not gated')

    # ---- 5. a pure alpha emitter must never be gated on -----------------
    alphas = [r for r, b, _ in rows if table[r][2] == 'no-electron-light']
    bad = [r for r in alphas if r in new_on]
    print(f'\n5. nuclides with no electron light: {len(alphas)}; '
          f'wrongly gated on: {len(bad)}')
    if bad:
        ok = False
        print('   FAIL ' + ' '.join(bad))

    # ---- 6. the eight above 1 % are the strong Auger emitters -----------
    top = sorted(((b, r) for r, b, _ in rows), reverse=True)[:8]
    print('\n6. the eight largest biases:')
    for b, r in top:
        print(f'   {r:10s} {b:7.3f} %')
    expect = {'Ar-37', 'Ca-41', 'Cr-51', 'Fe-55', 'Ni-59', 'Ge-68', 'Ti-44',
              'Co-57'}
    if {r for _, r in top} != expect:
        ok = False
        print(f'   FAIL: expected {sorted(expect)}')

    print(f'\n{"PASS" if ok else "FAIL"}')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
