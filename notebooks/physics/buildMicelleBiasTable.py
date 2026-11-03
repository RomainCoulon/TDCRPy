# -*- coding: utf-8 -*-
"""Build the shipped micelle-bias table from the catalogue screening.

Writes tdcrpy/Micelle/micelleBias.csv, the data v2.20.31's self-activation
gate reads instead of inferring from the electron-capture probability.

Why a table and not a formula: P(EC) explained less than half the ranking
(Spearman rho +0.482, against -0.821 for the droplet-crossing count and
-0.767 for the mean electron energy). It catches the Auger route to a slow
electron and misses the soft-beta route entirely -- Pu-241, H-3, Ra-228,
Pb-210, Nb-93m and Tc-99m all have a resolved bias with no EC branch at all.
The screening measured the real quantity for every nuclide in the catalogue,
so the gate can read it rather than correlate with it.

Columns
-------
nuclide      PenNuc label
bias_pct     light-weighted confinement loss, per cent of total light. This
             is the part a refitted L does NOT absorb; the (1 - fAq) bulk
             dilution, which it does, is excluded.
u_bias_pct   standard uncertainty on bias_pct
f_e          electron+positron share of the quenched light. Below the
             screening floor this is unmeasured and left empty, and bias_pct
             then carries pass 1's electron-only confinement, which bounds
             the bias from above since f_e <= 1.
conf_e_pct   confinement of the electron light alone
p_ec         electron-capture probability per decay, for reference only
E_mean_keV   mean quenched electron energy
n_cross_med  median droplet crossings per electron
source       "weighted" (pass 2) or "bound" (pass 1 upper bound)
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
P1 = os.path.join(HERE, 'micelleCatalogue.jsonl')
P2 = os.path.join(HERE, 'micelleCatalogueWeighted.jsonl')
OUT = os.path.join(ROOT, 'tdcrpy', 'Micelle', 'micelleBias.csv')

#: not nuclides -- entries in the PenNuc library that are not single decays
SKIP = ('Ja-1', 'AT211PO211')


def read(path):
    rows, key = {}, None
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            try:
                r = json.loads(line)
            except ValueError:
                continue
            rows[r['rad']] = r
            key = r.get('key', key)
    return rows, key


def main():
    p1, key1 = read(P1)
    p2, key2 = read(P2)

    # bias_pct comes from pass 2, so the header must quote pass 2's statistics
    # and version. Quoting pass 1's understated them by 5x once already. The
    # reference-only columns (p_ec, E_mean_keV, n_cross_med) still come from
    # pass 1, which the header says.
    phys = ('COCKTAIL', 'FAQ', 'DIAM', 'SIGMA', 'KB', 'V')
    differ = {k: (key1.get(k), key2.get(k)) for k in phys
              if key1.get(k) != key2.get(k)}
    if differ:
        sys.exit(f'the two passes were run under different conditions, so '
                 f'their columns cannot share a table: {differ}')
    key = dict(key2)

    out = []
    for rad, r in sorted(p1.items()):
        if rad in SKIP:
            continue
        if r.get('status') == 'no-electrons':
            # The correction reaches electrons and positrons only, so for
            # these it is inert whatever the gate decides. Recorded as a
            # measured zero, not omitted, so the gate never has to guess.
            out.append(dict(nuclide=rad, bias_pct='0.0', u_bias_pct='0.0',
                            f_e='0.0', conf_e_pct='', p_ec=f'{r["p_ec"]:.6g}',
                            E_mean_keV='', n_cross_med='',
                            source='no-electron-light'))
            continue
        if r.get('status') != 'ok':
            continue
        w = p2.get(rad)
        if w is not None and w.get('status') == 'ok':
            out.append(dict(
                nuclide=rad, bias_pct=f'{w["bias"]:.5g}',
                u_bias_pct=f'{w["u_bias"]:.4g}', f_e=f'{w["f_e"]:.5g}',
                conf_e_pct=f'{w["conf_e"]:.5g}', p_ec=f'{r["p_ec"]:.6g}',
                E_mean_keV=f'{r["E_mean"]:.4g}',
                n_cross_med=f'{r["n_med"]:.0f}', source='weighted'))
        else:
            out.append(dict(
                nuclide=rad, bias_pct=f'{r["conf"]:.5g}',
                u_bias_pct=f'{r["u_conf"]:.4g}', f_e='',
                conf_e_pct=f'{r["conf"]:.5g}', p_ec=f'{r["p_ec"]:.6g}',
                E_mean_keV=f'{r["E_mean"]:.4g}',
                n_cross_med=f'{r["n_med"]:.0f}', source='bound'))

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fields = ['nuclide', 'bias_pct', 'u_bias_pct', 'f_e', 'conf_e_pct',
              'p_ec', 'E_mean_keV', 'n_cross_med', 'source']
    with open(OUT, 'w', encoding='utf-8', newline='') as fh:
        fh.write('# Reverse-micelle confinement bias per nuclide, measured by\n'
                 '# notebooks/physics/micelleCatalogue.py (pass 1) and\n'
                 '# micelleCatalogueWeighted.py (pass 2).\n'
                 '#\n'
                 '# REFERENCE CONDITIONS -- the table is a screening at these,\n'
                 '# not a correction. The correction itself always uses the\n'
                 f'# configured values.\n'
                 f'#   cocktail   {key.get("COCKTAIL")}\n'
                 f'#   fAq        {key.get("FAQ")}\n'
                 f'#   diameter   {key.get("DIAM")} nm\n'
                 f'#   sigma      {key.get("SIGMA")} nm\n'
                 f'#   kB         {key.get("KB")} cm/keV\n'
                 f'#   V          {key.get("V")} mL\n'
                 f'#   N          {key.get("N_REC")} decays/nuclide\n'
                 f'#   TDCRPy     {key.get("tdcrpy")}\n'
                 '#\n'
                 '# bias_pct, u_bias_pct, f_e and conf_e_pct are measured at\n'
                 f'# the N above. The reference-only columns p_ec, E_mean_keV\n'
                 f'# and n_cross_med come from pass 1, at N = '
                 f'{key1.get("N_REC")} and TDCRPy {key1.get("tdcrpy")}.\n'
                 '#\n'
                 '# bias_pct is the light-weighted confinement loss: the part\n'
                 '# a refitted L does NOT absorb. The (1 - fAq) bulk dilution,\n'
                 '# which it does absorb, is excluded.\n'
                 '#\n'
                 '# A larger aqueous fraction shortens the organic mean free\n'
                 '# path and raises the bias, so at fAq > the reference these\n'
                 '# values UNDERSTATE it.\n')
        wr = csv.DictWriter(fh, fieldnames=fields)
        wr.writeheader()
        for row in out:
            wr.writerow(row)

    n_w = sum(1 for r in out if r['source'] == 'weighted')
    n_b = sum(1 for r in out if r['source'] == 'bound')
    n_z = sum(1 for r in out if r['source'] == 'no-electron-light')
    print(f'wrote {OUT}')
    print(f'  {len(out)} nuclides: {n_w} weighted, {n_b} bounded, '
          f'{n_z} with no electron light')
    for thr in (0.01, 0.05, 0.1, 1.0):
        n = sum(1 for r in out if float(r['bias_pct']) >= thr)
        print(f'  gate at {thr:>4} %: {n} nuclides activate')
    return 0


if __name__ == '__main__':
    sys.exit(main())
