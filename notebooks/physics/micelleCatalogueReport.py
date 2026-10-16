# -*- coding: utf-8 -*-
"""Final table and gate test for the catalogue-wide micelle screening.

Merges the two passes:

  pass 1  micelleCatalogue.py          confinement of the ELECTRON light,
                                       all 230 nuclides
  pass 2  micelleCatalogueWeighted.py  light-WEIGHTED bias, f_e * conf_e,
                                       for every nuclide pass 1 put above the
                                       floor

Pass 2 is the quantity that predicts an efficiency bias and pass 1 is not.
Pass 1 divides by the electron light alone, so for an alpha emitter it
describes a component that carries almost none of the light: Po-209 came out
third worst in the catalogue at 6.7 %, and carries 0.38 % of its light in
electrons, so its actual bias is 0.026 %. Below the floor pass 1 is still a
valid UPPER BOUND on the bias, since f_e <= 1, which is why those nuclides
need no second pass.

What this answers: is micCorrECThreshold = 0.01 -- a gate on the
electron-capture probability, shipped in v2.20.29 -- the right criterion for
self-activating the reverse-micelle correction? That threshold was chosen from
an argument, not a measurement. Here it is measured.

Usage
-----
    python micelleCatalogueReport.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
P1 = os.path.join(HERE, 'micelleCatalogue.jsonl')
P2 = os.path.join(HERE, 'micelleCatalogueWeighted.jsonl')

#: bias (%) at or above which leaving the correction off is worth caring about
MATTERS = 0.05
#: the shipped gate
EC_THRESHOLD = 0.01
#: not a nuclide -- a test entry in the PenNuc library
NOT_A_NUCLIDE = ('Ja-1', 'AT211PO211')


def read(path):
    rows = {}
    if not os.path.exists(path):
        return rows, None
    key = None
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
    p1, key = read(P1)
    p2, _ = read(P2)
    if not p1:
        sys.exit(f'no pass-1 results: {P1}')

    merged = {}
    for rad, r in p1.items():
        if rad in NOT_A_NUCLIDE:
            continue
        if r.get('status') != 'ok':
            merged[rad] = dict(rad=rad, status=r.get('status'))
            continue
        w = p2.get(rad)
        if w is not None and w.get('status') == 'ok':
            merged[rad] = dict(rad=rad, status='ok', p_ec=r['p_ec'],
                               E_mean=r['E_mean'], n_med=r['n_med'],
                               f_e=w['f_e'], conf_e=w['conf_e'],
                               bias=w['bias'], u_bias=w['u_bias'],
                               weighted=True)
        else:
            # below the floor: pass 1 is an upper bound, since f_e <= 1
            merged[rad] = dict(rad=rad, status='ok', p_ec=r['p_ec'],
                               E_mean=r['E_mean'], n_med=r['n_med'],
                               f_e=None, conf_e=r['conf'],
                               bias=r['conf'], u_bias=r['u_conf'],
                               weighted=False)

    ok = {k: v for k, v in merged.items() if v.get('status') == 'ok'}
    none = {k: v for k, v in merged.items()
            if v.get('status') == 'no-electrons'}
    bad = {k: v for k, v in merged.items()
           if str(v.get('status', '')).startswith('error')}

    print('REVERSE-MICELLE CONFINEMENT ACROSS THE TDCRPy CATALOGUE')
    if key:
        print(f'TDCRPy {key.get("tdcrpy")}, {key.get("COCKTAIL")} + '
              f'{key.get("FAQ", 0) * 100:.0f} % water, {key.get("DIAM")} nm '
              f'monodisperse, kB = {key.get("KB"):g} cm/keV, '
              f'{key.get("N_REC")} decays/nuclide')
    print(f'{len(merged)} nuclides: {len(ok)} measured, {len(none)} with no '
          f'electron light, {len(bad)} failed')
    print(f'{sum(1 for v in ok.values() if v["weighted"])} carry a '
          f'light-weighted bias from pass 2; the rest sat below '
          f'{MATTERS} % in pass 1, which bounds the bias from above.\n')
    print(f'Bulk dilution is {key.get("FAQ", 0) * 100:.0f} % for every '
          f'nuclide and is degenerate with L, so a refitted light yield '
          f'absorbs it.\nThe bias column is the part it cannot: '
          f'f_e x confinement of the electron light.\n')

    order = sorted(ok.values(), key=lambda r: -r['bias'])
    hdr = (f'{"nuclide":10s} {"P(EC)":>7s} {"f_e":>7s} {"<E>/keV":>9s} '
           f'{"bias %":>17s} {"cross":>7s} {"gate":>6s}')
    print(hdr)
    print('-' * len(hdr))
    for r in order:
        if r['bias'] < MATTERS / 2:
            continue
        gate = 'on' if r['p_ec'] >= EC_THRESHOLD else 'OFF'
        fe = f'{r["f_e"]:7.4f}' if r['f_e'] is not None else '      -'
        print(f'{r["rad"]:10s} {r["p_ec"]:7.4f} {fe} {r["E_mean"]:9.2f} '
              f'{r["bias"]:10.4f} +/-{r["u_bias"]:5.4f} '
              f'{r["n_med"]:7.0f} {gate:>6s}')
    print(f'(nuclides below {MATTERS / 2} % omitted; '
          f'{sum(1 for r in ok.values() if r["bias"] < MATTERS / 2)} of them)')

    # ------------------------------------------------------------ the verdict
    print(f'\n\nDOES micCorrECThreshold = {EC_THRESHOLD} GATE THE RIGHT '
          f'NUCLIDES?')
    print(f'"matters" = light-weighted bias >= {MATTERS} %, resolved at '
          f'2 sigma.\n')

    def matters(r):
        return r['bias'] >= MATTERS and r['bias'] > 2 * r['u_bias']

    fn = sorted((r for r in ok.values()
                 if matters(r) and r['p_ec'] < EC_THRESHOLD),
                key=lambda r: -r['bias'])
    fp = sorted((r for r in ok.values()
                 if r['p_ec'] >= EC_THRESHOLD and not matters(r)),
                key=lambda r: -r['bias'])
    hit = [r for r in ok.values() if matters(r) and r['p_ec'] >= EC_THRESHOLD]

    print(f'correctly gated on : {len(hit)}')
    print(f'missed (false neg) : {len(fn)}   <-- silent bias')
    print(f'needless (false +) : {len(fp)}   <-- costs time only\n')

    print('MISSED -- bias matters, correction stays off:')
    for r in fn:
        print(f'  {r["rad"]:10s} bias {r["bias"]:7.4f} +/-{r["u_bias"]:.4f} % '
              f'<E> {r["E_mean"]:7.2f} keV  P(EC) {r["p_ec"]:.4g}')
    print('\nNEEDLESS -- gated on, bias negligible (worst 8):')
    for r in fp[:8]:
        print(f'  {r["rad"]:10s} bias {r["bias"]:7.4f} +/-{r["u_bias"]:.4f} % '
              f'<E> {r["E_mean"]:7.2f} keV  P(EC) {r["p_ec"]:.4g}')

    print('\nWHY THE GATE MISSES WHAT IT MISSES')
    print('Every missed nuclide is a SOFT-ELECTRON emitter with no electron-')
    print('capture branch: a low-endpoint beta (Pu-241 20.8 keV, H-3 18.6 keV,')
    print('Ra-228 45.5 keV, Pb-210 63 keV) or a low-energy isomeric')
    print('transition (Nb-93m, Tc-99m). The physical driver is a sub-keV')
    print('electron that cannot leave its birth droplet, and P(EC) is only a')
    print('proxy for that -- it catches the Auger route and misses the soft-')
    print('beta route entirely.')

    if none:
        print(f'\nNO ELECTRON LIGHT ({len(none)}) -- the correction cannot '
              f'act, whatever the flags say:')
        print('  ' + ' '.join(sorted(none)))
        print('  TDCRPy applies the correction to electrons and positrons '
              'only. The alpha light\n  of these is physically confined too, '
              'but the model does not reach it.')
    if bad:
        print(f'\nFAILED ({len(bad)}):')
        for k, v in sorted(bad.items()):
            print(f'  {k:10s} {v["status"]}')

    # ---------------------------------------------------------------- drivers
    arr = [(r['E_mean'], r['bias'], r['p_ec'], r['n_med'])
           for r in ok.values()]
    e, b, p, nm = (np.array(x) for x in zip(*arr))
    print('\nWHAT PREDICTS THE BIAS (Spearman rho):')
    for name, x in (('mean electron energy', e), ('P(EC)', p),
                    ('median crossings', nm)):
        rho = np.corrcoef(np.argsort(np.argsort(x)),
                          np.argsort(np.argsort(b)))[0, 1]
        print(f'  {name:22s} rho = {rho:+.3f}')

    # ------------------------------------------------------------------ plot
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.2))
        sc = axes[0].scatter(e, np.maximum(b, 1e-3), c=p, cmap='viridis',
                             s=28, edgecolor='k', linewidth=0.3,
                             vmin=0, vmax=1)
        axes[0].set_xscale('log')
        axes[0].set_yscale('log')
        axes[0].axhline(MATTERS, color='r', ls='--', lw=0.9,
                        label=f'{MATTERS} % "matters"')
        for r in sorted(ok.values(), key=lambda r: -r['bias'])[:6]:
            axes[0].annotate(r['rad'], (r['E_mean'], max(r['bias'], 1e-3)),
                             fontsize=7, xytext=(4, 2),
                             textcoords='offset points')
        axes[0].set_xlabel('mean quenched electron energy  (keV)')
        axes[0].set_ylabel('light-weighted confinement bias  (%)')
        axes[0].legend(fontsize=8)
        axes[0].grid(alpha=0.3, which='both')
        fig.colorbar(sc, ax=axes[0], label='P(EC)')

        gated = np.array([r['p_ec'] >= EC_THRESHOLD for r in ok.values()])
        axes[1].scatter(np.clip(p[gated], 1e-4, None),
                        np.maximum(b[gated], 1e-3), s=28, alpha=0.8,
                        edgecolor='k', linewidth=0.3, label='gate on')
        axes[1].scatter(np.clip(p[~gated], 1e-4, None),
                        np.maximum(b[~gated], 1e-3), s=28, alpha=0.8,
                        marker='s', color='crimson', edgecolor='k',
                        linewidth=0.3, label='gate off')
        axes[1].set_xscale('log')
        axes[1].set_yscale('log')
        axes[1].axvline(EC_THRESHOLD, color='b', ls='--', lw=0.9,
                        label=f'micCorrECThreshold = {EC_THRESHOLD}')
        axes[1].axhline(MATTERS, color='r', ls='--', lw=0.9)
        for r in fn[:5]:
            axes[1].annotate(r['rad'], (1e-4, r['bias']), fontsize=7,
                             xytext=(5, 0), textcoords='offset points')
        axes[1].set_xlabel('P(EC)   (clipped at 1e-4 for the log axis)')
        axes[1].set_ylabel('light-weighted confinement bias  (%)')
        axes[1].legend(fontsize=8, loc='lower right')
        axes[1].grid(alpha=0.3, which='both')
        axes[1].set_title('top-left quadrant = missed by the gate',
                          fontsize=9)

        fig.suptitle(f'Reverse-micelle confinement bias, {len(ok)} nuclides '
                     f'-- {key.get("COCKTAIL")} + '
                     f'{key.get("FAQ", 0) * 100:.0f} % water, '
                     f'{key.get("DIAM")} nm monodisperse', fontsize=11)
        fig.tight_layout()
        out = os.path.join(HERE, f'micelleCatalogue_{key.get("tdcrpy")}.png')
        fig.savefig(out, dpi=150)
        print(f'\nwrote {out}')
    except Exception as exc:
        print(f'\nplot skipped: {exc}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
