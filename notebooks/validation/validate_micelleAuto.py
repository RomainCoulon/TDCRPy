# -*- coding: utf-8 -*-
"""Validation of the self-activating micelle correction (v2.20.29).

The unit tests prove the decision logic and that the kernel is entered for the
right nuclide. This checks the consequence: with the SHIPPED configuration
(micCorr = False, micCorrAuto = True) in a micellar cocktail, does TDCRPy land
on the corrected efficiency for an electron-capture nuclide and on the
uncorrected one for a pure beta emitter?

Three arms per nuclide, all at the same L so that nothing is hidden by a
refitted free parameter:

    off      micCorrAuto = False, micCorr = False   -- pre-2.20.29 behaviour
    auto     micCorrAuto = True,  micCorr = False   -- what now ships
    on       micCorr = True                         -- explicit

The claim under test: auto == on for an EC nuclide, auto == off for a beta
emitter, and off != on wherever the correction matters. Each arm is an
independent Monte-Carlo sample, so the comparison is in units of the combined
standard uncertainty.

Result on v2.20.29, N = 20000, L = 6 keV-1, fAq = 0.10, 4 nm monodisperse:

    nuclide  P(EC)   eff_D off  eff_D auto  eff_D on  auto-on  auto-off  cost
    H-3     0.0000     0.69560     0.69584   0.67189   +6.9 s    +0.1 s  0.9x
    Fe-55   0.9999     0.80376     0.75972   0.75982   -0.0 s   -15.2 s  1.2x
    Cr-51   1.0000     0.75302     0.68677   0.68738   -0.2 s   -26.2 s  2.0x
    Ni-63   0.0000     0.85811     0.86034   0.85053   +3.4 s    +0.8 s  1.0x
    C-14    0.0000     0.96360     0.96191   0.96113   +0.5 s    -1.1 s  1.0x
    Tc-99   0.0000     0.97976     0.97798   0.97567   +1.8 s    -1.5 s  1.0x

Self-activation lands on the explicit arm for Fe-55 and Cr-51 (0.0 s and
0.2 s away) while sitting 15 s and 26 s below the uncorrected one, and leaves
every pure beta emitter on the uncorrected arm. The cost of the default is
1.2x for Fe-55 and 2.0x for Cr-51 and nothing at all elsewhere.

Note that the off-to-on gap is NOT the micelle effect: at fixed L it is
dominated by the (1 - fAq) bulk dilution, which is degenerate with L and so
absorbed once L is refitted -- which is why H-3 moves by 3.4 % here and by
0.46 % in micelleCheck.ipynb, where L is fitted to a measured TDCR. This
script tests which arm self-activation reproduces, not the size of the
physics.
"""
import json
import math
import os
import sys
import time
import warnings

from importlib.metadata import version

import numpy as np
import tdcrpy
import tdcrpy.TDCRPy as T

lib = tdcrpy.TDCR_model_lib

N = int(os.environ.get("VALIDATE_N", 20000))
L = 6.0          # keV-1, the working point used throughout the micelle study
KB = 1.0e-5      # cm/keV
FAQ = 0.10
V = 10           # mL

RADS = ("H-3", "Fe-55", "Cr-51", "Ni-63", "C-14", "Tc-99")


def arm(rad, micCorr, micCorrAuto):
    lib.update_config_batch({"fAq": FAQ, "micCorr": micCorr,
                             "micCorrAuto": micCorrAuto,
                             "micCorrECThreshold": 0.01})
    T._micelle_advised.clear()
    T._micelle_auto_announced.clear()
    # Warm the kernel up before timing: a numba compile charged to one arm is
    # how a 56x inflation factor was once invented out of nothing.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        T.TDCRPy(L, rad, "1", 5, KB, V)
        t0 = time.perf_counter()
        out = T.TDCRPy(L, rad, "1", N, KB, V)
        dt = time.perf_counter() - t0
    lib.cleanup_record_files()
    return out[2], out[3], dt          # eff_D, u(eff_D), seconds


def main():
    saved = (lib.fAq, lib.micCorr, lib.micCorrAuto, lib.micCorrECThreshold)
    rows = []
    try:
        print(f"TDCRPy {version('TDCRPy')} "
              f"validation of micCorrAuto -- N = {N}, L = {L} keV-1, "
              f"fAq = {FAQ}, diam = {lib.diam_micelle} nm, "
              f"sigma = {lib.sigma_micelle} nm")
        print()
        hdr = (f"{'nuclide':8s} {'P(EC)':>7s} {'eff_D off':>11s} "
               f"{'eff_D auto':>11s} {'eff_D on':>11s} {'auto-on':>9s} "
               f"{'auto-off':>9s} {'cost':>6s}")
        print(hdr)
        print("-" * len(hdr))
        for rad in RADS:
            p_ec = lib.ec_branch_probability(rad)
            d_off, u_off, t_off = arm(rad, False, False)
            d_auto, u_auto, t_auto = arm(rad, False, True)
            d_on, u_on, t_on = arm(rad, True, True)

            z_on = (d_auto - d_on) / math.hypot(u_auto, u_on)
            z_off = (d_auto - d_off) / math.hypot(u_auto, u_off)
            print(f"{rad:8s} {p_ec:7.4f} {d_off:11.5f} {d_auto:11.5f} "
                  f"{d_on:11.5f} {z_on:8.2f}s {z_off:8.2f}s "
                  f"{t_auto / t_off:5.1f}x")
            rows.append(dict(rad=rad, p_ec=p_ec,
                             off=d_off, u_off=u_off,
                             auto=d_auto, u_auto=u_auto,
                             on=d_on, u_on=u_on,
                             z_on=z_on, z_off=z_off,
                             t_off=t_off, t_auto=t_auto, t_on=t_on))
    finally:
        lib.update_config_batch({"fAq": saved[0], "micCorr": saved[1],
                                 "micCorrAuto": saved[2],
                                 "micCorrECThreshold": saved[3]})

    # ------------------------------------------------------------- verdict --
    print()
    ok = True
    for r in rows:
        critical = r["p_ec"] >= 0.01
        want = "on" if critical else "off"
        agrees_with_target = abs(r["z_on"] if critical else r["z_off"]) < 4.0
        if not agrees_with_target:
            ok = False
            print(f"  FAIL {r['rad']}: self-activated result is not the "
                  f"{want} arm (z = "
                  f"{r['z_on'] if critical else r['z_off']:.2f})")
        else:
            print(f"  ok   {r['rad']:7s} self-activated = {want} arm "
                  f"(z = {r['z_on'] if critical else r['z_off']:+.2f})")
        if critical and abs(r["z_off"]) < 2.0:
            print(f"       note {r['rad']}: off and auto agree too -- the "
                  f"correction is below the noise at N = {N}")

    here = os.path.dirname(os.path.abspath(__file__))
    stem = os.path.join(here, f"validation_micelleAuto_{version('TDCRPy')}")
    with open(stem + ".json", "w", encoding="utf-8") as fh:
        json.dump(dict(N=N, L=L, kB=KB, fAq=FAQ,
                       diam=lib.diam_micelle, sigma=lib.sigma_micelle,
                       tdcrpy=version("TDCRPy"), rows=rows), fh, indent=1)
    print(f"\nwrote {stem}.json")

    # ---------------------------------------------------------------- plot --
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        x = np.arange(len(rows))
        fig, ax = plt.subplots(figsize=(8.2, 4.4))
        for key, lab, mk in (("off", "micCorrAuto = False (pre-2.20.29)", "s"),
                             ("auto", "shipped default (self-activating)", "o"),
                             ("on", "micCorr = True (explicit)", "x")):
            ax.errorbar(x, [100 * r[key] for r in rows],
                        yerr=[100 * r["u_" + key] for r in rows]
                        if "u_" + key in rows[0] else None,
                        marker=mk, ls="none", ms=8, capsize=3, label=lab)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{r['rad']}\nP(EC)={r['p_ec']:.2f}" for r in rows])
        ax.set_ylabel(r"$\varepsilon_D$  (%)")
        ax.set_title(f"Self-activating micelle correction, L = {L} "
                     r"keV$^{-1}$, $f_{Aq}$ = " f"{FAQ}, N = {N}")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(stem + ".png", dpi=150)
        print(f"wrote {stem}.png")
    except Exception as exc:
        print(f"plot skipped: {exc}")

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
