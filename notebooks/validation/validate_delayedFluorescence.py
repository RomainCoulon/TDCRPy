# -*- coding: utf-8 -*-
"""Validation of the prompt + delayed (TTA) fluorescence option (v2.20.30).

The unit tests pin the closed forms and reproduce the reference's published
Psi averages. This checks the one thing the model exists to produce: a light
yield, and hence an efficiency, that depends on the coincidence resolving
time W. The plain Birks model has no W-dependence at all.

PAIRED MEASUREMENT. The effect at the shipped parameters is a few tenths of a
percent in eff_D, which is at or below the Monte-Carlo noise of any affordable
N -- a first version of this script compared independent samples at N = 5000
and could not resolve it, reporting a spurious +1.8 % drift in the arm that
has none by construction. So every arm here is run from the SAME seed, which
gives every arm the same decay sample: the delayed model adds no random draws
and changing tau does not change the draw sequence, so a difference between
two arms at one seed is the model, exactly, with no MC noise in it at all.
That also makes the off arm's W-independence an exact check rather than a
statistical one.

Three arms per nuclide at fixed L:

    off        delayedFluorescence = False          -- pre-2.20.30
    UG         eta = 0.03, tau_d =  85 ns           -- Ultima Gold family
    toluene    eta = 0.05, tau_d = 138 ns           -- neat toluene

Fixed L is deliberate. Under a light yield refitted independently at each W
the prompt factor f_p is unidentifiable and the delayed factor is largely
absorbed too, since both are multiplicative rescalings -- which is why the
reference pipeline fits eta and tau_d against the residual trend in
Activity(W) rather than against a single efficiency.

Result on v2.20.30, N = 2000, 2 seeds, L = 6 keV-1, kB = 1e-5 cm/keV, Ultima
Gold stopping power. eff_D gain over the plain Birks model, per cent:

    nuclide  arm        W=10ns   30ns    100ns   300ns   1us
    H-3      UG        -1.0337 +0.4636 +0.9551 +1.3663 +1.6086
             toluene   -1.0199 +0.5285 +1.2310 +1.9785 +2.5124
    Fe-55    UG        -0.9514 +0.3911 +0.7968 +1.1291 +1.3219
             toluene   -0.9391 +0.4452 +1.0205 +1.6118 +2.0210
    Ni-63    UG        -0.3472 +0.1561 +0.3216 +0.4601 +0.5418
             toluene   -0.3426 +0.1780 +0.4146 +0.6665 +0.8465
    S-35     UG        -0.1128 +0.0495 +0.1018 +0.1456 +0.1713
             toluene   -0.1113 +0.0564 +0.1312 +0.2106 +0.2672
    C-14     UG        -0.0824 +0.0359 +0.0739 +0.1056 +0.1243
             toluene   -0.0813 +0.0409 +0.0952 +0.1527 +0.1937
    Tc-99    UG        -0.0461 +0.0206 +0.0424 +0.0607 +0.0714
             toluene   -0.0455 +0.0235 +0.0547 +0.0879 +0.1115

The off arm came out flat to 0.00e+00 -- bit-identical across all five
windows at a fixed seed -- for every one of the six, so none of them has an
isomeric cascade in the range where tau's other role would confound this.

Three things the numbers say:

1. The gain is NEGATIVE at W = 10 ns, for every nuclide. f_p(10, 3) = 0.964,
   so 3.6 % of the prompt light falls outside the window, and that outweighs
   the 11 % of the delayed reservoir the window has collected by then. The
   model is not a light BONUS; it redistributes light in time.
2. The ordering is H-3 > Fe-55 > Ni-63 > S-35 ~ C-14 > Tc-99 -- lowest mean
   energy first. Two causes compound: a soft spectrum draws more on the
   high-LET part of Psi, and a nuclide low on the efficiency curve converts a
   light change into an efficiency change, whereas Tc-99 at eff_D = 0.975 is
   saturated. This reproduces the reference's own finding that the W-trend it
   removes is largest for H-3.
3. Toluene (eta 0.05, tau_d 138 ns) gains more than Ultima Gold (0.03, 85 ns)
   at every window beyond the crossover, by about 1.5x at 1 us -- as it must,
   with the larger amplitude and the slower capture constant.

Reference: R. Coulon, "TDCR Measurements -- Lecture Notes", ch. 3, and
newTDCRmodel/DELAYED_FLUORESCENCE_MODEL.md.
"""
import json
import os
import statistics
import sys
import time
import warnings
from importlib.metadata import version

import numpy as np
import tdcrpy
import tdcrpy.TDCRPy as T

lib = tdcrpy.TDCR_model_lib

N = int(os.environ.get("VALIDATE_N", 4000))
SEEDS = tuple(int(x) for x in os.environ.get("VALIDATE_SEEDS", "1,2,3").split(","))
L = 6.0            # keV-1
KB = 1.0e-5        # cm/keV
V = 10             # mL
TAUS = (10.0, 30.0, 100.0, 300.0, 1000.0)
RADS = ("H-3", "C-14", "S-35", "Tc-99", "Ni-63", "Fe-55")
ARMS = (("off", None, None),
        ("UG", 0.03, 85.0),
        ("toluene", 0.05, 138.0))


def eff_d(rad, tau, eta, tau_d, seed):
    """eff_D at one (nuclide, window, arm, seed)."""
    if eta is None:
        lib.update_config_batch({"delayedFluorescence": False, "tau": tau})
    else:
        lib.update_config_batch({"delayedFluorescence": True, "tau": tau,
                                 "eta_delayed": eta, "tau_delayed": tau_d,
                                 "tau_prompt": 3.0})
    np.random.seed(seed)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = T.TDCRPy(L, rad, "1", N, KB, V)
    lib.cleanup_record_files()
    return out[2]


def main():
    saved = {"delayedFluorescence": lib.delayedFluorescence,
             "tau": lib.tau, "eta_delayed": lib.eta_delayed,
             "tau_delayed": lib.tau_delayed, "tau_prompt": lib.tau_prompt,
             "micCorr": lib.micCorr, "fAq": lib.fAq}
    rows = {}
    ok = True
    try:
        lib.update_config_batch({"micCorr": False, "fAq": 0.0})
        print(f"TDCRPy {version('TDCRPy')} -- delayed-fluorescence validation")
        print(f"N = {N}, seeds = {SEEDS}, L = {L} keV-1, kB = {KB:g} cm/keV, "
              f"cocktail = {lib.lsCocktail()}, "
              f"S_a = {lib.sa_tta():.1f} keV/cm")
        print(f"windows: {', '.join(f'{t:g}' for t in TAUS)} ns")
        print("\nPaired: each arm shares its decay sample with the others at "
              "the same seed,\nso a difference between arms carries no "
              "Monte-Carlo noise.\n")
        t_start = time.perf_counter()
        for rad in RADS:
            rows[rad] = {}
            print(f"{rad}:")
            per_arm = {}
            for name, eta, tau_d in ARMS:
                # mean over seeds of eff_D at each window
                by_tau = []
                for tau in TAUS:
                    vals = [eff_d(rad, tau, eta, tau_d, s) for s in SEEDS]
                    by_tau.append(vals)
                per_arm[name] = by_tau
                means = [statistics.fmean(v) for v in by_tau]
                rows[rad][name] = dict(eff=by_tau, mean=means)
                print(f"   {name:8s} eff_D " +
                      " ".join(f"{m:.6f}" for m in means))

            # ---- the off arm must be EXACTLY flat in W
            off = per_arm["off"]
            worst = 0.0
            for i_seed in range(len(SEEDS)):
                ref = off[0][i_seed]
                for i_tau in range(1, len(TAUS)):
                    worst = max(worst, abs(off[i_tau][i_seed] - ref))
            # The Birks light yield has no W-dependence, but tau has a
            # SECOND role in TDCRPy: it splits an isomeric transition slower
            # than tau into its own event, and that branch consumes random
            # draws. So for a nuclide with a cascade in that range the off
            # arm legitimately moves, and the two effects are entangled.
            rows[rad]["off"]["w_spread"] = worst
            note = ("exactly flat" if worst <= 1e-12 else
                    "moves -- tau also splits isomeric transitions here")
            print(f"   off arm spread over W, same seed: {worst:.2e}  ({note})")

            # ---- the paired gain, arm minus off, at each W
            for name in ("UG", "toluene"):
                gains = []
                for i_tau in range(len(TAUS)):
                    g = [100.0 * (per_arm[name][i_tau][k]
                                  / off[i_tau][k] - 1.0)
                         for k in range(len(SEEDS))]
                    gains.append(statistics.fmean(g))
                rows[rad][name]["gain_pct"] = gains
                print(f"   {name:8s} gain vs off (%) " +
                      " ".join(f"{g:+.4f}" for g in gains))
                if not all(b > a for a, b in zip(gains, gains[1:])):
                    ok = False
                    print(f"   FAIL: {name} gain is not monotone in W: "
                          f"{gains}")
                if gains[0] >= 0.0:
                    ok = False
                    print(f"   FAIL: {name} at W = {TAUS[0]:g} ns should LOSE "
                          f"light (f_p < 1), got {gains[0]:+.4f} %")
                # the gain has to stand clear of whatever the off arm did
                span = abs(gains[-1] - gains[0])
                if worst > 1e-12:
                    off_pct = 100.0 * worst / off[0][0]
                    if off_pct > 0.1 * span:
                        ok = False
                        print(f"   FAIL: {name}: the off arm's own W-drift "
                              f"({off_pct:.4f} %) is not negligible against "
                              f"the gain span ({span:.4f} %) -- the "
                              f"measurement is confounded")
            u_tol = rows[rad]["toluene"]["gain_pct"][-1]
            u_ug = rows[rad]["UG"]["gain_pct"][-1]
            if u_tol <= u_ug:
                ok = False
                print(f"   FAIL: toluene (eta 0.05) gains less than Ultima "
                      f"Gold (eta 0.03) at 1 us: {u_tol:+.4f} vs "
                      f"{u_ug:+.4f} %")
        print(f"\ntotal {time.perf_counter() - t_start:.0f} s")
    finally:
        lib.update_config_batch(saved)
        lib.cleanup_record_files()

    here = os.path.dirname(os.path.abspath(__file__))
    stem = os.path.join(here,
                        f"validation_delayedFluorescence_{version('TDCRPy')}")
    with open(stem + ".json", "w", encoding="utf-8") as fh:
        json.dump(dict(N=N, seeds=list(SEEDS), L=L, kB=KB, taus=list(TAUS),
                       tdcrpy=version("TDCRPy"), S_a=lib.sa_tta(),
                       rows=rows), fh, indent=1)
    print(f"wrote {stem}.json")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 3, figsize=(12.5, 6.8), sharex=True)
        for ax, rad in zip(axes.ravel(), RADS):
            for name in ("UG", "toluene"):
                ax.plot(TAUS, rows[rad][name]["gain_pct"],
                        marker="o", ms=4, label=name)
            ax.axhline(0.0, color="k", lw=0.6)
            ax.set_xscale("log")
            ax.set_title(rad, fontsize=10)
            ax.grid(alpha=0.3)
        for ax in axes[:, 0]:
            ax.set_ylabel(r"$\varepsilon_D$ gain over plain Birks  (%)")
        for ax in axes[1, :]:
            ax.set_xlabel("coincidence resolving time W  (ns)")
        axes[0, 0].legend(fontsize=8)
        fig.suptitle(f"Prompt + delayed (TTA) fluorescence, paired against "
                     f"the same decay sample\n"
                     f"L = {L} " r"keV$^{-1}$, " f"N = {N}, "
                     f"{len(SEEDS)} seeds, {lib.lsCocktail()} stopping power",
                     fontsize=11)
        fig.tight_layout()
        fig.savefig(stem + ".png", dpi=150)
        print(f"wrote {stem}.png")
    except Exception as exc:
        print(f"plot skipped: {exc}")

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
