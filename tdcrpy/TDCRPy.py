# -*- coding: utf-8 -*-
"""
TDCRPy — Monte Carlo code for TDCR measurement efficiency estimation.

This module exposes four public functions:

* :func:`TDCRPy`   — core stochastic/analytical efficiency calculation.
* :func:`eff`      — stochastic efficiency from a measured TDCR ratio.
* :func:`effA`     — analytical efficiency from a measured TDCR ratio.
* :func:`objectFct`— objective function used internally by :func:`eff`.

The internal helper :func:`_relax_atom` handles atomic-shell relaxation
cascades and is not intended for direct use.

Authors
-------
Romain Coulon, Jialin Hu
Bureau International des Poids et Mesures (BIPM)

References
----------
Coulon et al., *Applied Radiation and Isotopes* (2024),
https://doi.org/10.1016/j.apradiso.2024.111518
"""

# ---------------------------------------------------------------------------
# Standard-library and third-party imports
# ---------------------------------------------------------------------------
import os
import tempfile
import configparser
import importlib.resources
from importlib.resources import files

import numpy as np
from tqdm import tqdm
import scipy.optimize as opt

import tdcrpy.TDCR_model_lib as tl

# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

#: Radionuclides for which the fast analytical model is available.
_PURE_BETA_NUCLIDES = [
    "H-3", "C-14", "S-35", "Ca-45", "Ni-63",
    "Sr-89", "Sr-90", "Tc-99", "Pm-147", "Pu-241",
]

#: Matching discretisation counts (energy bins) for each pure-beta nuclide.
_PURE_BETA_NE = [7000, 1000, 1000, 500, 2000, 500, 200, 500, 1000, 7000]


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------
def _relax_atom(daughter_relax, particle_vec, energy_vec, rad,
                display=False, unc_data=False):
    """
    Propagate atomic-shell relaxation cascades in-place.
    Now correctly processes all secondary vacancies by dynamically checking list length.
    """
    i_part = 0
    while i_part < len(particle_vec):
        shell = particle_vec[i_part]
        
        # Only process shell-vacancy tokens.
        if not ("Atom_K" in shell or "Atom_L" in shell or "Atom_M" in shell):
            i_part += 1
            continue

        relaxation = True
        while relaxation:
            tf, ef = tl.relaxation_atom(
                daughter_relax, rad, particle_vec[i_part], uncData=unc_data
            )

            if tf == "XKA":
                # K→L transition: vacancy moves to L shell.
                particle_vec[i_part] = "Atom_L"
                particle_vec.append(tf)
                energy_vec.append(ef)

            elif tf == "XKB":
                # K→M transition: vacancy moves to M shell.
                particle_vec[i_part] = "Atom_M"
                particle_vec.append(tf)
                energy_vec.append(ef)
                relaxation = False

            elif tf == "XL":
                particle_vec[i_part] = "Atom_M"
                particle_vec.append(tf)
                energy_vec.append(ef)
                relaxation = False

            elif tf == "Auger KLL":
                # Two L-shell vacancies remain after KLL Auger.
                particle_vec[i_part] = "Atom_L"
                tf1, ef1 = tl.relaxation_atom(
                    daughter_relax, rad, "Atom_L", uncData=unc_data
                )
                particle_vec.append(tf)
                energy_vec.append(ef)
                particle_vec.append(tf1)
                energy_vec.append(ef1)
                if tf1 == "Auger L":
                    particle_vec.extend(["Atom_M", "Atom_M"])
                    energy_vec.extend([0, 0])
                else:
                    particle_vec.append("Atom_M")
                    energy_vec.append(0)
                # Loop again: the new Atom_L vacancy still needs relaxing.

            elif tf == "Auger KLX":
                particle_vec[i_part] = "Atom_L"
                particle_vec.extend(["Atom_M", tf])
                energy_vec.extend([0, ef])

            elif tf == "Auger KXY":
                particle_vec[i_part] = "Atom_M"
                particle_vec.extend(["Atom_M", tf])
                energy_vec.extend([0, ef])
                relaxation = False
                
            elif tf == "Auger L":
                particle_vec[i_part] = "Atom_M"
                particle_vec.extend(["Atom_M", tf])
                energy_vec.extend([0, ef])
                relaxation = False

            elif tf == "Auger M" or tf == "XM":
                # Relaxation concludes
                particle_vec[i_part] = "Atom_N"
                particle_vec.append(tf)
                energy_vec.append(ef)
                relaxation = False

            else:
                # Unresolved transition edge-case
                particle_vec[i_part] = "Atom_N"
                particle_vec.append(tf)
                energy_vec.append(ef)
                relaxation = False
                if display:
                    print("Unresolved transition:", tf)
                    
        i_part += 1

    return particle_vec, energy_vec


# def _relax_atom(daughter_relax, particle_vec, energy_vec, rad,
#                 display=False, unc_data=False):
#     """
#     Propagate atomic-shell relaxation cascades in-place.

#     After a vacancy is created (electron capture or internal conversion),
#     this function iteratively samples X-ray emissions and Auger electrons
#     until the atom has fully de-excited, updating *particle_vec* and
#     *energy_vec* accordingly.

#     Parameters
#     ----------
#     daughter_relax : str
#         Symbol of the daughter nucleus (e.g. ``"Co-60"``).
#     particle_vec : list of str
#         Mutable list of particle labels for the current decay event.
#         Shell-vacancy tokens (``"Atom_K"``, ``"Atom_L"``, ``"Atom_M"``)
#         are replaced or removed as relaxation proceeds.
#     energy_vec : list of float
#         Mutable list of particle energies (keV) paired with *particle_vec*.
#     rad : str
#         Parent radionuclide label, forwarded to
#         :func:`tdcrpy.TDCR_model_lib.relaxation_atom`.
#     display : bool, optional
#         If ``True``, print unresolved transition labels to stdout.
#         Default is ``False``.
#     unc_data : bool, optional
#         If ``True``, propagate nuclear-data uncertainties through
#         :func:`tdcrpy.TDCR_model_lib.relaxation_atom`. Default is ``False``.

#     Returns
#     -------
#     particle_vec : list of str
#         Updated particle-label list (same object as input).
#     energy_vec : list of float
#         Updated energy list (same object as input).

#     Notes
#     -----
#     The function modifies *particle_vec* and *energy_vec* **in place** and
#     also returns them for convenience.  Entries marked ``"Atom_K"``,
#     ``"Atom_L"``, or ``"Atom_M"`` signal shell vacancies; entries without
#     an ``"Atom_"`` prefix are physical particles.
#     """
#     for i_part in range(len(particle_vec)):
#         shell = particle_vec[i_part]
#         # Only process shell-vacancy tokens.
#         if not ("Atom_K" in shell or "Atom_L" in shell or "Atom_M" in shell):
#             continue

#         relaxation = True
#         while relaxation:
#             tf, ef = tl.relaxation_atom(
#                 daughter_relax, rad, particle_vec[i_part], uncData=unc_data
#             )

#             if tf == "XKA":
#                 # K→L transition: vacancy moves to L shell.
#                 particle_vec[i_part] = "Atom_L"
#                 particle_vec.append(tf)
#                 energy_vec.append(ef)

#             elif tf == "XKB":
#                 # K→M transition: vacancy moves to M shell.
#                 particle_vec[i_part] = "Atom_M"
#                 particle_vec.append(tf)
#                 energy_vec.append(ef)
#                 relaxation = False

#             elif tf == "XL":
#                 particle_vec[i_part] = "Atom_M"
#                 particle_vec.append(tf)
#                 energy_vec.append(ef)
#                 relaxation = False

#             elif tf == "Auger KLL":
#                 # Two L-shell vacancies remain after KLL Auger.
#                 particle_vec[i_part] = "Atom_L"
#                 tf1, ef1 = tl.relaxation_atom(
#                     daughter_relax, rad, "Atom_L", uncData=unc_data
#                 )
#                 particle_vec.append(tf)
#                 energy_vec.append(ef)
#                 particle_vec.append(tf1)
#                 energy_vec.append(ef1)
#                 if tf1 == "Auger L":
#                     particle_vec.extend(["Atom_M", "Atom_M"])
#                     energy_vec.extend([0, 0])
#                 else:
#                     particle_vec.append("Atom_M")
#                     energy_vec.append(0)
#                 # Loop again: the new Atom_L vacancy still needs relaxing.

#             elif tf == "Auger KLX":
#                 particle_vec[i_part] = "Atom_L"
#                 particle_vec.extend(["Atom_M", tf])
#                 energy_vec.extend([0, ef])
#                 # Loop again.

#             elif tf == "Auger KXY":
#                 particle_vec[i_part] = "Atom_M"
#                 particle_vec.extend(["Atom_M", tf])
#                 energy_vec.extend([0, ef])
#                 relaxation = False

#             elif tf == "Auger L":
#                 particle_vec[i_part] = "Atom_M"
#                 particle_vec.extend(["Atom_M", tf])
#                 energy_vec.extend([0, ef])
#                 relaxation = False

#             else:
#                 if display:
#                     print(f"\t\t undetermined X or Auger type = {tf}")
#                 relaxation = False

#     return particle_vec, energy_vec


def _read_config():
    """
    Load the TDCRPy configuration file and return a
    :class:`configparser.ConfigParser` instance together with the most
    commonly needed scalar parameters.

    Returns
    -------
    config : configparser.ConfigParser
    tau : float
        Coincidence resolving time (ns).
    ext_dt : float
        Extended dead time (µs).
    meas_time : float
        Measurement time (min).
    mic_corr : bool
        Whether the micelle correction is active.
    ne_electron : int
        Number of energy bins for electron quenching.
    ne_alpha : int
        Number of energy bins for alpha quenching.
    """
    config = configparser.ConfigParser()
    with importlib.resources.as_file(
        files("tdcrpy").joinpath("config.toml")
    ) as cfg_path:
        config.read(cfg_path)
    inp = config["Inputs"]
    return (
        config,
        inp.getfloat("tau"),
        inp.getfloat("extDT"),
        inp.getfloat("measTime"),
        inp.getboolean("micCorr"),
        inp.getint("nE_electron"),
        inp.getint("nE_alpha"),
    )


def _open_record_files(temp_dir):
    """
    Create and initialise the four temporary record files used when
    ``record=True``.

    Parameters
    ----------
    temp_dir : str
        Directory in which the temporary files are created.

    Returns
    -------
    recfile1, recfile2, recfile3, recfile4 : str
        Absolute paths to the four record files.
    """
    headers = {
        "Temp_E0.txt": (
            "# TDCRPy output: initial energies from nuclear decays\n"
            "# Column 1: KPAR (1=electron, 2=photon, 3=positron, 4=alpha)\n"
            "# Column 2: Energy in eV\n"
            "# Column 3: Decay number mod 100\n"
            "# Column 4: Cascade number mod 10\n"
            "# Column 5: Particle age (s)\n"
        ),
        "Temp_E1.txt": (
            "# TDCRPy output: deposited energies from nuclear decays\n"
            "# Column 1: KPAR (1=electron, 2=photon, 3=positron, 4=alpha)\n"
            "# Column 2: Energy in eV\n"
            "# Column 3: Decay number mod 100\n"
            "# Column 4: Cascade number mod 10\n"
            "# Column 5: Particle age (s)\n"
        ),
        "Temp_E2.txt": (
            "# TDCRPy output: quenched deposited energies from nuclear decays\n"
            "# Column 1: KPAR (1=electron, 2=photon, 3=positron, 4=alpha)\n"
            "# Column 2: Energy in eV\n"
            "# Column 3: Decay number mod 100\n"
            "# Column 4: Cascade number mod 10\n"
            "# Column 5: Particle age (s)\n"
        ),
        "Temp_E3.txt": (
            "# TDCRPy output: detection probabilities\n"
            "# Column 1: Decay number mod 100\n"
            "# Column 2: detection probability — single events\n"
            "# Column 3: detection probability — double coincidences\n"
            "# Column 4: detection probability — triple coincidences\n"
        ),
    }
    paths = []
    for fname, header in headers.items():
        p = os.path.join(temp_dir, fname)
        with open(p, "w") as f:
            f.write(header)
        paths.append(p)
    return tuple(paths)


def _particle_kpar(p):
    """
    Map a particle label to its KPAR integer code for record files.

    Returns
    -------
    kpar : str or None
        ``"1"`` (electron/Auger), ``"2"`` (photon/X-ray),
        ``"3"`` (positron), ``"4"`` (alpha), or ``None`` for unknown types.
    """
    if p == "electron" or "Auger" in p:
        return "1"
    if p == "gamma" or "X" in p:
        return "2"
    if p == "positron":
        return "3"
    if p == "alpha":
        return "4"
    return None


def _write_particles(file_handle, particle_vec, energy_vec, idec, t1=0):
    """
    Write one event's particles to an open record file.

    Parameters
    ----------
    file_handle : file-like
        Open file object (text mode, append).
    particle_vec : list of str
        Particle labels.
    energy_vec : list of float
        Energies in keV; written in eV.
    idec : int
        Decay index (used for formatting spacing).
    t1 : float, optional
        Particle age in seconds. Default is 0.
    """
    # Use a fixed-width decay index: always two characters when idec < 100.
    idx_str = f"{idec:2d}"
    for irec, p in enumerate(particle_vec):
        kpar = _particle_kpar(p)
        if kpar is None:
            continue
        file_handle.write(
            f"{kpar} {energy_vec[irec] * 1e3:.6E} {idx_str} 1 {t1:.6E}\n"
        )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def TDCRPy(
    L, Rad, pmf_1, N, kB, V,
    mode="eff",
    display=False,
    barp=False,
    Smodel=True,
    record=False,
    readRecHist=False,
    uncData=False,
    fullMC=False,
    # Legacy capitalisation aliases kept for backward compatibility.
    Display=None,
):
    """
    Core TDCRPy calculation: stochastic or analytical detection-efficiency
    estimation for a liquid scintillation counter.

    Each Monte Carlo trial simulates one complete nuclear decay and the
    subsequent scintillation / photon-detection chain, producing per-event
    detection probabilities that are averaged over *N* trials.

    Parameters
    ----------
    L : float or tuple of float
        Free parameter (keV⁻¹).  A scalar value implies a symmetric PMT
        configuration; a 3-tuple ``(L_A, L_B, L_C)`` enables asymmetric
        modelling of the three photomultiplier tubes.
    Rad : str
        Comma-separated list of radionuclide labels, e.g. ``"Co-60"`` or
        ``"Co-60, H-3"``.
    pmf_1 : str
        Comma-separated relative activity fractions matching *Rad*, e.g.
        ``"1"`` or ``"0.8, 0.2"``.  Values are normalised internally.
    N : int
        Number of Monte Carlo trials (= simulated decays).  The JCGM 101
        guide recommends N ≥ 10 000 for reliable uncertainty estimates.
        Ignored when the analytical model is selected via ``Smodel=False``.
    kB : float
        Birks constant (cm keV⁻¹), governing ionisation quenching.
    V : float
        Volume of the liquid scintillator (mL).
    mode : {"eff", "dis"}, optional
        ``"eff"`` (default) returns efficiency estimates and their
        uncertainties; ``"dis"`` returns the raw per-event efficiency lists.
    display : bool, optional
        Print step-by-step decay details to stdout.  Default is ``False``.
    barp : bool, optional
        Show a ``tqdm`` progress bar during the MC loop.  Default is
        ``False``.
    Smodel : bool, optional
        ``True`` (default) runs the full stochastic model.  ``False``
        selects the faster analytical model, which is only accurate for
        pure beta emitters listed in :data:`_PURE_BETA_NUCLIDES`.
    record : bool, optional
        Write intermediate decay histories to temporary files (useful for
        debugging and post-processing).  Default is ``False``.
    readRecHist : bool, optional
        Read previously recorded decay histories from the temporary files
        instead of running a new MC simulation.  Default is ``False``.
    uncData : bool, optional
        If ``True``, sample nuclear-data uncertainties for each trial
        (propagation of nuclear-data uncertainty).  Default is ``False``.
    fullMC : bool, optional
        If ``True``, use a fully stochastic photon-transport model
        (:func:`~tdcrpy.TDCR_model_lib.detectProbabilitiesMC`) instead of
        the semi-analytical default.  Default is ``False``.
    Display : bool or None, optional
        *Deprecated alias* for *display*.  If not ``None``, overrides
        *display*.

    Returns
    -------
    When ``mode="eff"``
        A 14-tuple of floats ``(mean, std)`` pairs for seven efficiency
        types in the following order:

        * ``mean_efficiency_S, std_efficiency_S`` — single events.
        * ``mean_efficiency_D, std_efficiency_D`` — logical sum of doubles.
        * ``mean_efficiency_T, std_efficiency_T`` — triple coincidences.
        * ``mean_efficiency_AB, std_efficiency_AB`` — AB doubles.
        * ``mean_efficiency_BC, std_efficiency_BC`` — BC doubles.
        * ``mean_efficiency_AC, std_efficiency_AC`` — AC doubles.
        * ``mean_efficiency_D2, std_efficiency_D2`` — C/N doubles.

    When ``mode="dis"``
        A 4-tuple of lists
        ``(efficiency_S, efficiency_D, efficiency_T, efficiency_D2)``,
        one entry per simulated decay.

    Raises
    ------
    ValueError
        If *pmf_1* does not sum to 1 (a warning is printed instead of
        raising to preserve backward compatibility).

    Notes
    -----
    The ``record`` → ``readRecHist`` two-step workflow can speed up
    repeated evaluations at different *L* values because the expensive
    nuclear-decay + quenching simulation is performed only once.

    See Also
    --------
    eff : Determine *L* and efficiencies from a measured TDCR ratio.
    effA : Analytical variant of :func:`eff`.
    """
    # ------------------------------------------------------------------ #
    # 0. Handle backward-compatible keyword alias                          #
    # ------------------------------------------------------------------ #
    if Display is not None:
        display = Display

    symm = not isinstance(L, (tuple, list))

    # ------------------------------------------------------------------ #
    # 1. Load configuration                                                #
    # ------------------------------------------------------------------ #
    config, tau, ext_dt, meas_time, mic_corr, ne_electron, ne_alpha = _read_config()

    # ------------------------------------------------------------------ #
    # 2. Analytical model short-circuit                                    #
    # ------------------------------------------------------------------ #
    if not Smodel:
        if Rad in _PURE_BETA_NUCLIDES:
            ne = _PURE_BETA_NE[_PURE_BETA_NUCLIDES.index(Rad)]
        else:
            ne = 1000
        out = tl.modelAnalytical(L, 1, 1, 1, 1, Rad, kB, V, mode, symm, ne)
        if mode == "eff":
            return out[0], 0, out[1], 0, out[2], 0
        return out

    # ------------------------------------------------------------------ #
    # 3. Record-file initialisation                                        #
    # ------------------------------------------------------------------ #
    if record:
        temp_dir = tempfile.gettempdir()
        recfile1, recfile2, recfile3, recfile4 = _open_record_files(temp_dir)

    # ------------------------------------------------------------------ #
    # 4. Optional header display                                           #
    # ------------------------------------------------------------------ #
    if barp:
        tl.display_header()

    # ------------------------------------------------------------------ #
    # 5. Read-from-history branch                                          #
    # ------------------------------------------------------------------ #
    if readRecHist:
        return _run_from_history(
            L, N, tau, ext_dt, meas_time, fullMC, mode
        )

    # ------------------------------------------------------------------ #
    # 6. Parse radionuclide list and fractions                             #
    # ------------------------------------------------------------------ #
    Rad = [r.strip() for r in Rad.replace(" ", "").split(",")]
    pmf_1 = [float(x) for x in pmf_1.split(",")]
    if len(pmf_1) > 1 and abs(sum(pmf_1) - 1.0) > 1e-9:
        print(f"Warning: pmf_1 sums to {sum(pmf_1):.6f}, expected 1.")

    # ------------------------------------------------------------------ #
    # 7. Load nuclear-decay data for all requested radionuclides           #
    # ------------------------------------------------------------------ #
    decay_data = _load_decay_data(Rad)
    (
        particle, p_branch, e_branch, LevelDaughter,
        levelNumber, prob_trans, prob_branch, levelEnergy,
        transitionType, e_trans, next_level, Q_value,
        DaughterVec, Pdaughter, Transition_prob_sum,
        u_prob_trans, trans_halfLife,
    ) = decay_data

    # ------------------------------------------------------------------ #
    # 8. Monte Carlo main loop                                             #
    # ------------------------------------------------------------------ #
    eff_lists = {
        "S": [], "D": [], "T": [],
        "AB": [], "BC": [], "AC": [], "D2": [],
    }

    iterator = (
        tqdm(range(N), desc="Processing", unit=" decays")
        if (barp and not display)
        else range(N)
    )

    for idec in iterator:
        particle_vec = []
        energy_vec = []

        # -------------------------------------------------------------- #
        # 8.0  Sample radionuclide from mixture                           #
        # -------------------------------------------------------------- #
        index_rad = tl.sampling(pmf_1)
        rad_i = Rad[index_rad]
        if display:
            print(f"\n\n Trial {idec + 1} — Sampled radionuclide: {rad_i}")

        # -------------------------------------------------------------- #
        # 8.1  Nuclear decay                                               #
        # -------------------------------------------------------------- #
        (
            particle_vec, energy_vec,
            particle_vec2, energy_vec2,
            evenement, t1,
            level_before_trans, iDaughter,
        ) = _sample_nuclear_decay(
            index_rad, Rad, DaughterVec, Pdaughter, particle, p_branch,
            e_branch, LevelDaughter, levelNumber, prob_trans, prob_branch,
            levelEnergy, transitionType, e_trans, next_level,
            Transition_prob_sum, trans_halfLife, u_prob_trans,
            tau, uncData, display,
        )

        # -------------------------------------------------------------- #
        # 8.2  Atomic relaxation                                           #
        # -------------------------------------------------------------- #
        daughter_relax = DaughterVec[index_rad][iDaughter]
        particle_vec, energy_vec = _relax_atom(
            daughter_relax, particle_vec, energy_vec,
            Rad[index_rad], display=display, unc_data=uncData,
        )
        if display:
            _print_relaxation("Prompt", particle_vec, energy_vec)

        if evenement != 1:
            particle_vec2, energy_vec2 = _relax_atom(
                daughter_relax, particle_vec2, energy_vec2,
                Rad[index_rad], display=display, unc_data=uncData,
            )
            if display:
                _print_relaxation("Delayed", particle_vec2, energy_vec2)

        # -------------------------------------------------------------- #
        # 8.3  Beta-particle spectrum sampling                             #
        # -------------------------------------------------------------- #
        particle_vec, energy_vec = _sample_beta_spectra(
            particle_vec, energy_vec, rad_i, level_before_trans, display
        )
        if evenement != 1:
            particle_vec2, energy_vec2 = _sample_beta_spectra(
                particle_vec2, energy_vec2, rad_i, level_before_trans, display
            )

        energy_vec_initial = energy_vec.copy()

        # -------------------------------------------------------------- #
        # 8.3a  Optional: record initial energies                         #
        # -------------------------------------------------------------- #
        if record:
            with open(recfile1, "a") as f:
                _write_particles(f, particle_vec, energy_vec, idec)
                if evenement != 1:
                    _write_particles(f, particle_vec2, energy_vec2, idec, t1)

        # -------------------------------------------------------------- #
        # 8.4  Radiation–matter interaction (energy deposition)           #
        # -------------------------------------------------------------- #
        particle_vec, energy_vec = _interact_prompt(particle_vec, energy_vec, V)
        if display:
            _print_interaction("Prompt", particle_vec, energy_vec)

        energy_vec_initial2 = None
        if evenement != 1:
            energy_vec_initial2 = energy_vec2.copy()
            particle_vec2, energy_vec2 = _interact_prompt(
                particle_vec2, energy_vec2, V
            )
            if display:
                _print_interaction("Delayed", particle_vec2, energy_vec2)

        # -------------------------------------------------------------- #
        # 8.4a  Optional: record deposited energies                       #
        # -------------------------------------------------------------- #
        if record:
            with open(recfile2, "a") as f:
                _write_particles(f, particle_vec, energy_vec, idec)
                if evenement != 1:
                    _write_particles(f, particle_vec2, energy_vec2, idec, t1)

        # -------------------------------------------------------------- #
        # 8.5  Scintillation quenching (Birks model)                      #
        # -------------------------------------------------------------- #
        e_quenching = _quench(
            particle_vec, energy_vec, energy_vec_initial,
            kB, ne_electron, ne_alpha, mic_corr, display,
        )

        e_quenching2 = 0
        if evenement != 1:
            e_quenching2 = _quench(
                particle_vec2, energy_vec2, energy_vec_initial2,
                kB, ne_electron, ne_alpha, mic_corr, display,
            )

        # -------------------------------------------------------------- #
        # 8.5a  Optional: record quenched energies                        #
        # -------------------------------------------------------------- #
        if record:
            with open(recfile3, "a") as f:
                _write_particles(f, particle_vec, e_quenching, idec)
                if evenement != 1:
                    _write_particles(f, particle_vec2, e_quenching2, idec, t1)

        # -------------------------------------------------------------- #
        # 8.6  TDCR detection probabilities                               #
        # -------------------------------------------------------------- #
        if evenement == 1:
            t1 = 0

        detect_fn = (
            tl.detectProbabilitiesMC if fullMC else tl.detectProbabilities
        )
        (
            eff0_S, eff0_D, eff0_T,
            eff0_AB, eff0_BC, eff0_AC, eff0_D2,
        ) = detect_fn(L, e_quenching, e_quenching2, t1, evenement, ext_dt, meas_time)

        eff_lists["S"].append(eff0_S)
        eff_lists["D"].append(eff0_D)
        eff_lists["T"].append(eff0_T)
        eff_lists["AB"].append(eff0_AB)
        eff_lists["BC"].append(eff0_BC)
        eff_lists["AC"].append(eff0_AC)
        eff_lists["D2"].append(eff0_D2)

        # -------------------------------------------------------------- #
        # 8.6a  Optional: record detection probabilities                  #
        # -------------------------------------------------------------- #
        if record:
            with open(recfile4, "a") as f:
                line = f"{idec} {eff0_S} {eff0_D} {eff0_T}"
                if not symm:
                    line += f" {eff0_AB} {eff0_BC} {eff0_AC}"
                f.write(line + "\n")

    # ------------------------------------------------------------------ #
    # 9. Final estimators                                                  #
    # ------------------------------------------------------------------ #
    out_eff = tl.efficienciesEstimates(
        eff_lists["S"], eff_lists["D"], eff_lists["T"],
        eff_lists["AB"], eff_lists["BC"], eff_lists["AC"],
        eff_lists["D2"], N,
    )

    if mode == "eff":
        return out_eff
    if mode == "dis":
        return (
            eff_lists["S"], eff_lists["D"],
            eff_lists["T"], eff_lists["D2"],
        )


# ---------------------------------------------------------------------------
# Private sub-routines extracted from the MC loop
# ---------------------------------------------------------------------------

def _load_decay_data(rad_list):
    """
    Read PenNuc decay data for each radionuclide in *rad_list*.

    Parameters
    ----------
    rad_list : list of str
        Radionuclide labels.

    Returns
    -------
    tuple
        17-tuple of lists extracted from
        :func:`tdcrpy.TDCR_model_lib.readPenNuc2`, one entry per nuclide.
    """
    particle, p_branch, e_branch, LevelDaughter = [], [], [], []
    levelNumber, prob_trans, prob_branch, levelEnergy = [], [], [], []
    transitionType, e_trans, next_level, Q_value = [], [], [], []
    DaughterVec, Pdaughter = [], []
    Transition_prob_sum, u_prob_trans, trans_halfLife = [], [], []

    for rad_i in rad_list:
        pn = tl.readPenNuc2(rad_i)
        DaughterVec.append(pn[0])
        Pdaughter.append(pn[1])
        Q_value.append(pn[2])
        particle.append(pn[3])
        e_branch.append(pn[4])
        p_branch.append(pn[5])
        LevelDaughter.append(pn[6])
        prob_branch.append(pn[7])
        transitionType.append(pn[8])
        e_trans.append(pn[9])
        prob_trans.append(pn[10])
        levelNumber.append(pn[11])
        next_level.append(pn[12])
        levelEnergy.append(pn[13])
        Transition_prob_sum.append(pn[14])
        trans_halfLife.append(pn[15])
        u_prob_trans.append(pn[16])

    return (
        particle, p_branch, e_branch, LevelDaughter,
        levelNumber, prob_trans, prob_branch, levelEnergy,
        transitionType, e_trans, next_level, Q_value,
        DaughterVec, Pdaughter, Transition_prob_sum,
        u_prob_trans, trans_halfLife,
    )


def _sample_nuclear_decay(
    index_rad, Rad, DaughterVec, Pdaughter, particle, p_branch,
    e_branch, LevelDaughter, levelNumber, prob_trans, prob_branch,
    levelEnergy, transitionType, e_trans, next_level,
    Transition_prob_sum, trans_halfLife, u_prob_trans,
    tau, unc_data, display,
):
    """
    Sample one complete nuclear decay event.

    Performs daughter selection, branch selection, and the full isomeric
    transition cascade for both prompt and (if applicable) delayed events.

    Returns
    -------
    particle_vec, energy_vec : list
        Prompt particles and energies.
    particle_vec2, energy_vec2 : list
        Delayed particles and energies (empty if ``evenement == 1``).
    evenement : int
        Number of separate coincidence windows (1 = prompt only).
    t1 : float
        Decay time (s) of the last isomeric transition sampled.
    level_before_trans : float
        Daughter energy level just after the primary particle emission.
    iDaughter : int
        Sampled daughter index.
    """
    particle_vec = []
    energy_vec = []
    particle_vec2 = []
    energy_vec2 = []
    evenement = 1
    t1 = 0

    # ---- Daughter sampling ----
    iDaughter = tl.sampling(
        np.asarray(Pdaughter[index_rad]) / sum(np.asarray(Pdaughter[index_rad]))
    )

    if display:
        print(f"\t Sampled daughter: {DaughterVec[index_rad][iDaughter]}")

    # ---- Branch sampling ----
    branch_i = tl.sampling(tl.normalise(prob_branch[index_rad][iDaughter]))
    level_before_trans = 0.0

    if p_branch[index_rad][iDaughter][branch_i] != []:
        branch_proba = tl.normalise(p_branch[index_rad][iDaughter][branch_i])
        isub = tl.sampling(branch_proba)
        particle_branch = particle[index_rad][iDaughter][branch_i][isub]
        energy_branch = e_branch[index_rad][iDaughter][branch_i][isub]
        levelOftheDaughter = LevelDaughter[index_rad][iDaughter][branch_i][isub]
        level_before_trans = levelOftheDaughter

        if display:
            _print_branch(particle_branch, energy_branch, levelOftheDaughter)

        particle_vec.append(particle_branch)
        energy_vec.append(energy_branch)
    else:
        # Isomeric transition with no primary particle emission.
        tp = tl.normalise(Transition_prob_sum[index_rad][iDaughter])
        it = tl.sampling(tp)
        levelOftheDaughter = levelNumber[index_rad][iDaughter][it][0]
        if display:
            print("\t Sampled branch: isomeric transition — no particle emitted")
            print(f"\t\t Level: {levelOftheDaughter}")

    # ---- Isomeric transition cascade ----
    while levelOftheDaughter > 0:
        i_level = levelNumber[index_rad][iDaughter].index([levelOftheDaughter])
        t1 = np.random.exponential(
            trans_halfLife[index_rad][iDaughter][i_level][0] / np.log(2)
        )

        if t1 > tau * 1e-9:
            evenement += 1
            if display:
                print(
                    f"\t\t Delayed transition: t = {t1 * 1e9:.2f} ns "
                    f"(τ½ = {trans_halfLife[index_rad][iDaughter][i_level][0] * 1e9:.2f} ns)"
                )

        if transitionType[index_rad][iDaughter][i_level] != []:
            if unc_data:
                prob_s = [
                    np.random.normal(xpt, u_prob_trans[index_rad][iDaughter][i_level][ipt])
                    for ipt, xpt in enumerate(prob_trans[index_rad][iDaughter][i_level])
                ]
                prob_t = tl.normalise(prob_s)
            else:
                prob_t = tl.normalise(prob_trans[index_rad][iDaughter][i_level])

            index_t = tl.sampling(prob_t)
            t_type = transitionType[index_rad][iDaughter][i_level][index_t]
            e_t = e_trans[index_rad][iDaughter][i_level][index_t]

            if display:
                _print_transition(
                    levelEnergy[index_rad][iDaughter][i_level][0], t_type, e_t,
                    next_level[index_rad][iDaughter][i_level][index_t],
                )

            # Append to the correct event window.
            target_pv = particle_vec2 if evenement != 1 else particle_vec
            target_ev = energy_vec2 if evenement != 1 else energy_vec
            _score_transition(t_type, e_t, target_pv, target_ev)

            levelOftheDaughter = next_level[index_rad][iDaughter][i_level][index_t]

        else:
            print(
                "Warning: no transition data for "
                f"{DaughterVec[index_rad][iDaughter]}, "
                f"level {levelOftheDaughter}"
            )
            levelOftheDaughter = 0

    return (
        particle_vec, energy_vec,
        particle_vec2, energy_vec2,
        evenement, t1,
        level_before_trans, iDaughter,
    )


def _score_transition(t_type, e_t, particle_vec, energy_vec):
    """
    Append the particle(s) produced by one isomeric transition to
    *particle_vec* / *energy_vec*.

    Parameters
    ----------
    t_type : str
        Transition type code (``"GA"`` = gamma, ``"EK"`` … ``"EN"`` =
        internal conversion from shell K … N).
    e_t : float
        Transition energy (keV).
    particle_vec, energy_vec : list
        Lists updated in-place.
    """
    if t_type == "GA":
        particle_vec.append("gamma")
        energy_vec.append(e_t)
    else:
        # Internal conversion electron.
        particle_vec.append("electron")
        energy_vec.append(e_t)
        # Record the resulting shell vacancy.
        shell_map = {
            "EK": "Atom_K", "EL": "Atom_L",
            "EL1": "Atom_L1", "EL2": "Atom_L2", "EL3": "Atom_L3",
            "EM": "Atom_M", "EN": "Atom_N",
        }
        vacancy = shell_map.get(t_type)
        if vacancy:
            particle_vec.append(vacancy)
            energy_vec.append(0)


def _sample_beta_spectra(particle_vec, energy_vec, rad_i, level_before_trans, display):
    """
    Replace ``"beta"`` / ``"beta+"`` tokens with sampled kinetic energies
    from the BetaShape spectra.

    For ``"beta+"``, two 511 keV annihilation photons are also appended.

    Parameters
    ----------
    particle_vec : list of str
        Particle labels (modified in-place).
    energy_vec : list of float
        Energies in keV (modified in-place).
    rad_i : str
        Radionuclide label (e.g. ``"Co-60"``).
    level_before_trans : float
        Energy level of the daughter just before the transition.
    display : bool
        If ``True``, print sampled energies.

    Returns
    -------
    particle_vec, energy_vec : list
        The same lists, now with beta tokens replaced.
    """
    for ipart, p in enumerate(particle_vec):
        if p == "beta":
            e_b, p_b = tl.readBetaShape(rad_i, "beta-", level_before_trans)
            idx = tl.sampling(p_b)
            particle_vec[ipart] = "electron"
            energy_vec[ipart] = e_b[idx]
            if display:
                print(f"\t\t beta- sampled energy = {energy_vec[ipart]:.3f} keV")

        elif p == "beta+":
            e_b, p_b = tl.readBetaShape(rad_i, "beta+", level_before_trans)
            idx = tl.sampling(p_b)
            particle_vec[ipart] = "positron"
            energy_vec[ipart] = e_b[idx]
            # Two 511 keV annihilation photons.
            particle_vec.extend(["gamma", "gamma"])
            energy_vec.extend([511.0, 511.0])
            if display:
                print(f"\t\t beta+ sampled energy = {energy_vec[ipart]:.3f} keV")

    return particle_vec, energy_vec


def _interact_prompt(particle_vec, energy_vec, V):
    """
    Apply radiation–matter interactions (photoelectric, pair production,
    Compton scattering, beta/alpha energy deposition) to the particle list.

    Parameters
    ----------
    particle_vec : list of str
        Particle labels (modified in-place).
    energy_vec : list of float
        Energies in keV (modified in-place).
    V : float
        Scintillator volume (mL).

    Returns
    -------
    particle_vec, energy_vec : list
        Updated lists.
    """
    for ipart, p in enumerate(particle_vec):
        if p in ("electron", "beta+"):
            energy_vec[ipart] = tl.energie_dep_beta2(energy_vec[ipart], v=V)
            particle_vec[ipart] = "electron"

        elif p in ("gamma", "XKA", "XKB", "XL"):
            Ei = energy_vec[ipart]
            Ed = tl.energie_dep_gamma2(Ei, v=V)

            if Ei == Ed:
                # Photoelectric effect.
                e_pe, lacune, element = tl.interaction_scintillation(Ed)
                _, e_sec, _, par_sec = tl.relaxation_atom_ph(lacune, element, v=V)
                particle_vec[ipart] = "electron"
                energy_vec[ipart] = e_pe
                # Append secondary particles from atomic relaxation.
                particle_vec.extend(par_sec)
                energy_vec.extend(e_sec)
            elif Ed == Ei - 1022:
                # Pair production.
                E_e = (Ei - 1022) / 2.0
                energy_vec[ipart] = E_e
                particle_vec.append("positron")
                energy_vec.append(E_e)
            else:
                # Compton scattering.
                energy_vec[ipart] = Ed

            particle_vec[ipart] = "electron"

        elif "Auger" in p:
            energy_vec[ipart] = tl.energie_dep_beta2(energy_vec[ipart], v=V)
            particle_vec[ipart] = "electron"

    return particle_vec, energy_vec


def _quench(particle_vec, energy_vec, energy_vec_initial,
            kB, ne_electron, ne_alpha, mic_corr, display):
    """
    Apply the Birks quenching model to each particle.

    Parameters
    ----------
    particle_vec : list of str
    energy_vec : list of float
        Deposited energies (keV).
    energy_vec_initial : list of float
        Initial (pre-interaction) energies (keV), needed for the electron
        quenching integral.
    kB : float
        Birks constant (cm keV⁻¹).
    ne_electron : int
        Number of integration bins for electron quenching.
    ne_alpha : int
        Number of integration bins for alpha quenching.
    mic_corr : bool
        Apply reverse-micelle correction when ``True``.
    display : bool

    Returns
    -------
    e_quenching : list of float
        Quenched energies in keV (zero for non-scintillating particles).
    """
    e_quenching = []
    for ipart, p in enumerate(particle_vec):
        if p == "alpha":
            energy_vec[ipart] = tl.Em_a(energy_vec[ipart], kB, ne_alpha)
            e_quenching.append(energy_vec[ipart])

        elif p in ("electron", "positron"):
            eq = tl.Em_e(
                energy_vec_initial[ipart] * 1e3,
                energy_vec[ipart] * 1e3,
                kB * 1e3,
                ne_electron,
            ) * 1e-3
            if mic_corr:
                eq = tl.pure_mc_efficient_energy_numba(eq)
            energy_vec[ipart] = eq
            e_quenching.append(eq)

        else:
            e_quenching.append(0)

    if display:
        for ipart, p in enumerate(particle_vec):
            if not p.startswith("Atom"):
                print(
                    f"\t\t quenched energy of {p} = "
                    f"{np.round(e_quenching[ipart], 3)} keV"
                )

    return e_quenching


def _run_from_history(L, N, tau, ext_dt, meas_time, fullMC, mode):
    """
    Replay detection-probability computation from a previously recorded
    quenched-energy history file (``Temp_E2.txt``).

    Parameters
    ----------
    L : float or tuple
        Free parameter(s) forwarded to the detection model.
    N : int
        Number of recorded decays (used for uncertainty normalisation).
    tau : float
        Coincidence resolving time (ns).
    ext_dt : float
        Extended dead time (µs).
    meas_time : float
        Measurement time (min).
    fullMC : bool
        Use the fully stochastic detection model when ``True``.
    mode : str
        ``"eff"`` or ``"dis"``.

    Returns
    -------
    Same return type as :func:`TDCRPy` for the selected *mode*.
    """
    temp_dir = tempfile.gettempdir()
    recfile3 = os.path.join(temp_dir, "Temp_E2.txt")

    eff_lists = {k: [] for k in ("S", "D", "T", "AB", "BC", "AC", "D2")}
    detect_fn = (
        tl.detectProbabilitiesMC if fullMC else tl.detectProbabilities
    )

    with open(recfile3, "r") as f:
        decaym = -1
        e_quenching = []
        e_quenching2 = []
        evenement = 1
        t1 = 0.0

        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            decay = int(parts[2])

            if decay != decaym:
                if decay > 0:
                    # Flush previous event.
                    _append_efficiencies(
                        detect_fn, L, e_quenching, e_quenching2,
                        t1, evenement, ext_dt, meas_time, eff_lists,
                    )

                # Start new event.
                energy = float(parts[1]) * 1e-3
                t1 = float(parts[4])
                decaym = decay
                e_quenching = [energy]
                e_quenching2 = []
                evenement = 1

            else:
                energy = float(parts[1]) * 1e-3
                t1 = float(parts[4])
                if t1 > tau * 1e-9:
                    evenement += 1
                    e_quenching2.append(energy)
                else:
                    e_quenching.append(energy)

        # Flush last event.
        _append_efficiencies(
            detect_fn, L, e_quenching, e_quenching2,
            t1, evenement, ext_dt, meas_time, eff_lists,
        )

    out_eff = tl.efficienciesEstimates(
        eff_lists["S"], eff_lists["D"], eff_lists["T"],
        eff_lists["AB"], eff_lists["BC"], eff_lists["AC"],
        eff_lists["D2"], N,
    )
    if mode == "eff":
        return out_eff
    return (eff_lists["S"], eff_lists["D"], eff_lists["T"], eff_lists["D2"])


def _append_efficiencies(detect_fn, L, e_q, e_q2,
                          t1, evenement, ext_dt, meas_time, eff_lists):
    """Compute and store one event's detection probabilities."""
    eff0 = detect_fn(L, e_q, e_q2, t1, evenement, ext_dt, meas_time)
    keys = ("S", "D", "T", "AB", "BC", "AC", "D2")
    for k, v in zip(keys, eff0):
        eff_lists[k].append(v)


# ---------------------------------------------------------------------------
# Display helpers (only active when display=True)
# ---------------------------------------------------------------------------

def _print_relaxation(label, particle_vec, energy_vec):
    print(f"\n\t ATOMIC RECOMBINATION — {label}")
    for i, p in enumerate(particle_vec):
        if p.startswith("Atom"):
            print(f"\t\t electron left the {p[5:]} shell")
        elif p in ("beta", "beta+"):
            print(f"\t\t {p} energy = {energy_vec[i]:.3f} keV")
        else:
            print(f"\t\t emitted {p} energy = {energy_vec[i]:.3f} keV")


def _print_interaction(label, particle_vec, energy_vec):
    print(f"\n\t INTERACTION — {label}")
    for i, p in enumerate(particle_vec):
        if not p.startswith("Atom") and energy_vec[i] != 0:
            print(f"\t\t {p} deposited energy = {energy_vec[i]:.3f} keV")


def _print_branch(particle_branch, energy_branch, level):
    print("\t Sampled decay branch:")
    if particle_branch.startswith("Atom"):
        shell = particle_branch[5:]
        print(f"\t\t Electron capture on {shell} shell")
    else:
        print(f"\t\t Particle: {particle_branch}")
        print(f"\t\t Energy = {energy_branch} keV")
    print(f"\t\t Daughter level: {level}")


def _print_transition(level_e, t_type, e_t, next_lev):
    label_map = {
        "GA": "gamma ray",
        "EK": "conversion electron (K shell)",
        "EL": "conversion electron (L shell)",
        "EL1": "conversion electron (L1 shell)",
        "EL2": "conversion electron (L2 shell)",
        "EL3": "conversion electron (L3 shell)",
        "EM": "conversion electron (M shell)",
        "EN": "conversion electron (N shell)",
        "EO": "conversion electron (O shell)",
    }
    print(f"\t\t Level energy = {level_e} keV")
    print(f"\t\t {label_map.get(t_type, t_type)} energy = {e_t} keV")
    print(f"\t\t Next level = {next_lev}")


# ---------------------------------------------------------------------------
# Objective function and efficiency-from-measurement wrappers
# ---------------------------------------------------------------------------

def objectFct(L, TD, Rad, pmf_1, N, kB, V):
    """
    Objective function minimised by :func:`eff` to fit *L* to a measured
    TDCR ratio.

    The function reads pre-recorded decay histories (written by a prior
    call to :func:`TDCRPy` with ``record=True``) so that the expensive
    nuclear-decay simulation is executed only once.

    Parameters
    ----------
    L : float or tuple of float
        Free parameter(s) (keV⁻¹).
    TD : float or tuple of float
        Measured TDCR ratio(s). A scalar implies the global T/D ratio; a
        4-tuple ``(TD, T/AB, T/BC, T/AC)`` enables asymmetric fitting.
    Rad : str
        Radionuclide list (comma-separated).
    pmf_1 : str
        Activity fractions (comma-separated).
    N : int
        Number of MC trials (passed through to :func:`TDCRPy`).
    kB : float
        Birks constant (cm keV⁻¹).
    V : float
        Scintillator volume (mL).

    Returns
    -------
    res : float
        Sum of squared residuals between the modelled and measured
        coincidence ratios.
    """
    symm = not isinstance(TD, (tuple, list))
    eff_model = TDCRPy(L, Rad, pmf_1, N, kB, V, readRecHist=True)

    if symm:
        tdcr_calc = eff_model[4] / eff_model[2]  # T / D
        return (tdcr_calc - TD) ** 2

    # Asymmetric: minimise sum-of-squares over the three partial ratios.
    tab = eff_model[4] / eff_model[6]   # T / AB
    tbc = eff_model[4] / eff_model[8]   # T / BC
    tac = eff_model[4] / eff_model[10]  # T / AC
    return (TD[1] - tab) ** 2 + (TD[2] - tbc) ** 2 + (TD[3] - tac) ** 2


def eff(TD, Rad, pmf_1, kB, V,
        N=10000, L=1, maxiter=20, xatol=1e-7,
        disp=False, Lbounds=(0.1, 10)):
    """
    Determine the free parameter *L* and detection efficiencies from a
    measured TDCR ratio, using the full stochastic Monte Carlo model.

    The procedure is:

    1. Run one full MC simulation (``record=True``) to cache the decay
       histories.
    2. Optimise *L* by minimising :func:`objectFct` via a bounded
       scalar search (symmetric case) or Nelder–Mead simplex
       (asymmetric case).
    3. Re-evaluate efficiencies at the optimal *L* using the cached
       histories (``readRecHist=True``).

    Parameters
    ----------
    TD : float or tuple of float
        Measured TDCR ratio(s).  A scalar ``T/D`` implies a symmetric
        PMT configuration.  A 4-tuple ``(T/D, T/AB, T/BC, T/AC)`` is
        required for the asymmetric case.
    Rad : str
        Comma-separated radionuclide labels.
    pmf_1 : str
        Comma-separated activity fractions.
    kB : float
        Birks constant (cm keV⁻¹).
    V : float
        Scintillator volume (mL).
    N : int, optional
        Number of MC trials for the initial simulation.  Default is 10 000.
    L : float, optional
        Initial guess for the free parameter (keV⁻¹).  Default is 1.
    maxiter : int, optional
        Maximum number of optimiser iterations.  Default is 20.
    xatol : float, optional
        Absolute tolerance on *L* for the Nelder–Mead convergence criterion.
        Default is ``1e-7``.
    disp : bool, optional
        Pass ``True`` to print optimiser diagnostics.  Default is ``False``.
    Lbounds : tuple of float, optional
        ``(lower, upper)`` bounds on *L* for the bounded scalar search.
        Default is ``(0.1, 10)``.

    Returns
    -------
    L0 : float
        Optimised global free parameter (keV⁻¹).
    L : tuple of float
        Triplet of free parameters ``(L_A, L_B, L_C)``.
        Equal to ``(L0, L0, L0)`` for the symmetric case.
    eff_S, u_eff_S : float
        Single-event counting efficiency and its standard uncertainty.
    eff_D, u_eff_D : float
        Double-coincidence counting efficiency and its standard uncertainty.
    eff_T, u_eff_T : float
        Triple-coincidence counting efficiency and its standard uncertainty.
    eff_AB, u_eff_AB : float
        AB coincidence efficiency and standard uncertainty.
    eff_BC, u_eff_BC : float
        BC coincidence efficiency and standard uncertainty.
    eff_AC, u_eff_AC : float
        AC coincidence efficiency and standard uncertainty.
    eff_D2, u_eff_D2 : float
        C/N coincidence efficiency and standard uncertainty.

    See Also
    --------
    effA : Faster analytical alternative.
    TDCRPy : Underlying MC engine.
    """
    symm = not isinstance(TD, (tuple, list))

    # Step 1 — cache decay histories.
    TDCRPy(L, Rad, pmf_1, N, kB, V, record=True)

    # Step 2 — scalar optimisation for L (symmetric).
    td0 = TD if symm else TD[0]
    result = opt.minimize_scalar(
        objectFct,
        args=(td0, Rad, pmf_1, N, kB, V),
        method="bounded",
        bounds=(Lbounds[0], Lbounds[1]),
        options={"disp": disp, "maxiter": maxiter},
    )
    L0 = result.x
    L_opt = (L0, L0, L0)

    # Step 2b — Nelder–Mead refinement for asymmetric case.
    if not symm:
        result2 = opt.minimize(
            objectFct,
            L_opt,
            args=(TD, Rad, pmf_1, N, kB, V),
            method="nelder-mead",
            options={"xatol": xatol, "disp": disp, "maxiter": maxiter},
        )
        L_opt = tuple(result2.x)

    # Step 3 — final efficiency evaluation.
    L_final = L0 if symm else L_opt
    out = TDCRPy(L_final, Rad, pmf_1, N, kB, V, readRecHist=True)

    (
        eff_S, u_eff_S, eff_D, u_eff_D,
        eff_T, u_eff_T,
        eff_AB, u_eff_AB, eff_BC, u_eff_BC, eff_AC, u_eff_AC,
        eff_D2, u_eff_D2,
    ) = out

    return (
        L0, L_opt,
        eff_S, u_eff_S,
        eff_D, u_eff_D,
        eff_T, u_eff_T,
        eff_AB, u_eff_AB,
        eff_BC, u_eff_BC,
        eff_AC, u_eff_AC,
        eff_D2, u_eff_D2,
    )


def effA(TD, Rad, pmf_1, kB, V,
         L=1, maxiter=20, xatol=1e-7,
         disp=False, Lbounds=(0.1, 10),
         cerenkov=False):
    """
    Determine the free parameter *L* and detection efficiencies from a
    measured TDCR ratio, using the **analytical** model.

    The analytical model is significantly faster than the stochastic MC
    model but is only valid for pure beta emitters.  For other nuclides,
    the result is an approximation.

    Parameters
    ----------
    TD : float or tuple of float
        Measured TDCR ratio(s).  See :func:`eff` for the expected format.
    Rad : str
        Radionuclide list (comma-separated).
    pmf_1 : str
        Activity fractions (comma-separated).
    kB : float
        Birks constant (cm keV⁻¹).
    V : float
        Scintillator volume (mL).
    L : float, optional
        Initial guess for the free parameter (keV⁻¹).  Default is 1.
    maxiter : int, optional
        Maximum optimiser iterations.  Default is 20.
    xatol : float, optional
        Nelder–Mead convergence tolerance.  Default is ``1e-7``.
    disp : bool, optional
        Print optimiser diagnostics.  Default is ``False``.
    Lbounds : tuple of float, optional
        ``(lower, upper)`` bounds for the scalar search.  Default is
        ``(0.1, 10)``.
    cerenkov : bool, optional
        Use the Čerenkov model (:func:`~tdcrpy.TDCR_model_lib.modelCerenkov`)
        instead of the standard analytical model.  Default is ``False``.

    Returns
    -------
    L0 : float
        Optimised global free parameter (keV⁻¹).
    L : tuple of float
        Triplet ``(L_A, L_B, L_C)``.
    eff_S : float
        Single-event counting efficiency.
    eff_D : float
        Double-coincidence counting efficiency.
    eff_T : float
        Triple-coincidence counting efficiency.

    Notes
    -----
    Unlike :func:`eff`, this function does not return per-channel
    uncertainties because the analytical model does not produce MC
    sampling noise.

    See Also
    --------
    eff : Full stochastic variant with uncertainty output.
    """
    symm = not isinstance(TD, (tuple, list))
    model_fn = tl.modelCerenkov if cerenkov else tl.modelAnalytical

    def _args_scalar(td_scalar):
        if cerenkov:
            return (td_scalar, td_scalar, td_scalar, td_scalar, Rad, "res")
        return (td_scalar, td_scalar, td_scalar, td_scalar, Rad, kB, V, "res", 1e3)

    def _args_vec(td_tuple):
        if cerenkov:
            return (td_tuple[0], td_tuple[1], td_tuple[2], td_tuple[3], Rad, "res")
        return (td_tuple[0], td_tuple[1], td_tuple[2], td_tuple[3],
                Rad, kB, V, "res", 1e3)

    # Scalar bounded search.
    td0 = TD if symm else TD[0]
    result = opt.minimize_scalar(
        model_fn,
        args=_args_scalar(td0),
        method="bounded",
        bounds=(Lbounds[0], Lbounds[1]),
        options={"disp": disp, "maxiter": maxiter},
    )
    L0 = result.x
    L_opt = (L0, L0, L0)

    # Nelder–Mead refinement for asymmetric case.
    if not symm:
        result2 = opt.minimize(
            model_fn,
            L_opt,
            args=_args_vec(TD),
            method="nelder-mead",
            options={"xatol": xatol, "disp": disp, "maxiter": maxiter},
        )
        L_opt = tuple(result2.x)

    # Final efficiency evaluation.
    L_final = L0 if symm else L_opt
    if cerenkov:
        if symm:
            out = model_fn(L_final, TD, TD, TD, TD, Rad, "eff")
        else:
            out = model_fn(L_final, TD[0], TD[1], TD[2], TD[3], Rad, "eff")
    else:
        if symm:
            out = model_fn(L_final, TD, TD, TD, TD, Rad, kB, V, "eff", 1e3)
        else:
            out = model_fn(L_final, TD[0], TD[1], TD[2], TD[3],
                           Rad, kB, V, "eff", 1e3)

    eff_S, eff_D, eff_T = out[0], out[1], out[2]

    return L0, L_opt, eff_S, eff_D, eff_T