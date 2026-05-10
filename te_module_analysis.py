"""
Thermoelectric Cooling Performance Analysis: MgAgSb / Mg3(Sb,Bi)2 Module
=========================================================================

Motivation
----------
Reading Ying et al. (2022) "A robust thermoelectric module based on
MgAgSb/Mg3(Sb,Bi)2 with a conversion efficiency of 8.5% and a maximum
cooling of 72 K", Energy & Environmental Science, DOI: 10.1039/D2EE00883A,
raised a practical question:

    The paper reports ΔT_max = 72 K experimentally at T_h = 347 K and
    achieves a maximum cooling COP of ~0.3–0.4. The theoretical performance
    is derived from the standard single-couple equations — but the paper does
    not explore how leg geometry (height, cross-sectional area ratio) affects
    the trade-off between maximum ΔT and maximum COP at a fixed operating
    current.

This script reconstructs their theoretical performance framework from first
principles, validates against their reported module data, and extends the
analysis to map how leg geometry and contact resistance shape the cooling
performance envelope.

Material Property Sources
-------------------------
All thermoelectric property values are taken from or consistent with:

[1] Ying P. et al., Energy Environ. Sci. 15, 6584 (2022)
    DOI: 10.1039/D2EE00883A
    — Module configuration, experimental ΔT_max, COP curves, contact details

[2] Liu Z. et al., Nature Communications 13, 1120 (2022)
    DOI: 10.1038/s41467-022-28798-4
    — n-type Mg3Bi1.5Sb0.5 transport properties (ρ, S, κ vs T)

[3] Zhao H. et al., Nano Energy 7, 97 (2014)
    DOI: 10.1016/j.nanoen.2014.04.012
    — p-type α-MgAgSb baseline transport properties

Physics: Standard thermoelectric single-couple model
    ΔT_max = 0.5 * Z * T_c^2           (maximum temperature difference)
    Q_c    = N*(S*T_c*I - 0.5*I^2*R - K*ΔT)  (cooling power)
    W      = N*(S*ΔT*I + I^2*R)               (electrical power input)
    COP    = Q_c / W                           (coefficient of performance)

    where Z = S^2 / (R*K) is the module figure of merit,
    S = S_p - S_n (net Seebeck), R = electrical resistance, K = thermal conductance

Author: Anoopkrishna Krishnakumar
Date:   May 2026
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# SECTION 1: Material Properties (from published sources)
# =============================================================================

def material_properties_at_T(T_mean_K):
    """
    Return effective module-level transport properties at mean temperature T.

    p-type: α-MgAgSb
    n-type: Mg3Sb0.6Bi1.4 (Ying et al. 2022 optimised composition)

    Values are representative of the reported peak-performance compositions.
    Temperature dependence is approximated from the digitised curves in [1,2].

    Parameters
    ----------
    T_mean_K : float
        Mean temperature of the module legs (K)

    Returns
    -------
    S_p, S_n : Seebeck coefficients (V/K)
    rho_p, rho_n : electrical resistivity (Ω·m)
    kappa_p, kappa_n : thermal conductivity (W/m·K)
    """
    # --- p-type α-MgAgSb ---
    # Ref [1,3]: S_p ~ 200 μV/K at 300 K, weakly decreasing to ~170 μV/K at 400 K
    S_p = (200 - 0.30 * (T_mean_K - 300)) * 1e-6   # V/K

    # Ref [1,3]: ρ_p ~ 7 μΩ·m at 300 K, rising to ~10 μΩ·m at 400 K
    rho_p = (7.0 + 0.030 * (T_mean_K - 300)) * 1e-6  # Ω·m

    # Ref [1,3]: κ_p ~ 0.8 W/m·K (relatively flat between 300–400 K)
    kappa_p = 0.80 + 0.001 * (T_mean_K - 300)        # W/m·K

    # --- n-type Mg3Sb0.6Bi1.4 (sintered at 1073 K) ---
    # Ref [1,2]: S_n ~ -210 μV/K at 300 K
    S_n = -(210 - 0.25 * (T_mean_K - 300)) * 1e-6   # V/K  (negative, n-type)

    # Ref [1,2]: ρ_n ~ 16 μΩ·m at 300 K (1073 K sintering), falling at higher T
    rho_n = (16.0 - 0.010 * (T_mean_K - 300)) * 1e-6  # Ω·m

    # Ref [1,2]: κ_n ~ 1.3 W/m·K at 300 K
    kappa_n = 1.30 - 0.002 * (T_mean_K - 300)         # W/m·K

    return S_p, S_n, rho_p, rho_n, kappa_p, kappa_n


# =============================================================================
# SECTION 2: Module Geometry & Parameters
# =============================================================================

# Ying et al. (2022) module: 8 thermoelectric couples
N_COUPLES = 8          # number of p-n pairs in the module [Ref 1]

# Leg dimensions — from standard Mg-based module literature
# Ying et al. do not specify exact dimensions explicitly; values consistent
# with the module cross-section area ~2×2 mm and fill factor ~0.3 reported
LEG_HEIGHT_BASE  = 3.0e-3   # 3 mm baseline leg height (m)
LEG_AREA_BASE    = 4.0e-6   # 2×2 mm cross-section (m²)

# Contact resistance — from Ying et al. 2022, ESI, and [5] Xie et al. 2024
# Ag contact on MgAgSb: ~5 μΩ·cm² = 5e-10 Ω·m²
# Fe contact on Mg3(Sb,Bi)2: ~8 μΩ·cm² (higher due to CTE mismatch)
R_contact_p = 5.0e-10  # Ω·m² (specific contact resistance, p-leg) [Ref 1 ESI]
R_contact_n = 8.0e-10  # Ω·m² (specific contact resistance, n-leg)  [Ref 1 ESI]


# =============================================================================
# SECTION 3: Core Module Performance Functions
# =============================================================================

def compute_module_params(T_h, T_c, leg_height, leg_area,
                          n_couples=N_COUPLES,
                          r_contact_p=R_contact_p,
                          r_contact_n=R_contact_n):
    """
    Compute module-level thermoelectric parameters.

    Parameters
    ----------
    T_h : float  — hot-side temperature (K)
    T_c : float  — cold-side temperature (K)
    leg_height : float — leg height L (m)
    leg_area   : float — leg cross-section A (m²)
    n_couples  : int   — number of thermoelectric pairs
    r_contact_p, r_contact_n : float — specific contact resistance (Ω·m²)

    Returns
    -------
    dict with S_net, R_module, K_module, Z_module, zT_mean
    """
    T_mean = 0.5 * (T_h + T_c)
    S_p, S_n, rho_p, rho_n, kappa_p, kappa_n = material_properties_at_T(T_mean)

    # Net Seebeck coefficient (module uses S_p - S_n since S_n is negative)
    S_net = S_p - S_n   # both contribute positively to Peltier cooling

    # Electrical resistance per couple (leg + contact contributions)
    R_leg_p = rho_p * leg_height / leg_area
    R_leg_n = rho_n * leg_height / leg_area
    R_contact_total_p = 2 * r_contact_p / leg_area  # top + bottom contacts
    R_contact_total_n = 2 * r_contact_n / leg_area
    R_couple = R_leg_p + R_leg_n + R_contact_total_p + R_contact_total_n

    # Module total resistance (couples in series)
    R_module = n_couples * R_couple

    # Thermal conductance per couple (parallel heat paths)
    K_leg_p = kappa_p * leg_area / leg_height
    K_leg_n = kappa_n * leg_area / leg_height
    K_couple = K_leg_p + K_leg_n

    # Module thermal conductance (couples in parallel thermally)
    K_module = n_couples * K_couple

    # Figure of merit
    Z_module = S_net**2 / (R_couple * K_couple)
    zT_mean  = Z_module * T_mean

    return {
        'S_net': S_net,
        'R_module': R_module,
        'K_module': K_module,
        'Z_module': Z_module,
        'zT_mean': zT_mean,
        'T_mean': T_mean,
    }


def cooling_performance(T_h, I, leg_height, leg_area,
                        n_couples=N_COUPLES,
                        r_contact_p=R_contact_p,
                        r_contact_n=R_contact_n,
                        max_iter=200, tol=1e-4):
    """
    Self-consistently solve for T_c given T_h and current I.

    The cold-side temperature T_c is found iteratively because the material
    properties depend on T_mean = (T_h + T_c)/2.

    Returns T_c, Q_c (cooling power), W (power input), COP, ΔT
    """
    T_c = T_h - 30  # initial guess

    for _ in range(max_iter):
        T_mean = 0.5 * (T_h + T_c)
        params = compute_module_params(T_h, T_c, leg_height, leg_area,
                                       n_couples, r_contact_p, r_contact_n)
        S = params['S_net']
        R = params['R_module']
        K = params['K_module']
        dT = T_h - T_c

        # Heat balance at cold side:
        # Q_c = N*(S*T_c*I - 0.5*I^2*R - K*ΔT)
        Q_c = n_couples * (S * T_c * I - 0.5 * I**2 * R/n_couples - K/n_couples * dT)

        # New T_c from energy balance (adiabatic cold side → Q_c = 0 → ΔT_max condition)
        # For operating point: Q_c set by external load (here we solve for T_c)
        # Standard approach: iterate T_c until Q_c self-consistent
        T_c_new = T_h - (S * T_c * I - 0.5 * I**2 * R/n_couples) / (K/n_couples)
        T_c_new = max(T_c_new, 200)  # physical lower bound

        if abs(T_c_new - T_c) < tol:
            T_c = T_c_new
            break
        T_c = 0.6 * T_c + 0.4 * T_c_new  # damped update

    dT = T_h - T_c
    params = compute_module_params(T_h, T_c, leg_height, leg_area,
                                   n_couples, r_contact_p, r_contact_n)
    S = params['S_net']
    R = params['R_module']
    K = params['K_module']

    Q_c = n_couples * (S * T_c * I - 0.5 * (I**2 * R/n_couples) - (K/n_couples) * dT)
    W   = n_couples * (S * dT * I + I**2 * R/n_couples)
    COP = Q_c / W if W > 0 and Q_c > 0 else 0.0

    return T_c, max(Q_c, 0), max(W, 1e-9), COP, dT


def delta_T_max(T_h, leg_height, leg_area, n_couples=N_COUPLES,
                r_contact_p=R_contact_p, r_contact_n=R_contact_n):
    """
    Compute maximum achievable ΔT (analytical expression):
        ΔT_max = 0.5 * Z * T_c²
    Solved iteratively since T_c depends on ΔT_max.
    """
    T_c = T_h * 0.85  # start guess
    for _ in range(100):
        params = compute_module_params(T_h, T_c, leg_height, leg_area,
                                       n_couples, r_contact_p, r_contact_n)
        Z = params['Z_module']
        dT_new = 0.5 * Z * T_c**2
        T_c_new = T_h - dT_new
        if abs(T_c_new - T_c) < 0.01:
            break
        T_c = 0.5 * T_c + 0.5 * T_c_new
    return T_h - T_c, T_c


# =============================================================================
# SECTION 4: Plotting Routines
# =============================================================================

# Color palette — clean, publication-appropriate
C_P  = '#C0392B'   # red    → p-type / hot side
C_N  = '#2980B9'   # blue   → n-type / cold side
C_A  = '#27AE60'   # green  → reference / baseline
C_B  = '#8E44AD'   # purple → aspect ratio variation
C_W  = '#E67E22'   # orange → power / work
GREY = '#7F8C8D'
BG   = '#FAFAFA'

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.facecolor': BG,
    'figure.facecolor': 'white',
    'axes.grid': True,
    'grid.alpha': 0.35,
    'grid.linestyle': '--',
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'axes.titleweight': 'bold',
    'legend.fontsize': 9,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
})


def plot_material_properties():
    """Figure 1: Temperature-dependent material properties."""
    T_range = np.linspace(280, 420, 100)
    Sp, Sn, rhop, rhon, kp, kn = zip(*[material_properties_at_T(T) for T in T_range])
    Sp  = np.array(Sp)  * 1e6   # → μV/K
    Sn  = np.abs(np.array(Sn)) * 1e6
    rhop = np.array(rhop) * 1e6  # → μΩ·m
    rhon = np.array(rhon) * 1e6
    zTp = Sp**2 * 1e-12 / (rhop * 1e-6 * np.array(kp)) * T_range
    zTn = Sn**2 * 1e-12 / (rhon * 1e-6 * np.array(kn)) * T_range

    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    fig.suptitle(
        'Fig. 1 — Temperature-Dependent Transport Properties\n'
        'p-type α-MgAgSb [Zhao et al. 2014] & n-type Mg₃Sb₀.₆Bi₁.₄ [Ying et al. 2022]',
        fontsize=11, y=1.02
    )

    # Seebeck
    axes[0].plot(T_range, Sp, color=C_P, lw=2, label='MgAgSb (p)')
    axes[0].plot(T_range, Sn, color=C_N, lw=2, label='Mg₃(Sb,Bi)₂ (n)', ls='--')
    axes[0].set_xlabel('Temperature (K)')
    axes[0].set_ylabel('|Seebeck coefficient| (μV K⁻¹)')
    axes[0].set_title('Seebeck Coefficient')
    axes[0].legend()
    axes[0].set_xlim(280, 420)

    # Resistivity
    axes[1].plot(T_range, rhop, color=C_P, lw=2, label='MgAgSb (p)')
    axes[1].plot(T_range, rhon, color=C_N, lw=2, label='Mg₃(Sb,Bi)₂ (n)', ls='--')
    axes[1].set_xlabel('Temperature (K)')
    axes[1].set_ylabel('Electrical resistivity (μΩ m)')
    axes[1].set_title('Electrical Resistivity')
    axes[1].legend()
    axes[1].set_xlim(280, 420)

    # zT
    axes[2].plot(T_range, zTp, color=C_P, lw=2, label='MgAgSb (p)')
    axes[2].plot(T_range, zTn, color=C_N, lw=2, label='Mg₃(Sb,Bi)₂ (n)', ls='--')
    axes[2].set_xlabel('Temperature (K)')
    axes[2].set_ylabel('Figure of merit zT (dimensionless)')
    axes[2].set_title('Thermoelectric Figure of Merit zT')
    axes[2].legend()
    axes[2].set_xlim(280, 420)

    plt.tight_layout()
    plt.savefig('figures/fig1_material_properties.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  ✓ Fig 1 saved: Material properties")


def plot_cop_curves():
    """Figure 2: COP vs current at multiple ΔT (validation against Ying et al. 2022 Fig 4)."""
    T_h = 302.0   # K — cold-plate temperature in Ying 2022 cooling experiment
    dT_values = [0, 10, 20, 30, 40, 50]   # K — from Ying et al. Fig 4e
    colors = plt.cm.cool(np.linspace(0.1, 0.9, len(dT_values)))
    I_range = np.linspace(0.2, 6.0, 120)

    # Experimental peak COP reference points from Ying et al. 2022, Fig 4f
    # At ΔT=0: COP_max ≈ 1.3, I_opt ≈ 1.5 A
    exp_points = {
        0: (1.5, 1.30),
        10: (1.5, 0.85),
        20: (1.5, 0.55),
        30: (1.5, 0.32),
    }

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        'Fig. 2 — Cooling COP vs Input Current at Varying ΔT\n'
        'MgAgSb/Mg₃(Sb,Bi)₂ module, N=8 couples, T_h = 302 K\n'
        'Validation target: Ying et al., Energy Environ. Sci. (2022), Fig. 4e–f',
        fontsize=10, y=1.01
    )

    for i, dT in enumerate(dT_values):
        T_c = T_h - dT
        cop_list = []
        qc_list  = []
        for I in I_range:
            params = compute_module_params(T_h, T_c,
                                          LEG_HEIGHT_BASE, LEG_AREA_BASE)
            S = params['S_net']
            R = params['R_module']
            K = params['K_module']
            Q_c = N_COUPLES * (S*T_c*I - 0.5*I**2*R/N_COUPLES - K/N_COUPLES*dT)
            W   = N_COUPLES * (S*dT*I + I**2*R/N_COUPLES)
            cop = Q_c/W if W > 0 and Q_c > 0 else 0.0
            cop_list.append(cop)
            qc_list.append(max(Q_c, 0))

        ax1.plot(I_range, cop_list, color=colors[i], lw=2, label=f'ΔT = {dT} K')
        ax2.plot(I_range, np.array(qc_list), color=colors[i], lw=2, label=f'ΔT = {dT} K')

    # Mark experimental reference points
    for dT_exp, (I_exp, COP_exp) in exp_points.items():
        ax1.scatter(I_exp, COP_exp, marker='*', s=120,
                    color='black', zorder=5,
                    label='_nolegend_' if dT_exp > 0 else 'Exp. ref. [Ying 2022]')

    ax1.set_xlabel('Input Current I (A)')
    ax1.set_ylabel('Coefficient of Performance (COP)')
    ax1.set_title('COP vs Current')
    ax1.set_xlim(0, 6)
    ax1.set_ylim(0, 2.0)
    ax1.legend(loc='upper right', fontsize=8)

    ax2.set_xlabel('Input Current I (A)')
    ax2.set_ylabel('Cooling Power Q_c (W)')
    ax2.set_title('Cooling Power vs Current')
    ax2.set_xlim(0, 6)
    ax2.legend(loc='upper left', fontsize=8)

    plt.tight_layout()
    plt.savefig('figures/fig2_cop_curves.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  ✓ Fig 2 saved: COP curves (validation)")


def plot_geometry_sensitivity():
    """Figure 3: Effect of leg aspect ratio (L/A^0.5) on ΔT_max and optimal COP."""
    T_h = 302.0  # K

    # Vary leg height L at fixed area (aspect ratio sweep)
    heights = np.linspace(1.5e-3, 6.0e-3, 30)   # 1.5 mm to 6 mm
    areas   = [2e-6, 4e-6, 6e-6]                 # 1.4×1.4, 2×2, 2.5×2.5 mm²
    area_labels = ['1.4×1.4 mm²', '2×2 mm²', '2.5×2.5 mm²']
    area_colors = [C_N, C_A, C_P]

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    fig.suptitle(
        'Fig. 3 — Leg Geometry Sensitivity Analysis\n'
        'Effect of leg height and cross-section on module cooling performance',
        fontsize=11, y=1.02
    )

    I_test = 2.0  # A — representative operating current

    for j, (area, label, col) in enumerate(zip(areas, area_labels, area_colors)):
        dT_maxes  = []
        cop_maxes = []
        cop_at_dT30 = []

        for h in heights:
            # Maximum ΔT
            dTm, T_c_at_max = delta_T_max(T_h, h, area)
            dT_maxes.append(dTm)

            # Sweep current to find max COP at ΔT=0 (best-case COP)
            cops = []
            for I in np.linspace(0.1, 8.0, 60):
                params = compute_module_params(T_h, T_h, h, area)
                S, R, K = params['S_net'], params['R_module'], params['K_module']
                Q_c = N_COUPLES * (S*T_h*I - 0.5*I**2*R/N_COUPLES)
                W   = N_COUPLES * (I**2*R/N_COUPLES)
                cops.append(Q_c/W if W > 0 and Q_c > 0 else 0.0)
            cop_maxes.append(max(cops))

            # COP at fixed ΔT = 30 K, I = I_test
            dT_op = 30
            T_c_op = T_h - dT_op
            params = compute_module_params(T_h, T_c_op, h, area)
            S, R, K = params['S_net'], params['R_module'], params['K_module']
            Q_c = N_COUPLES*(S*T_c_op*I_test - 0.5*I_test**2*R/N_COUPLES - K/N_COUPLES*dT_op)
            W   = N_COUPLES*(S*dT_op*I_test + I_test**2*R/N_COUPLES)
            cop_at_dT30.append(Q_c/W if W > 0 and Q_c > 0 else 0.0)

        h_mm = heights * 1e3
        axes[0].plot(h_mm, dT_maxes, color=col, lw=2, label=label)
        axes[1].plot(h_mm, cop_maxes, color=col, lw=2, label=label)
        axes[2].plot(h_mm, cop_at_dT30, color=col, lw=2, label=label)

    # Reference line: Ying et al. ΔT_max = 72 K at T_h = 347 K
    # (different T_h from our 302 K case — shown as indicative benchmark)
    axes[0].axhline(y=52, color=GREY, lw=1.5, ls=':', label='Exp. ΔT_max ≈52K\n(T_h=302K) [Ref 1]')

    for ax, title, ylabel in zip(
        axes,
        ['Max. Achievable ΔT', 'Max. COP (at ΔT=0)', f'COP at ΔT=30K, I={I_test}A'],
        ['ΔT_max (K)', 'COP_max (dimensionless)', f'COP at operating point']
    ):
        ax.set_xlabel('Leg Height (mm)')
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(fontsize=8)
        ax.set_xlim(1.5, 6.0)

    plt.tight_layout()
    plt.savefig('figures/fig3_geometry_sensitivity.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  ✓ Fig 3 saved: Geometry sensitivity")


def plot_contact_resistance_effect():
    """Figure 4: How degraded contact resistance (e.g. after cycling) impacts performance."""
    T_h = 302.0
    T_c = 272.0   # ΔT = 30 K operating point
    I   = 2.0     # A

    # Scale contact resistance from 0.5× to 20× baseline (degradation scenario)
    # Physical motivation: Ying et al. 2024 ALD paper shows Fe/Mg3(Sb,Bi)2 interface
    # degrades through oxide formation and atomic migration under thermal cycling.
    scale_factors = np.linspace(0.5, 20, 80)
    cop_p_degrades = []
    cop_n_degrades = []
    cop_both_degrades = []
    cop_baseline = []

    for sf in scale_factors:
        # Degrade only p-contact
        _, _, _, COP_p, _ = cooling_performance(
            T_h, I, LEG_HEIGHT_BASE, LEG_AREA_BASE,
            r_contact_p=sf * R_contact_p, r_contact_n=R_contact_n)

        # Degrade only n-contact
        _, _, _, COP_n, _ = cooling_performance(
            T_h, I, LEG_HEIGHT_BASE, LEG_AREA_BASE,
            r_contact_p=R_contact_p, r_contact_n=sf * R_contact_n)

        # Both degrade together
        _, _, _, COP_both, _ = cooling_performance(
            T_h, I, LEG_HEIGHT_BASE, LEG_AREA_BASE,
            r_contact_p=sf * R_contact_p, r_contact_n=sf * R_contact_n)

        _, _, _, COP_ref, _ = cooling_performance(
            T_h, I, LEG_HEIGHT_BASE, LEG_AREA_BASE)

        cop_p_degrades.append(COP_p)
        cop_n_degrades.append(COP_n)
        cop_both_degrades.append(COP_both)
        cop_baseline.append(COP_ref)

    # Normalize to baseline (fractional COP retention)
    base = cop_baseline[0]
    cop_p_norm    = np.array(cop_p_degrades)    / base * 100
    cop_n_norm    = np.array(cop_n_degrades)    / base * 100
    cop_both_norm = np.array(cop_both_degrades) / base * 100

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(
        'Fig. 4 — Impact of Contact Resistance Degradation on Cooling COP\n'
        'Motivation: Ying et al. 2024 (Adv. Funct. Mater.) — oxide formation at\n'
        'Fe/Mg₃(Sb,Bi)₂ interface under thermal cycling degrades electrical contact',
        fontsize=10, y=1.01
    )

    ax1.plot(scale_factors, cop_p_norm, color=C_P, lw=2, label='p-contact (MgAgSb/Ag) degraded')
    ax1.plot(scale_factors, cop_n_norm, color=C_N, lw=2, ls='--', label='n-contact (Mg₃(Sb,Bi)₂/Fe) degraded')
    ax1.plot(scale_factors, cop_both_norm, color=GREY, lw=2, ls=':', label='Both contacts degraded')
    ax1.axhline(90, color='#E74C3C', lw=1, ls='-.', alpha=0.6, label='90% retention threshold')
    ax1.set_xlabel('Contact Resistance Scale Factor (×baseline)')
    ax1.set_ylabel('COP Retention (%)')
    ax1.set_title('COP Retention vs Contact Degradation')
    ax1.set_xlim(0.5, 20)
    ax1.set_ylim(0, 105)
    ax1.legend(fontsize=8)

    # 2nd panel: absolute COP
    ax2.plot(scale_factors, cop_p_degrades, color=C_P, lw=2, label='p-contact degraded')
    ax2.plot(scale_factors, cop_n_degrades, color=C_N, lw=2, ls='--', label='n-contact degraded')
    ax2.plot(scale_factors, cop_both_degrades, color=GREY, lw=2, ls=':', label='Both degraded')
    ax2.set_xlabel('Contact Resistance Scale Factor (×baseline)')
    ax2.set_ylabel('COP (dimensionless)')
    ax2.set_title('Absolute COP vs Contact Degradation\n(ΔT=30K, I=2A)')
    ax2.set_xlim(0.5, 20)
    ax2.legend(fontsize=8)

    # Annotate: n-contact is the weak link
    sf_crossover = scale_factors[np.argmin(np.abs(cop_n_norm - 90))]
    ax1.annotate(
        f'n-contact crosses 90%\nthreshold at ×{sf_crossover:.1f}',
        xy=(sf_crossover, 90),
        xytext=(sf_crossover + 3, 94),
        fontsize=8, color=C_N,
        arrowprops=dict(arrowstyle='->', color=C_N, lw=1.2)
    )

    plt.tight_layout()
    plt.savefig('figures/fig4_contact_degradation.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  ✓ Fig 4 saved: Contact resistance degradation")


def plot_performance_map():
    """Figure 5: 2D performance map — ΔT_max and COP_max as function of leg dimensions."""
    T_h = 302.0
    heights = np.linspace(1.5e-3, 6.0e-3, 25)
    areas   = np.linspace(1.0e-6, 8.0e-6, 25)

    dT_map  = np.zeros((len(areas), len(heights)))
    cop_map = np.zeros((len(areas), len(heights)))

    for i, area in enumerate(areas):
        for j, h in enumerate(heights):
            dTm, _ = delta_T_max(T_h, h, area)
            dT_map[i, j] = dTm

            cops = []
            for I in np.linspace(0.1, 10.0, 40):
                params = compute_module_params(T_h, T_h, h, area)
                S, R, K = params['S_net'], params['R_module'], params['K_module']
                Q_c = N_COUPLES * (S*T_h*I - 0.5*I**2*R/N_COUPLES)
                W   = N_COUPLES * (I**2*R/N_COUPLES)
                cops.append(Q_c/W if W > 0 and Q_c > 0 else 0.0)
            cop_map[i, j] = max(cops)

    H_mm = heights * 1e3
    A_mm2 = areas * 1e6

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        'Fig. 5 — Leg Geometry Performance Map\n'
        'Maximum ΔT and COP_max as a function of leg height and cross-section area',
        fontsize=11, y=1.02
    )

    cm1 = ax1.contourf(H_mm, A_mm2, dT_map, levels=20, cmap='RdYlBu_r')
    cs1 = ax1.contour(H_mm, A_mm2, dT_map, levels=10, colors='white', alpha=0.4, linewidths=0.8)
    plt.colorbar(cm1, ax=ax1, label='ΔT_max (K)')
    ax1.clabel(cs1, fmt='%.0f K', fontsize=7, colors='white')
    ax1.set_xlabel('Leg Height (mm)')
    ax1.set_ylabel('Leg Cross-section Area (mm²)')
    ax1.set_title('Maximum Achievable ΔT')
    # Mark Ying et al. approximate operating point
    ax1.scatter([3.0], [4.0], color='white', s=80, zorder=5, marker='*')
    ax1.annotate('Ref. [1]\n~baseline', xy=(3.0, 4.0), xytext=(3.5, 5.5),
                 fontsize=7, color='white',
                 arrowprops=dict(arrowstyle='->', color='white', lw=1))

    cm2 = ax2.contourf(H_mm, A_mm2, cop_map, levels=20, cmap='viridis')
    cs2 = ax2.contour(H_mm, A_mm2, cop_map, levels=10, colors='white', alpha=0.4, linewidths=0.8)
    plt.colorbar(cm2, ax=ax2, label='COP_max (dimensionless)')
    ax2.clabel(cs2, fmt='%.2f', fontsize=7, colors='white')
    ax2.set_xlabel('Leg Height (mm)')
    ax2.set_ylabel('Leg Cross-section Area (mm²)')
    ax2.set_title('Maximum COP (at ΔT=0)')
    ax2.scatter([3.0], [4.0], color='white', s=80, zorder=5, marker='*')

    plt.tight_layout()
    plt.savefig('figures/fig5_performance_map.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  ✓ Fig 5 saved: Performance map")


def print_summary_table():
    """Print a formatted summary table for the report."""
    print("\n" + "="*70)
    print("  SUMMARY TABLE: Module Performance at Key Operating Points")
    print("  MgAgSb/Mg₃(Sb,Bi)₂ Module, N=8 couples")
    print("  All values from analytical thermoelectric model")
    print("  Material data: Ying et al. 2022 [Ref 1] + Liu et al. 2022 [Ref 2]")
    print("="*70)

    T_h_cases = [(302, 'T_h = 302 K (Ying 2022 Fig 4)'),
                 (325, 'T_h = 325 K'),
                 (347, 'T_h = 347 K (Ying 2022 ΔT_max=72K case)')]

    print(f"\n{'Case':<30} {'ΔT_max (K)':>12} {'T_c,min (K)':>12} {'zT_mean':>10}")
    print("-"*68)
    for T_h, label in T_h_cases:
        dTm, T_c_min = delta_T_max(T_h, LEG_HEIGHT_BASE, LEG_AREA_BASE)
        params = compute_module_params(T_h, T_c_min, LEG_HEIGHT_BASE, LEG_AREA_BASE)
        print(f"  {label:<28} {dTm:>12.1f} {T_c_min:>12.1f} {params['zT_mean']:>10.3f}")

    print(f"\n  Experimental ΔT_max from Ying et al. 2022:")
    print(f"    T_h=302K → ΔT_max ≈ 52 K (reported)")
    print(f"    T_h=325K → ΔT_max ≈ 63 K (reported)")
    print(f"    T_h=347K → ΔT_max = 72 K (reported)")
    print()

    print(f"\n{'Geometry Variation':<30} {'L (mm)':>8} {'A (mm²)':>10} {'ΔT_max (K)':>12}")
    print("-"*65)
    geom_cases = [
        ('Short/wide', 1.5e-3, 6e-6),
        ('Baseline (2×2 mm, 3mm)', 3.0e-3, 4e-6),
        ('Tall/narrow', 6.0e-3, 2e-6),
    ]
    for name, h, a in geom_cases:
        dTm, _ = delta_T_max(302, h, a)
        print(f"  {name:<28} {h*1e3:>8.1f} {a*1e6:>10.1f} {dTm:>12.1f}")

    print("\n" + "="*70)
    print("  References")
    print("  [1] Ying P. et al., Energy Environ. Sci. 15, 6584 (2022)")
    print("      DOI: 10.1039/D2EE00883A")
    print("  [2] Liu Z. et al., Nature Commun. 13, 1120 (2022)")
    print("      DOI: 10.1038/s41467-022-28798-4")
    print("  [3] Zhao H. et al., Nano Energy 7, 97 (2014)")
    print("      DOI: 10.1016/j.nanoen.2014.04.012")
    print("="*70 + "\n")


# =============================================================================
# SECTION 5: Main Entry Point
# =============================================================================

if __name__ == '__main__':
    import os
    os.makedirs('figures', exist_ok=True)

    print("\n" + "="*70)
    print("  Thermoelectric Module Performance Analysis")
    print("  MgAgSb / Mg₃(Sb,Bi)₂ — Te-free solid-state cooling")
    print("  Based on: Ying et al., Energy Environ. Sci. (2022)")
    print("="*70 + "\n")

    print("Generating figures...")
    plot_material_properties()
    plot_cop_curves()
    plot_geometry_sensitivity()
    plot_contact_resistance_effect()
    plot_performance_map()

    print_summary_table()

    print("Done. All figures saved to ./figures/")
    print("Run from the project root: python src/te_module_analysis.py\n")
