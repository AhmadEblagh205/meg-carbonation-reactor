"""Illustrative non-isothermal plug-flow reactor for EO carbonation.

    ethylene oxide (EO) + CO2 -> ethylene carbonate (EC)

This is a clean reconstruction of the reactor section from a UCL group design
project. The model is intentionally *illustrative*: kinetics and thermo values
are not a validated industrial parameter set. The value of the repository is
in the balances, numerical checks, and sensitivity analysis rather than in
claiming plant prediction.

Assumptions
-----------
- steady-state pseudo-homogeneous plug flow
- constant volumetric flow, density*heat-capacity, heat of reaction and UA/V
- second-order rate: r = k(T) C_EO C_CO2
- Arrhenius k(T)
- constant coolant temperature
- one irreversible reaction with 1:1:1 stoichiometry

The state is reaction extent plus temperature. Using extent makes exact
stoichiometry explicit and lets the solver stop reaction cleanly when the
limiting reactant is exhausted.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp

R_GAS = 8.31446261815324  # J/mol/K


@dataclass(frozen=True)
class ReactorParams:
    # Feed
    F_EO0: float = 100.0          # mol/s
    ratio_CO2: float = 1.1        # CO2:EO molar feed ratio
    T0: float = 380.0             # inlet temperature, K
    v_flow: float = 0.05          # m3/s, assumed constant

    # Kinetics (illustrative)
    A: float = 2.0e4              # m3/(mol s)
    Ea: float = 68_000.0          # J/mol

    # Thermo / heat transfer (illustrative)
    dH: float = -95_000.0         # J/mol EO, exothermic
    rho_cp: float = 2.4e6         # J/(m3 K)
    Ua: float = 60_000.0          # W/(m3 K)
    T_cool: float = 375.0         # K

    # Geometry
    V_total: float = 8.0          # m3

    @property
    def F_CO20(self) -> float:
        return self.ratio_CO2 * self.F_EO0

    @property
    def limiting_extent(self) -> float:
        """Maximum physically possible reaction extent, mol/s."""
        return min(self.F_EO0, self.F_CO20)

    @property
    def heat_capacity_rate(self) -> float:
        """rho*Cp*volumetric flow, W/K."""
        return self.rho_cp * self.v_flow


def validate_params(p: ReactorParams) -> None:
    positives = {
        "F_EO0": p.F_EO0,
        "ratio_CO2": p.ratio_CO2,
        "T0": p.T0,
        "v_flow": p.v_flow,
        "A": p.A,
        "Ea": p.Ea,
        "rho_cp": p.rho_cp,
        "Ua": p.Ua,
        "T_cool": p.T_cool,
        "V_total": p.V_total,
    }
    bad = [name for name, value in positives.items() if not np.isfinite(value) or value <= 0]
    if bad:
        raise ValueError(f"Parameters must be finite and positive: {', '.join(bad)}")
    if not np.isfinite(p.dH) or p.dH >= 0:
        raise ValueError("dH must be finite and negative for the exothermic model")


def arrhenius_k(T: float | np.ndarray, p: ReactorParams) -> float | np.ndarray:
    """Second-order kinetic constant in m3/(mol s)."""
    return p.A * np.exp(-p.Ea / (R_GAS * np.asarray(T)))


def rate(C_EO: float, C_CO2: float, T: float, p: ReactorParams) -> float:
    """Reaction rate, mol/(m3 s), with non-negative concentrations."""
    k = float(arrhenius_k(T, p))
    return k * max(C_EO, 0.0) * max(C_CO2, 0.0)


def _rhs(V: float, y: np.ndarray, p: ReactorParams) -> list[float]:
    """d/dV of [extent, temperature]."""
    xi, T = y
    F_EO = max(p.F_EO0 - xi, 0.0)
    F_CO2 = max(p.F_CO20 - xi, 0.0)
    C_EO = F_EO / p.v_flow
    C_CO2 = F_CO2 / p.v_flow
    r = rate(C_EO, C_CO2, T, p)
    dxi_dV = r
    dT_dV = ((-p.dH) * r - p.Ua * (T - p.T_cool)) / p.heat_capacity_rate
    return [dxi_dV, dT_dV]


def _depletion_event(V: float, y: np.ndarray, p: ReactorParams) -> float:
    return p.limiting_extent - y[0]


_depletion_event.terminal = True
_depletion_event.direction = -1


def simulate(p: ReactorParams, n_points: int = 300) -> dict[str, np.ndarray | float | bool]:
    """Integrate the PFR and return profiles plus balance diagnostics.

    If a limiting reactant is exhausted before the specified reactor volume,
    reaction stops and the remaining length is cooling-only. This avoids small
    negative molar flows from numerical overshoot in aggressive scenarios.
    """
    validate_params(p)
    if n_points < 2:
        raise ValueError("n_points must be at least 2")

    sol = solve_ivp(
        _rhs,
        (0.0, p.V_total),
        [0.0, p.T0],
        args=(p,),
        events=_depletion_event,
        dense_output=True,
        max_step=p.V_total / 300.0,
        rtol=1e-9,
        atol=[1e-9, 1e-9],
    )
    if not sol.success:
        raise RuntimeError(f"PFR integration failed: {sol.message}")

    V = np.linspace(0.0, p.V_total, n_points)
    depleted = bool(sol.t_events and len(sol.t_events[0]))

    if depleted:
        V_dep = float(sol.t_events[0][0])
        xi_dep, T_dep = sol.sol(V_dep)
        pre = V <= V_dep
        xi = np.empty_like(V)
        T = np.empty_like(V)
        xi[pre], T[pre] = sol.sol(V[pre])
        xi[~pre] = p.limiting_extent
        cooling_rate = p.Ua / p.heat_capacity_rate  # 1/m3
        T[~pre] = p.T_cool + (float(T_dep) - p.T_cool) * np.exp(
            -cooling_rate * (V[~pre] - V_dep)
        )
    else:
        V_dep = np.nan
        xi, T = sol.sol(V)

    # Remove sub-numerical event overshoot before deriving flows.
    xi = np.clip(np.asarray(xi, dtype=float), 0.0, p.limiting_extent)
    T = np.asarray(T, dtype=float)
    F_EO = p.F_EO0 - xi
    F_CO2 = p.F_CO20 - xi
    F_EC = xi.copy()
    X = xi / p.F_EO0

    # Independent energy-balance audit: integrate heat removed by the wall and
    # compare with generated reaction heat minus sensible enthalpy accumulation.
    q_wall_density = p.Ua * (T - p.T_cool)          # W/m3
    Q_removed = cumulative_trapezoid(q_wall_density, V, initial=0.0)  # W
    Q_generated = (-p.dH) * xi                      # W
    Q_sensible = p.heat_capacity_rate * (T - p.T0) # W
    energy_residual = Q_generated - Q_removed - Q_sensible
    scale = np.maximum(np.maximum(np.abs(Q_generated), np.abs(Q_removed)), 1.0)
    energy_residual_rel = energy_residual / scale

    return {
        "V": V,
        "extent": xi,
        "F_EO": F_EO,
        "F_CO2": F_CO2,
        "F_EC": F_EC,
        "T": T,
        "X": X,
        "Q_generated_W": Q_generated,
        "Q_removed_W": Q_removed,
        "Q_sensible_W": Q_sensible,
        "energy_residual_W": energy_residual,
        "energy_residual_rel": energy_residual_rel,
        "depleted": depleted,
        "depletion_volume_m3": V_dep,
    }


def outlet_conversion(p: ReactorParams) -> float:
    return float(simulate(p, n_points=200)["X"][-1])


def peak_temperature(p: ReactorParams) -> float:
    return float(np.max(simulate(p)["T"]))


def full_conversion_adiabatic_rise(p: ReactorParams) -> float:
    """Maximum reaction-only temperature rise if the limiting reactant fully reacts.

    This is a reference scale, not a universal upper bound when the coolant is
    hotter than the feed.
    """
    return (-p.dH) * p.limiting_extent / p.heat_capacity_rate


def thermal_ntu(p: ReactorParams) -> float:
    """Dimensionless wall heat-removal strength Ua*V/(rho*Cp*v_flow)."""
    return p.Ua * p.V_total / p.heat_capacity_rate
