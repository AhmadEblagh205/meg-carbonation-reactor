"""Run the illustrative EO carbonation PFR and generate core diagnostics.

Usage:
    python src/run.py
    python src/run.py --T0 390 --ratio 1.3
"""
from __future__ import annotations

import argparse
import json
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from reactor import (
    ReactorParams,
    full_conversion_adiabatic_rise,
    outlet_conversion,
    peak_temperature,
    simulate,
    thermal_ntu,
)

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "reports" / "figures"


def profiles_figure(res, p):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.5, 7), sharex=True)
    ax1.plot(res["V"], res["X"] * 100, lw=2)
    ax1.set(ylabel="EO conversion (%)", title="Axial profiles along the carbonation reactor")
    ax1.grid(alpha=0.3)
    ax2.plot(res["V"], res["T"] - 273.15, lw=2, label="reaction mixture")
    ax2.axhline(p.T_cool - 273.15, ls="--", lw=1.2, label="coolant")
    ax2.set(xlabel="Reactor volume (m$^3$)", ylabel="Temperature (°C)")
    ax2.legend()
    ax2.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGURES / "axial_profiles.png", dpi=160)
    plt.close(fig)


def sensitivity_figure(p):
    T0s = np.linspace(360, 400, 9)
    conv_T, peak_T = [], []
    for T0 in T0s:
        q = replace(p, T0=float(T0))
        conv_T.append(outlet_conversion(q) * 100)
        peak_T.append(peak_temperature(q) - 273.15)

    ratios = np.linspace(1.0, 1.6, 9)
    conv_r = [outlet_conversion(replace(p, ratio_CO2=float(r))) * 100 for r in ratios]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))
    ax1.plot(T0s - 273.15, conv_T, "o-", lw=1.8)
    ax1.set(xlabel="Inlet temperature (°C)", ylabel="Outlet conversion (%)",
            title="Sensitivity to inlet temperature")
    ax1.grid(alpha=0.3)
    ax1b = ax1.twinx()
    ax1b.plot(T0s - 273.15, peak_T, "s--", lw=1.5)
    ax1b.set_ylabel("Peak temperature (°C)")

    ax2.plot(ratios, conv_r, "o-", lw=1.8)
    ax2.set(xlabel="CO$_2$:EO feed ratio", ylabel="Outlet conversion (%)",
            title="Sensitivity to CO$_2$ excess")
    ax2.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGURES / "sensitivity.png", dpi=160)
    plt.close(fig)


def energy_audit_figure(res, p):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.5))
    scale = 1e-6
    ax1.plot(res["V"], scale * res["Q_generated_W"], label="reaction heat generated")
    ax1.plot(res["V"], scale * res["Q_removed_W"], label="wall heat removed")
    ax1.plot(res["V"], scale * res["Q_sensible_W"], label="sensible enthalpy rise")
    ax1.set(xlabel="Reactor volume (m$^3$)", ylabel="Cumulative power-equivalent (MW)",
            title="Integrated energy terms")
    ax1.grid(alpha=0.25)
    ax1.legend(fontsize=8)

    ax2.semilogy(res["V"], np.maximum(np.abs(res["energy_residual_rel"]), 1e-16))
    ax2.set(xlabel="Reactor volume (m$^3$)", ylabel="Relative balance residual",
            title="Independent energy-balance audit")
    ax2.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIGURES / "energy_balance_audit.png", dpi=160)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--T0", type=float, default=None, help="inlet temperature in K")
    ap.add_argument("--ratio", type=float, default=None, help="CO2:EO molar ratio")
    args = ap.parse_args()

    p = ReactorParams()
    if args.T0 is not None:
        p = replace(p, T0=args.T0)
    if args.ratio is not None:
        p = replace(p, ratio_CO2=args.ratio)

    FIGURES.mkdir(parents=True, exist_ok=True)
    res = simulate(p, n_points=600)
    profiles_figure(res, p)
    sensitivity_figure(p)
    energy_audit_figure(res, p)

    adiabatic_rise = full_conversion_adiabatic_rise(p)
    actual_peak_rise = float(np.max(res["T"]) - p.T0)
    max_energy_resid = float(np.max(np.abs(res["energy_residual_rel"])))
    summary = {
        "outlet_conversion_pct": round(float(res["X"][-1]) * 100, 3),
        "peak_temperature_C": round(float(np.max(res["T"])) - 273.15, 3),
        "outlet_temperature_C": round(float(res["T"][-1]) - 273.15, 3),
        "EC_production_mol_s": round(float(res["F_EC"][-1]), 3),
        "full_limiting_reactant_adiabatic_rise_K": round(float(adiabatic_rise), 3),
        "actual_peak_rise_K": round(actual_peak_rise, 3),
        "peak_rise_fraction_of_reaction_only_reference": round(actual_peak_rise / adiabatic_rise, 4),
        "thermal_NTU": round(float(thermal_ntu(p)), 4),
        "max_relative_energy_balance_residual": max_energy_resid,
        "limiting_reactant_depleted_inside_reactor": bool(res["depleted"]),
        "note": "Illustrative kinetics/thermodynamics; not a validated industrial reactor prediction.",
    }
    out = ROOT / "reports" / "reactor_summary.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"Figures -> {FIGURES}")


if __name__ == "__main__":
    main()
