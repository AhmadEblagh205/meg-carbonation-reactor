"""Global *scenario* sensitivity for illustrative model assumptions.

The ranges below are deliberate stress-test ranges, not calibrated parameter
uncertainty distributions. Latin-hypercube sampling plus Spearman rank
correlations answers a limited question: which assumptions drive model outputs
most strongly across this chosen scenario box?
"""
from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import qmc, spearmanr

from reactor import R_GAS, ReactorParams, arrhenius_k, simulate

ROOT = Path(__file__).resolve().parents[1]


def build_scenarios(n: int = 128, seed: int = 205):
    base = ReactorParams()
    sampler = qmc.LatinHypercube(d=6, seed=seed)
    u = sampler.random(n)

    # Independent scenario axes. Ea is varied while A is re-anchored so that
    # k(380 K) is controlled separately by k_ref_multiplier.
    T_ref = 380.0
    k_ref_base = float(arrhenius_k(T_ref, base))
    k_mult = np.exp(np.log(0.5) + u[:, 0] * (np.log(2.0) - np.log(0.5)))
    Ea = 60_000.0 + u[:, 1] * 16_000.0
    Ua_mult = 0.5 + u[:, 2] * 1.0
    dH_mult = 0.85 + u[:, 3] * 0.30
    T0 = 365.0 + u[:, 4] * 30.0
    ratio = 1.0 + u[:, 5] * 0.50

    names = ["k_ref_multiplier", "Ea_J_mol", "Ua_multiplier", "dH_magnitude_multiplier", "T0_K", "CO2_to_EO"]
    X = np.column_stack([k_mult, Ea, Ua_mult, dH_mult, T0, ratio])
    conv = np.empty(n)
    peak = np.empty(n)

    for i in range(n):
        k_ref = k_ref_base * k_mult[i]
        A = k_ref * np.exp(Ea[i] / (R_GAS * T_ref))
        p = replace(
            base,
            A=float(A),
            Ea=float(Ea[i]),
            Ua=float(base.Ua * Ua_mult[i]),
            dH=float(base.dH * dH_mult[i]),
            T0=float(T0[i]),
            ratio_CO2=float(ratio[i]),
        )
        res = simulate(p, n_points=140)
        conv[i] = float(res["X"][-1])
        peak[i] = float(np.max(res["T"]))

    return names, X, conv, peak


def rank_correlations(names, X, y):
    return {name: float(spearmanr(X[:, i], y).statistic) for i, name in enumerate(names)}


def main():
    names, X, conv, peak = build_scenarios()
    rho_conv = rank_correlations(names, X, conv)
    rho_peak = rank_correlations(names, X, peak)

    display = {
        "k_ref_multiplier": "rate at 380 K",
        "Ea_J_mol": "activation energy",
        "Ua_multiplier": "heat-transfer strength",
        "dH_magnitude_multiplier": "reaction heat magnitude",
        "T0_K": "inlet temperature",
        "CO2_to_EO": "CO2:EO ratio",
    }
    order = np.arange(len(names))
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), sharey=True)
    axes[0].barh(order, [rho_conv[n] for n in names])
    axes[0].set(xlabel="Spearman rank correlation", title="Drivers of conversion")
    axes[0].axvline(0, color="black", lw=0.8)
    axes[0].grid(axis="x", alpha=0.2)
    axes[1].barh(order, [rho_peak[n] for n in names])
    axes[1].set(xlabel="Spearman rank correlation", title="Drivers of peak temperature")
    axes[1].axvline(0, color="black", lw=0.8)
    axes[1].grid(axis="x", alpha=0.2)
    axes[0].set_yticks(order, labels=[display[n] for n in names])
    fig.suptitle("Scenario sensitivity — chosen ranges, not calibrated uncertainty")
    fig.tight_layout()
    out = ROOT / "reports" / "figures" / "scenario_sensitivity.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160)
    plt.close(fig)

    summary = {
        "n_latin_hypercube_scenarios": len(conv),
        "ranges_are": "illustrative stress-test ranges, not probability distributions",
        "spearman_conversion": {k: round(v, 4) for k, v in rho_conv.items()},
        "spearman_peak_temperature": {k: round(v, 4) for k, v in rho_peak.items()},
        "conversion_range_pct": [round(float(100 * conv.min()), 3), round(float(100 * conv.max()), 3)],
        "peak_temperature_range_C": [round(float(peak.min() - 273.15), 3), round(float(peak.max() - 273.15), 3)],
    }
    (ROOT / "reports" / "scenario_sensitivity_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
