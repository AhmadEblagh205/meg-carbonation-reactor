"""Map the conversion / peak-temperature trade-off without inventing a safety limit."""
from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from reactor import ReactorParams, simulate

ROOT = Path(__file__).resolve().parents[1]


def operating_grid(base: ReactorParams | None = None):
    base = base or ReactorParams()
    T0s = np.linspace(360.0, 400.0, 17)
    ratios = np.linspace(1.0, 1.6, 13)
    conversion = np.zeros((len(T0s), len(ratios)))
    peak = np.zeros_like(conversion)
    for i, T0 in enumerate(T0s):
        for j, ratio in enumerate(ratios):
            p = replace(base, T0=float(T0), ratio_CO2=float(ratio))
            res = simulate(p, n_points=140)
            conversion[i, j] = float(res["X"][-1])
            peak[i, j] = float(np.max(res["T"]))
    return T0s, ratios, conversion, peak


def pareto_mask(conversion: np.ndarray, peak: np.ndarray) -> np.ndarray:
    """Nondominated points when maximising conversion and minimising peak T."""
    c = conversion.ravel()
    t = peak.ravel()
    keep = np.ones(len(c), dtype=bool)
    for i in range(len(c)):
        # j dominates i if it is no worse on both objectives and strictly better on at least one.
        dom = (c >= c[i]) & (t <= t[i]) & ((c > c[i]) | (t < t[i]))
        dom[i] = False
        if np.any(dom):
            keep[i] = False
    return keep.reshape(conversion.shape)


def main():
    T0s, ratios, conversion, peak = operating_grid()
    pareto = pareto_mask(conversion, peak)
    TT, RR = np.meshgrid(T0s, ratios, indexing="ij")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    c1 = axes[0].contourf(ratios, T0s, 100 * conversion, levels=16)
    axes[0].set(xlabel="CO2:EO feed ratio", ylabel="Inlet temperature (K)", title="Outlet EO conversion (%)")
    fig.colorbar(c1, ax=axes[0], label="Conversion (%)")

    c2 = axes[1].contourf(ratios, T0s, peak - 273.15, levels=16)
    axes[1].set(xlabel="CO2:EO feed ratio", ylabel="Inlet temperature (K)", title="Peak temperature (°C)")
    fig.colorbar(c2, ax=axes[1], label="Peak temperature (°C)")

    axes[2].scatter(peak.ravel() - 273.15, 100 * conversion.ravel(), s=18, alpha=0.45, label="grid")
    axes[2].scatter((peak[pareto] - 273.15), 100 * conversion[pareto], s=30, label="nondominated")
    axes[2].set(xlabel="Peak temperature (°C)", ylabel="Outlet conversion (%)",
                title="Conversion / thermal-state trade-off")
    axes[2].grid(alpha=0.25)
    axes[2].legend(fontsize=8)

    out = ROOT / "reports" / "figures" / "operating_window.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)

    idx_low_peak = np.unravel_index(np.argmin(peak), peak.shape)
    idx_high_conv = np.unravel_index(np.argmax(conversion), conversion.shape)
    summary = {
        "grid_points": int(conversion.size),
        "pareto_points": int(pareto.sum()),
        "highest_conversion_pct": round(float(100 * conversion[idx_high_conv]), 3),
        "highest_conversion_condition": {
            "T0_K": float(T0s[idx_high_conv[0]]),
            "CO2_to_EO": float(ratios[idx_high_conv[1]]),
            "peak_temperature_C": round(float(peak[idx_high_conv] - 273.15), 3),
        },
        "lowest_peak_condition": {
            "T0_K": float(T0s[idx_low_peak[0]]),
            "CO2_to_EO": float(ratios[idx_low_peak[1]]),
            "conversion_pct": round(float(100 * conversion[idx_low_peak]), 3),
            "peak_temperature_C": round(float(peak[idx_low_peak] - 273.15), 3),
        },
        "interpretation": (
            "The map exposes the conversion-versus-thermal-state trade-off without inventing a safety threshold. "
            "A real design decision would need validated kinetics, materials limits and dynamic safety analysis."
        ),
    }
    (ROOT / "reports" / "operating_window_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
