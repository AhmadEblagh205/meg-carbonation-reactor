"""Thermal-sensitivity screen for the illustrative carbonation PFR.

This is deliberately *not* called a runaway/stability calculation. A steady
PFR sensitivity sweep cannot establish dynamic runaway, multiplicity, or a
stability boundary. It can, however, measure thermal amplification: how much
peak reactor temperature changes per kelvin of inlet-temperature change.
"""
from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from reactor import ReactorParams, peak_temperature

ROOT = Path(__file__).resolve().parents[1]


def peak_curve(T0s: np.ndarray, T_cool: float, base: ReactorParams | None = None) -> np.ndarray:
    base = base or ReactorParams()
    return np.array([
        peak_temperature(replace(base, T0=float(T0), T_cool=float(T_cool)))
        for T0 in T0s
    ])


def thermal_amplification(T0s: np.ndarray, T_cool: float, base: ReactorParams | None = None):
    peaks = peak_curve(T0s, T_cool, base)
    amp = np.gradient(peaks, T0s)  # dT_peak / dT_in
    return peaks, amp


def build_screen(base: ReactorParams | None = None):
    base = base or ReactorParams()
    T0s = np.linspace(370.0, 400.0, 13)
    coolants = [365.0, 370.0, 375.0, 380.0]
    records = {}
    for tc in coolants:
        peaks, amp = thermal_amplification(T0s, tc, base)
        i = int(np.argmax(amp))
        records[str(int(tc))] = {
            "max_dTpeak_dTin": float(amp[i]),
            "at_inlet_K": float(T0s[i]),
            "peak_temperature_K_at_that_point": float(peaks[i]),
            "peak_temperature_K_at_400K_feed": float(peaks[-1]),
        }
    return T0s, coolants, records


def main():
    base = ReactorParams()
    T0s, coolants, records = build_screen(base)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.6))
    for tc in coolants:
        peaks, amp = thermal_amplification(T0s, tc, base)
        ax1.plot(T0s, peaks, marker="o", ms=3.5, lw=1.6, label=f"coolant {tc:.0f} K")
        ax2.plot(T0s, amp, marker="o", ms=3.5, lw=1.6, label=f"coolant {tc:.0f} K")

    ax1.plot(T0s, T0s, ls="--", lw=1, color="grey", label="peak = inlet")
    ax1.set(xlabel="Inlet temperature (K)", ylabel="Peak reactor temperature (K)",
            title="Peak temperature state")
    ax1.grid(alpha=0.25)
    ax1.legend(fontsize=8)

    ax2.axhline(1.0, ls="--", lw=1, color="grey")
    ax2.set(xlabel="Inlet temperature (K)", ylabel=r"Thermal amplification $dT_{peak}/dT_{in}$",
            title="Sensitivity, not a dynamic-runaway criterion")
    ax2.grid(alpha=0.25)

    out = ROOT / "reports" / "figures" / "thermal_amplification.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)

    strongest = max(records.items(), key=lambda kv: kv[1]["max_dTpeak_dTin"])
    warmest_at_400 = max(records.items(), key=lambda kv: kv[1]["peak_temperature_K_at_400K_feed"])
    summary = {
        "definition": "finite-difference dT_peak/dT_in over a steady PFR parameter sweep",
        "not_a_claim_of": ["dynamic stability", "thermal runaway boundary", "reactor multiplicity"],
        "by_coolant_K": {k: {kk: round(vv, 4) for kk, vv in v.items()} for k, v in records.items()},
        "strongest_amplification_coolant_K": int(strongest[0]),
        "warmest_peak_at_400K_feed_coolant_K": int(warmest_at_400[0]),
        "finding": (
            "The coolant giving the highest absolute peak temperature is not the same as the coolant "
            "giving the strongest local amplification. Thermal state and thermal sensitivity are different diagnostics."
        ),
    }
    (ROOT / "reports" / "thermal_sensitivity_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
