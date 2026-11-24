# Changelog

## 0.3.0 — defensibility and thermal-sensitivity audit

- Reframed the old “runaway” sweep as a **thermal-sensitivity screen**. A steady-state PFR parameter sweep does not establish dynamic runaway, multiplicity or a stability boundary.
- Rewrote the reactor state in reaction-extent form so 1:1 stoichiometry is explicit and limiting-reactant depletion cannot create negative molar flows.
- Added parameter-domain and solver-success guards.
- Added post-depletion cooling through unused reactor volume.
- Added an independent integrated energy-balance audit.
- Added full-conversion reaction-only adiabatic temperature-rise and thermal-NTU reference scales.
- Added `dT_peak/dT_in` thermal-amplification maps. Stronger cooling lowers absolute peak temperature but can increase local amplification.
- Added a 221-point conversion/peak-temperature operating-window map with a nondominated Pareto set rather than an invented safety threshold.
- Added a 128-point Latin-hypercube **scenario sensitivity** study with reference-rate and activation-energy effects separated; ranges are explicitly not calibrated uncertainty distributions.
- Added `docs/METHODOLOGY.md` and expanded limitations/provenance language.
- Expanded tests from 4 to 15.

## 0.2.0 — update round 4

Added the original parametric hot-spot sweep over inlet and coolant temperature. Historical wording described the finite-difference slope as a runaway proxy; release 0.3.0 tightens that interpretation to thermal sensitivity.

## 0.1.0

Initial release.
