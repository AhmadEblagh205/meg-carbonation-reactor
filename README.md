# Non-Isothermal Carbonation Reactor — Thermal Sensitivity & Operating Window

[![CI](https://github.com/AhmadEblagh205/meg-carbonation-reactor/actions/workflows/ci.yml/badge.svg)](https://github.com/AhmadEblagh205/meg-carbonation-reactor/actions/workflows/ci.yml)
![Python](https://img.shields.io/badge/python-3.10%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)

This is a tested reconstruction of the ethylene-oxide carbonation reactor from my UCL mono-ethylene-glycol group design project. The full plant was a team exercise, and the reactor modelling here was my individual section. The kinetics and thermodynamic values are illustrative placeholders, so the model is not a validated industrial reactor model or a plant prediction.

## What it does

- Solves coupled mole and energy balances for a wall-cooled, pseudo-homogeneous PFR with a single second-order Arrhenius reaction (EO + CO2 to ethylene carbonate).
- Audits the integrated energy balance by comparing reaction heat, wall heat removal and sensible enthalpy rise.
- Computes dT_peak/dT_in as a thermal-amplification diagnostic across coolant and feed temperatures.
- Maps 221 combinations of feed temperature and CO2:EO ratio to show the conversion and peak-temperature trade-off.
- Runs a 128-point Latin-hypercube scenario screen over six model assumptions.

## Results

All values are for the illustrative parameters.

| Quantity | Value |
|---|---:|
| EO conversion | 92.36% |
| Peak temperature | 142.84 °C |
| Outlet temperature | 105.15 °C |
| EC production | 92.36 mol/s |
| Reaction-only full-conversion adiabatic rise | 79.17 K |
| Actual peak rise above feed | 35.99 K |
| Thermal NTU, `Ua V / (rho Cp v)` | 4.00 |

Stronger cooling lowers the absolute peak temperature but can raise the local sensitivity to feed temperature. Coolant at 365 K gives the largest amplification, 1.83 K/K near a 392.5 K feed, while at a 400 K feed coolant at 380 K gives the hottest peak, 448.5 K. In the operating-window grid the highest-conversion corner approaches 100% at a peak of about 188 °C, and the lowest-peak corner is about 123 °C at 81.4% conversion.

In the scenario screen, conversion is most strongly tied to CO2 excess (Spearman rho about +0.59), reference reaction rate (+0.57) and heat-transfer strength (-0.51). The ranges are chosen scenarios, not probability distributions. Figures are in `reports/figures/`.

## How to run

```bash
pip install -r requirements.txt

python src/run.py
python src/thermal_sensitivity.py
python src/operating_window.py
python src/scenario_sensitivity.py
pytest -q
```

`python src/runaway.py` still works and now runs the thermal-sensitivity screen. It does not compute dynamic runaway.

## Notes

The model assumes constant volumetric flow, constant rho*Cp, constant heat of reaction, a constant wall heat-transfer coefficient per reactor volume, a constant coolant temperature and one irreversible reaction. It leaves out gas-liquid mass transfer, phase equilibrium, catalyst deactivation, pressure drop and reactor dynamics. A steady-state sweep cannot establish dynamic stability, so I describe the sweep as thermal sensitivity. The balances, numerical safeguards and results for each analysis are in `docs/METHODOLOGY.md`.

## Licence

MIT, see [LICENSE](LICENSE).
