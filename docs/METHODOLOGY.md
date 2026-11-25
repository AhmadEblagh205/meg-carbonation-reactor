# Methodology and validation notes

## 1. Scope

This repository reconstructs the carbonation-reactor section of a UCL process-design project as a transparent, testable numerical model. Kinetic and thermodynamic parameters are illustrative. Conclusions are therefore about **model behaviour under stated assumptions**, not about a particular industrial reactor.

## 2. Reaction and state choice

The reaction is represented as

```text
EO + CO2 -> EC
```

with reaction extent `xi(V)` in mol/s. The molar flows are

```text
F_EO  = F_EO,0  - xi
F_CO2 = F_CO2,0 - xi
F_EC  = xi
```

so the stoichiometry is exact by construction.

The maximum physically possible extent is

```text
xi_max = min(F_EO,0, F_CO2,0)
```

and integration stops reaction at that limit. If the limiting reactant is exhausted before the end of the reactor, the remaining section is cooling-only.

## 3. Rate law

The illustrative second-order rate is

```text
r = k(T) C_EO C_CO2
k(T) = A exp(-Ea / RT)
```

with constant volumetric flow, so concentrations are `F_i / v_flow`.

`A` and `Ea` are placeholders rather than literature-fitted kinetic parameters.

## 4. PFR balances

The extent balance is

```text
dxi/dV = r
```

and the energy balance is

```text
dT/dV = [(-dH) r - Ua (T - T_cool)] / (rho Cp v_flow)
```

where `Ua` is heat-transfer conductance per reactor volume.

## 5. Energy-balance audit

At any reactor volume, cumulative reaction heat is

```text
Q_generated = (-dH) xi
```

and sensible enthalpy change is

```text
Q_sensible = (rho Cp v_flow) (T - T0)
```

Wall heat removal is evaluated independently by numerical integration:

```text
Q_removed = integral Ua (T - T_cool) dV
```

The residual

```text
Q_generated - Q_removed - Q_sensible
```

is used as a numerical audit rather than being forced to zero by definition.

## 6. Reference scales

Two simple engineering scales help interpret the base case.

### Reaction-only full-conversion adiabatic rise

```text
DeltaT_ad,ref = (-dH) xi_max / (rho Cp v_flow)
```

For the default parameters this is about 79.17 K. It is a reference scale, not a universal upper bound when the coolant is hotter than the feed.

### Thermal NTU

```text
NTU_th = Ua V_total / (rho Cp v_flow)
```

The default value is 4.0, indicating strong wall heat exchange relative to the sensible heat-capacity flow.

## 7. Thermal amplification

The steady-state screen computes

```text
dT_peak / dT_in
```

by finite differences over a feed-temperature grid. A value above one means the peak temperature rises by more than one kelvin for each kelvin increase in feed temperature over that local region.

This is **not** a dynamic stability criterion. The model has no transient energy storage, control loop, relief system or multiplicity calculation, so the analysis must not be described as establishing a thermal-runaway boundary.

## 8. Operating-window/Pareto analysis

A grid over `T0` and `CO2:EO` ratio is evaluated for two objectives:

- maximise outlet conversion;
- minimise peak reactor temperature.

A point is labelled nondominated when no other grid point gives at least as much conversion at no higher peak temperature, with a strict improvement in at least one objective.

No arbitrary safety threshold is imposed.

## 9. Global scenario sensitivity

The scenario screen uses a 128-point Latin hypercube over six deliberately chosen axes:

- reference rate at 380 K: 0.5x to 2x;
- activation energy: 60 to 76 kJ/mol;
- heat-transfer strength: 0.5x to 1.5x;
- reaction-heat magnitude: 0.85x to 1.15x;
- inlet temperature: 365 to 395 K;
- CO2:EO ratio: 1.0 to 1.5.

When activation energy is changed, `A` is re-anchored so the reference rate at 380 K is controlled separately. This avoids conflating baseline activity with Arrhenius temperature sensitivity.

Spearman rank correlation is then used as a monotonic global-sensitivity diagnostic.

The ranges are stress-test choices, not posterior distributions or measured parameter uncertainty. Therefore the resulting coefficients are scenario sensitivities, not probabilistic uncertainty attribution.

## 10. Main limitations

The current model omits several effects that can matter materially in a real EO/CO2 reactor:

- validated catalytic kinetics;
- gas-liquid mass transfer and solubility;
- non-ideal phase behaviour;
- pressure drop;
- property variation with composition/temperature;
- catalyst deactivation;
- side reactions/selectivity;
- reactor dynamics, multiplicity and control;
- equipment/material/relief constraints.

A defensible next step would be to replace the illustrative kinetics with published/experimental kinetics, then re-run the same balance, sensitivity and validation framework before making process claims.

## Notes moved from the README

Numerical safeguards added in the v2 audit: parameter-domain validation, explicit solver-success checking, a reaction-extent formulation for exact stoichiometry, handling of limiting-reactant depletion without negative outlet flows, post-depletion cooling through the remaining reactor volume, an independent integrated energy-balance residual, and tests that separate absolute thermal state from local thermal amplification.

Naming. An earlier version called the steady-state parameter sweep a runaway study. A static PFR sweep cannot establish dynamic stability, multiplicity or a runaway boundary, so the analysis now reports dT_peak/dT_in as a thermal-amplification diagnostic. The legacy command `python src/runaway.py` still runs and prints that it is the thermal-sensitivity screen.

Thermal sensitivity results from the controlled sweep. Coolant at 365 K gives the strongest local amplification, 1.83 K/K near a 392.5 K feed. At a 400 K feed, coolant at 380 K gives the hottest absolute peak, 448.5 K. Stronger cooling lowers the absolute peak but can raise the local sensitivity to feed temperature.

Operating window. The grid has 221 combinations of feed temperature and CO2:EO ratio. The highest-conversion corner approaches 100% conversion at a peak of about 188 C. The lowest-peak corner is about 123 C at 81.4% conversion. No "safe" temperature is defined, because that would need validated kinetics, material limits and dynamic safety analysis.

Scenario screen. The 128-point Latin-hypercube ranges are not probability distributions, so the result is not calibrated uncertainty quantification. Spearman rank correlations: conversion with CO2 excess +0.59, reference reaction rate +0.57 and heat-transfer strength -0.51; peak temperature with feed temperature +0.53, reference reaction rate +0.46, heat-transfer strength -0.43 and reaction-heat magnitude +0.38.

Not established: validated kinetics for a commercial catalyst, plant conversion or selectivity, dynamic stability or relief sizing, gas-liquid mass transfer, phase equilibrium and pressure drop, and parameter uncertainty from experimental data.
