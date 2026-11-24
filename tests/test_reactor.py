import numpy as np
import pytest

from reactor import (
    ReactorParams,
    full_conversion_adiabatic_rise,
    simulate,
    thermal_ntu,
    validate_params,
)
from operating_window import operating_grid, pareto_mask
from thermal_sensitivity import thermal_amplification


def test_stoichiometry_is_exact_along_reactor():
    p = ReactorParams()
    res = simulate(p)
    assert np.allclose(p.F_EO0 - res["F_EO"], res["F_EC"], atol=1e-8)
    assert np.allclose(p.F_CO20 - res["F_CO2"], res["F_EC"], atol=1e-8)


def test_conversion_is_monotone_bounded_and_flows_nonnegative():
    res = simulate(ReactorParams())
    assert (np.diff(res["X"]) >= -1e-10).all()
    assert 0.0 <= res["X"][-1] <= 1.0
    assert np.min(res["F_EO"]) >= -1e-10
    assert np.min(res["F_CO2"]) >= -1e-10


def test_limiting_reactant_depletion_is_handled_without_negative_flow():
    p = ReactorParams(ratio_CO2=0.5, A=2.0e6, T0=400.0)
    res = simulate(p)
    assert res["depleted"] is True
    assert 0.0 < res["depletion_volume_m3"] < p.V_total
    assert res["X"][-1] == pytest.approx(0.5, abs=1e-10)
    assert np.min(res["F_CO2"]) >= -1e-12
    # After reaction stops, temperature relaxes toward the coolant.
    dep_idx = np.searchsorted(res["V"], res["depletion_volume_m3"])
    assert res["T"][-1] < np.max(res["T"][max(dep_idx - 1, 0):])


def test_energy_balance_closes_numerically():
    res = simulate(ReactorParams(), n_points=600)
    assert np.max(np.abs(res["energy_residual_rel"])) < 1e-5


def test_reference_scales_are_reproducible():
    p = ReactorParams()
    assert full_conversion_adiabatic_rise(p) == pytest.approx(79.1666667, rel=1e-8)
    assert thermal_ntu(p) == pytest.approx(4.0, rel=1e-12)


def test_hotter_feed_raises_peak_temperature():
    cold = np.max(simulate(ReactorParams(T0=370.0))["T"])
    hot = np.max(simulate(ReactorParams(T0=390.0))["T"])
    assert hot > cold


def test_warmer_coolant_raises_absolute_peak_at_fixed_feed():
    cold = np.max(simulate(ReactorParams(T0=400.0, T_cool=365.0))["T"])
    warm = np.max(simulate(ReactorParams(T0=400.0, T_cool=380.0))["T"])
    assert warm > cold


def test_thermal_state_and_sensitivity_are_not_same_diagnostic():
    T0s = np.linspace(370.0, 400.0, 13)
    peaks_cold, amp_cold = thermal_amplification(T0s, 365.0)
    peaks_warm, amp_warm = thermal_amplification(T0s, 380.0)
    assert peaks_warm[-1] > peaks_cold[-1]  # warmer coolant => hotter absolute state
    assert np.max(amp_cold) > np.max(amp_warm)  # but stronger local amplification under cold coolant
    assert np.max(amp_cold) > 1.0


def test_pareto_mask_contains_only_nondominated_points():
    conversion = np.array([[0.80, 0.90], [0.85, 0.92]])
    peak = np.array([[400.0, 430.0], [390.0, 450.0]])
    mask = pareto_mask(conversion, peak)
    c = conversion.ravel()
    t = peak.ravel()
    for i, keep in enumerate(mask.ravel()):
        if not keep:
            continue
        dominated = (c >= c[i]) & (t <= t[i]) & ((c > c[i]) | (t < t[i]))
        dominated[i] = False
        assert not dominated.any()


def test_operating_grid_conversion_and_peak_are_finite():
    T0s, ratios, conversion, peak = operating_grid()
    assert conversion.shape == (len(T0s), len(ratios))
    assert np.isfinite(conversion).all()
    assert np.isfinite(peak).all()
    assert ((conversion >= 0) & (conversion <= 1)).all()


@pytest.mark.parametrize(
    "kwargs",
    [
        {"ratio_CO2": 0.0},
        {"v_flow": -1.0},
        {"Ua": 0.0},
        {"dH": 95_000.0},
    ],
)
def test_invalid_parameters_fail_loudly(kwargs):
    with pytest.raises(ValueError):
        validate_params(ReactorParams(**kwargs))


def test_n_points_must_be_at_least_two():
    with pytest.raises(ValueError):
        simulate(ReactorParams(), n_points=1)
