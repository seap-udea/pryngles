"""Focused tests for the modern ``System`` interface and related regressions.

These cover the full pipeline (``initialize_simulation`` ->
``spangle_system`` -> ``compute_lightcurve``) for each effect, plus
regression tests for specific bugs found in the codebase.
"""

import numpy as np
import pandas as pd
import pytest

import pryngles as pr
from pryngles import Consts, Scatterer, SCATTERERS_CATALOGUE


def build_system(nspangles=200, ring=True):
    """Return a ready (initialized + spangled) Star/Planet[+/Ring] System."""
    system = pr.System()
    star = system.add(
        kind="Star", radius=Consts.rsun / system.ul, limb_coeffs=[0.65]
    )
    planet = system.add(
        kind="Planet",
        parent=star,
        a=0.2,
        e=0.0,
        radius=Consts.rsaturn / system.ul,
        nspangles=nspangles,
    )
    if ring:
        system.add(
            kind="Ring",
            parent=planet,
            fi=1.5,
            fe=2.5,
            i=30 * Consts.deg,
            nspangles=nspangles,
        )
    system.initialize_simulation()
    system.spangle_system()
    return system


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------


def test_initialize_and_spangle():
    system = build_system()
    assert system._simulated is True
    assert system._spangled is True
    assert isinstance(system.data, pd.DataFrame)
    assert len(system.data) > 0
    assert set(["Star", "Planet", "Ring"]).issubset(set(system.data.name.unique()))


def test_orbit_period_available():
    system = build_system()
    assert system.orbit.Ps[0] > 0


def test_compute_lightcurve_transit():
    system = build_system()
    times = np.array([0.0, 0.5]) * system.orbit.Ps[0]
    lc = system.compute_lightcurve(times=times, effects=["transit"])
    assert lc["total_flux"].shape == (2,)
    assert np.all(np.isfinite(lc["total_flux"]))
    # Transit removes light -> the transit contribution is <= 0.
    assert (lc["transit"].sum(axis=1) <= 0).all()


def test_compute_lightcurve_reflection():
    system = build_system()
    times = np.array([0.0, 0.5]) * system.orbit.Ps[0]
    lc = system.compute_lightcurve(times=times, effects=["reflection"])
    assert lc["total_flux"].shape == (2,)
    assert np.all(np.isfinite(lc["total_flux"]))
    # Reflected light is added on top of the stellar baseline (>= 0).
    assert (lc["reflection"].sum(axis=1) >= 0).all()


def test_compute_lightcurve_emission():
    system = build_system()
    planet = system.bodies["Planet"]
    planet.set_temperature_model(
        {"type": "Two Temperature", "params": {"T_day": 1800, "T_night": 400}}
    )
    times = np.array([0.0]) * system.orbit.Ps[0]
    lc = system.compute_lightcurve(times=times, effects=["emission"])
    assert lc["total_flux"].shape == (1,)
    assert np.isfinite(lc["total_flux"][0])
    assert (lc["emission"] >= 0).all().all()


def test_compute_lightcurve_unknown_effect_rejected():
    system = build_system()
    with pytest.raises(ValueError):
        system.compute_lightcurve(times=np.array([0.0]), effects=["bogus"])


@pytest.mark.slow
def test_compute_lightcurve_polarization():
    system = build_system(nspangles=120)
    times = np.array([0.0]) * system.orbit.Ps[0]
    lc = system.compute_lightcurve(times=times, effects=["polarization"])
    assert "scattering" in lc
    assert "polarization" in lc
    # Scattered Stokes I must be finite for every body.
    assert np.all(np.isfinite(lc["scattering"].values))


# ---------------------------------------------------------------------------
# Regression tests for specific bugs
# ---------------------------------------------------------------------------


def test_polarization_ringless_planet_raises():
    """Polarization currently requires a Ring child on each Planet; a ring-less
    planet must fail with a clear error rather than an IndexError."""
    system = build_system(ring=False)
    with pytest.raises(ValueError, match="requires a Ring"):
        system.update_Polarization()


def test_read_fourier_data_creates_scatterers():
    """_read_Fourier_data used trailing-comma tuples, which crashed when
    building StokesScatterer instances."""
    s = pr.System()
    s._read_Fourier_data()
    assert s.nmatp > 0
    assert s.nmatr > 0
    assert hasattr(s, "SCp") and hasattr(s, "SCr")


def test_scatterer_reset_catalogue():
    """Scatterer.reset_catalogue used to only rebind a local name."""
    pr.NeutralSurface()  # registers a new scatterer
    assert len(SCATTERERS_CATALOGUE) > 0
    Scatterer.reset_catalogue()
    assert len(SCATTERERS_CATALOGUE) == 0
