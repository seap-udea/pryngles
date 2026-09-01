"""Basic smoke tests for the pryngles package and core utilities."""

import pytest

import pryngles as pr


def test_package_version():
    assert pr.version == "1.0.2"


def test_consts_basic_constants():
    assert pr.Consts.ppm == 1e6
    assert pr.Consts.ppb == 1e9
    assert pr.Consts.rad * pr.Consts.deg == pytest.approx(1.0)
    # Physical constants from Rebound (SI)
    assert pr.Consts.au > 0
    assert pr.Consts.msun > 0


def test_flatten_treats_strings_as_scalars():
    """Misc.flatten must not iterate over strings (regression for Python-2
    `basestring` NameError)."""
    assert list(pr.Misc.flatten([["a", "b", ["c"]]])) == ["a", "b", "c"]
    assert list(pr.Misc.flatten(["ab"])) == ["ab"]
    assert list(pr.Misc.flatten([1, [2, [3]]])) == [1, 2, 3]


def test_sampler_circle_geometry():
    sp = pr.Sampler(N=200, seed=1)
    sp.gen_circle()
    assert sp.dim == 2
    assert sp.ss.shape == (200, 3)
    assert sp.A == pytest.approx(3.141592653589793)


def test_sampler_sphere_geometry():
    sp = pr.Sampler(N=200, seed=1)
    sp.gen_sphere()
    assert sp.dim == 3
    assert sp.ss.shape == (200, 3)
    # Unit sphere: every sample point should be ~1 from origin.
    norms = (sp.ss**2).sum(axis=1) ** 0.5
    assert norms == pytest.approx(1.0, abs=1e-6)
