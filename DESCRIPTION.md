# Pryngles — Package Description (Engineering & Science)

This document provides a **global, engineering- and science-oriented description** of the `pryngles` codebase for a team that will maintain, extend, and modernize it.

`pryngles` (PlanetaRY spaNGLES) is a Python package for **forward modeling the photometric and polarimetric signatures** of exoplanetary systems, with an emphasis on **ringed exoplanets** (exorings) in both **transiting** and **non-transiting** configurations.

The package is built around a discretization strategy: approximate extended surfaces (planetary spheres, rings/disks, atmospheres) using many small surface facets (“**spangles**”, historically “sequins”), compute illumination/visibility for each facet, and integrate the resulting contributions to generate time series (“light curves”) and Stokes observables.

<p align="center">
  <img src="https://raw.githubusercontent.com/seap-udea/pryngles/kiss/gallery/description-spangles.png" alt="Spangles (planet + ring discretization)" width="70%"/>
</p>

---

## Scientific scope

`pryngles` targets scenarios where the observed signal is shaped by:

- **Reflected light**: stellar light scattered by planetary/ring surfaces (including realistic phase functions through configurable “scatterers”).
- **Transit / occultation**: flux removed by the projected overlap of an extended object (planet + rings) with the stellar disk, including limb darkening.
- **Thermal emission**: planetary radiation integrated over a bandpass, driven by parametric temperature maps (including tidally locked prescriptions).
- **Polarization**: polarized reflected light expressed through **Stokes** parameters and derived quantities (e.g., degree of polarization).

The scientific motivation and validation of the approach is documented in the peer-reviewed literature already referenced in the repository:

- Veenstra, Zuluaga, Alvarado-Montes, Sucerquia & Stam (2025), [*A general polarimetric model for transiting and nontransiting ringed exoplanets* (PDF)](https://github.com/seap-udea/pryngles/blob/kiss/doc/papers/pdfs/2025-veenstra-zuluaga-alvarado-montes-sucerquia-stam_AA_693_A310.pdf), A&A 693, A310, doi:10.1051/0004-6361/202347194  
- Zuluaga, Sucerquia & Alvarado-Montes (2022), [*The bright side of the light curve: a general photometric model for non-transiting exorings* (PDF)](https://github.com/seap-udea/pryngles/blob/kiss/doc/papers/pdfs/2022-zuluaga-sucerquia-alvarado-montes_AsCom_40_100623_arXiv2207.08636.pdf), AsCom 40, 100623, doi:10.1016/j.ascom.2022.100623  
- Sucerquia, Alvarado-Montes, Zuluaga, Montesinos & Bayo (2020), [*Scattered light may reveal the existence of ringed exoplanets* (PDF)](https://github.com/seap-udea/pryngles/blob/kiss/doc/papers/pdfs/2020-sucerquia-alvarado-montes-zuluaga-montesinos-bayo_MNRASL_496_L85_arXiv2004.14121.pdf), MNRASL 496(1), L85-L90, doi:10.1093/mnrasl/slaa080

---

## Core modeling idea: spangles (surface discretization)

The package discretizes extended bodies into **spangles** that carry the minimal state needed to compute observables:

- **Geometry**: positions and normals in relevant frames, projected areas, limb angles, azimuths.
- **Visibility/illumination flags**: which spangles are visible to the observer and illuminated by the star(s), including mutual shadowing/occultation when modeled.
- **Optical properties**: assigned “scatterer” objects (surface/atmosphere/ring) used to evaluate albedo/phase functions and polarization behavior.

The figure below shows an example of spangle state identification for a specific configuration. The color of each spangle (point) is assigned based on whether it receives stellar light, is visible to the observer, is shadowed, and related geometry/illumination conditions.

<p align="center">
  <img src="https://raw.githubusercontent.com/seap-udea/pryngles/kiss/gallery/description-spangles2.png" alt="Legacy visualization: ringed planet spangles (RingedPlanet interface)" width="70%"/>
</p>

The figure below illustrates the geometric quantities involved in the flux computation for each spangle.

<p align="center">
  <img src="https://raw.githubusercontent.com/seap-udea/pryngles/kiss/gallery/spangles_directions.png" alt="Spangle geometric directions" width="70%"/>
</p>

Spangles are stored and updated in a **vectorized tabular representation** (pandas `DataFrame`) to support large \(N\) (thousands to millions) of facets with efficient numerical updates.

---

## Interfaces: legacy `RingedPlanet` vs modern `System`

Because `pryngles` evolved in the context of research, it currently exposes two high-level interfaces.

### `System` interface (recommended)

`System` is the modern, general interface for **arbitrary hierarchical systems**:

- Stars (including multiple light sources), planets, rings/disks, moons, observers.
- A unified computation pipeline that integrates:
  - orbital dynamics (via `rebound`),
  - spangling (surface discretization),
  - per-time-step perspective integration,
  - effect synthesis (transit/reflection/emission/polarization).

The canonical workflow is:

1. Build the body tree (`System.add(...)`)
2. Initialize the dynamical simulation (`System.initialize_simulation()`)
3. Generate spangles (`System.spangle_system()`)
4. Compute time series (`System.compute_lightcurve(...)`)

The most complete demonstration of this interface is the tutorial notebook:

- [tutorials/System-Interface-Tutorial.ipynb](https://github.com/seap-udea/pryngles/blob/kiss/tutorials/System-Interface-Tutorial.ipynb)

### `RingedPlanet` interface (legacy)

`RingedPlanet` is the original interface used in early papers and notebooks. It encodes a more specialized geometry and historical assumptions. It remains in the repository for backward compatibility and reproducibility of prior results, but **long-term development is expected to converge on the `System` interface**.

In the current codebase, the two interfaces can be partially bridged (e.g., constructing a `RingedPlanet` representation from a `System` description via `System.ensamble_system()`).

Complete interface examples are available in:

- [tutorials/Quickstart.ipynb](https://github.com/seap-udea/pryngles/blob/kiss/tutorials/Quickstart.ipynb)

---

## Package structure (high-level)

The modern code lives under `src/pryngles/`. The most relevant modules for engineering work are:

- `system.py`
  - **`System`**: orchestrates the entire modeling pipeline (body registry, simulation, spangling, per-time integration, effect synthesis, and lightcurve outputs).
- `body.py`
  - **`Body`** and specialized body kinds (Star/Planet/Ring/Observer), plus helpers such as **`Detector`** used to generate synthetic observational signals (noise/uncertainties).
- `spangler.py`
  - **`Spangler`**: generates and manages the spangle discretization, storing spangle state and providing geometry/plotting utilities.
- `scatterer.py`
  - **`Scatterer`** (abstract interface) and built-in scatterers (e.g., gray Lambertian surfaces/atmospheres, blackbody-like behavior). Scatterers encapsulate reflectance laws used by reflection/polarization computations.
- `orbit.py`
  - Orbital hierarchy primitives and integration utilities (wrapping `rebound`).
- `science.py`
  - Numerical and geometric utilities (coordinate transforms, angular relations, helper integrals), used throughout the modeling pipeline.
- `legacy.py`
  - Legacy `RingedPlanet` implementation retained for compatibility.
- `consts.py`
  - Physical constants, units, and package-wide enumerations/defaults.
- `plot.py`
  - Higher-level plotting and visualization utilities.

The documentation is built with Sphinx under `docs/` and is published on Read the Docs: https://pryngles.readthedocs.io

The figure below summarizes the package architecture.

<p align="center">
  <img src="https://raw.githubusercontent.com/seap-udea/pryngles/kiss/gallery/pryngles-architecture.png" alt="Pryngles architecture" width="100%"/>
</p>

---

## Data flow and computation pipeline (engineering view)

At a high level, `pryngles` turns a system definition into observables by repeatedly applying four stages:

1. **System definition**
   - A hierarchical tree of `Body` objects is created (primary/child relationships, orbital elements, optical/physical properties).

2. **Dynamics and geometry**
   - A `rebound.Simulation` evolves the orbital state.
   - For each time \(t\), the system computes the **observer perspective** (line-of-sight unit vector, projected areas, azimuths) and the **illumination state** (star direction(s), incidence angles).

3. **Effect updates (per time step)**
   - Transit: compute fractional stellar flux removed by occulting spangles.
   - Reflection: compute incident flux per spangle and apply the scatterer-defined albedo/phase behavior; integrate visible contributions.
   - Thermal emission: compute band-integrated emission from spangles given temperature prescriptions and viewing geometry.
   - Polarization: compute Stokes contributions from scattering geometry; integrate to total Stokes and derived degree of polarization.

4. **Aggregation and products**
   - `System.compute_lightcurve(...)` produces a time-indexed output table (multi-indexed by body and effect) and can optionally generate a detector-level synthetic signal (photon counting / noise model).

This design makes it possible to add new physics by introducing:

- a new per-spangle state field,
- a new update method in `System` (vectorized over spangles),
- and, when needed, a new `Scatterer` implementation for reflectance/polarization laws.

---

## Capabilities (current)

Based on the current tutorials and code, `pryngles` supports:

- **Photometric light curves** for star + planet + rings configurations, including:
  - transit/occultation and limb darkening,
  - diffuse reflection from surfaces and rings,
  - thermal emission over a wavelength band.
- **Polarimetric modeling**:
  - Stokes flux integration and degree of polarization time series for reflected light.
- **Multi-body systems**:
  - multiple planets, moons, rings/disks, and (in modern structure) multiple light sources.
- **Scientific visualization**:
  - geometric previews of the system and spangles,
  - plotting utilities and notebook-driven workflows.
- **Synthetic observation layer**:
  - optional detector simulation producing noisy signals and uncertainties.

The figure below shows an example simulated light curve for WASP-43b, including primary/secondary transit geometry, thermal emission, and polarization signal.

<p align="center">
  <img src="https://raw.githubusercontent.com/seap-udea/pryngles/kiss/gallery/wasp43_synthetic_lightcurve_combined.png" alt="Synthetic combined light curve for WASP-43b" width="70%"/>
</p>

<p align="center">
  <img src="https://raw.githubusercontent.com/seap-udea/pryngles/kiss/gallery/wasp43_polarization_lightcurve.png" alt="Polarization light curve for WASP-43b" width="70%"/>
</p>

---

## Practical usage logic (what users do)

In user-facing workflows (especially the `System` interface), usage typically follows:

### 0) Define constants and units

Before building the system, define physical parameters and convert them to the internal unit system used by `System`.

```python
# Orbital parameters
period    = 0.81347753           # Orbital period [days]
a_abs     = 0.01526              # Semi-major axis [AU]
a         = 4.857                # Semi-major axis [R_star]
inc       = 90                   # Inclination [deg]
ecc       = 0.0                  # Eccentricity
omega     = 0.0                  # Argument of periastron [deg]

# Stellar parameters
T_star    = 4520                 # Effective temperature [K]
M_star    = 0.717  * pr.Consts.msun
R_star_m  = (a_abs / a) * pr.Consts.au  # Stellar radius [m]

# Planet parameters
M_planet  = 0.0   * pr.Consts.mjupiter
R_planet  = 1.036 * pr.Consts.rjupiter

# Ring geometry
R_in      = 1.5                  # Inner ring radius [R_planet]
R_out     = 2.5                  # Outer ring radius [R_planet]
ring_inc  = 60                   # Ring inclination [deg]
ring_tau  = 0.4                  # Ring optical depth

# Observing band
lambda_min = 1.1e-6              # [m]
lambda_max = 1.7e-6              # [m]
```


### 1) Create the system

Define bodies (star, planet, ring) with geometry, orbital parameters, and optical models.

Minimal `System` setup (from `tutorials/System-Interface-Tutorial.ipynb`):

```python
# Create the system
system = pr.System()

# Add the central star
star = system.add(
    kind="Star",
    m=M_star / system.um,              # Mass in system units
    T_eff=T_star,                      # Effective temperature [K]
    radius=R_star_m / system.ul,       # Radius in system units
    limb_coeffs=[0.65],                # Limb-darkening coefficients
)
```

Add the planet and ring:

```python
planet = system.add(
  kind="Planet",
  parent=star,
  m=M_planet / system.um,
  radius=R_planet / system.ul,
  a=a_abs,
  inc=pr.DEG * (90 - inc),
  e=ecc,
  omega=pr.DEG * omega,
)

ring = system.add(
  kind='Ring',
  parent=planet,
  fi=R_in,
  fe=R_out,
  i=pr.DEG * ring_inc,
  taur=ring_tau,
)
```

Unit conversion rules used throughout the workflow:

- Divide SI masses by `system.um` when passing `m=...`.
- Divide SI lengths by `system.ul` when passing `radius=...`.
- Convert time arrays from days to internal units with `times_days * pr.Consts.day / system.ut`.

### 2) Initialize simulation and spangle surfaces

Initialize `rebound`, build the spangle discretization, and compute any precomputed geometry needed for fast time stepping.

```python
# Initialize the dynamical simulation with Rebound
system.initialize_simulation()

# Discretize the system (spangle)
system.spangle_system()
```

You can inspect the system configuration by plotting a specific orbital perspective:

```python
system.n_obs = spy.eul2m(
    np.deg2rad(omega),
    np.deg2rad(inc),
    0,
    3, 1, 3,
)[0]

system.integrate_perspective(0)
system.sg.plot2d()
```

<p align="center">
  <img src="https://raw.githubusercontent.com/seap-udea/pryngles/kiss/gallery/wasp43_system_view.png" alt="System geometry view (primary transit)" width="70%"/>
</p>

### 3) Define observer and compute light curve

Run `compute_lightcurve(times=..., effects=[...])` for the desired observables.
Define detector properties here so they can be reused when computing transit and when generating synthetic measurements later.

Example computing transit, emission, and polarization (excerpted from the same tutorial):

```python
DETECTOR_PROPERTIES = {
  'wavelength_min': lambda_min,
  'wavelength_max': lambda_max,
  'apperture': 0.5,
  'quantum_eff': 0.9,
  't_cadence': 15 * 60,
  'distance': 100 * pr.Consts.pc,
}

# Transit
system.compute_lightcurve(
    times=times_system,
    effects=['transit'],
    signal=DETECTOR_PROPERTIES,
)
transit_flux = system.lightcurve['transit']

# Emission
system.compute_lightcurve(
    times=times_system,
    bandwidth=(lambda_min, lambda_max),
    effects=['emission'],
)
emission_flux = system.lightcurve['emission']

# Polarization
system.compute_lightcurve(
    times=times_system,
    effects=['polarization'],
)
polarized_flux = system.lightcurve['scattering']
polarized_degree = system.lightcurve['polarization']
```

### 4) Analyze / visualize

Plot fluxes and polarization, compare planet-only vs planet+ring contributions, and iterate on parameters.

```python
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# Polarization degree
ax1.plot(times_days, polarized_degree.sum(axis=1), color='mediumvioletred', linewidth=2, label='Total')
ax1.plot(times_days, polarized_degree['Planet'], linewidth=2, label='Planet')
ax1.plot(times_days, polarized_degree['Ring'], linewidth=2, label='Ring')
ax1.set_ylabel('Degree of polarization [-]')
ax1.legend(loc='upper right')
ax1.grid(True, alpha=0.3, linestyle=':')

# Flux components
ax2.plot(times_days, total_flux, color='C0', linewidth=2, label='Light curve')
ax2.plot(times_days, 1 + transit_flux.sum(axis=1), color='C1', linewidth=2, label='Transit')
ax2.plot(times_days, 1 + emission_flux.sum(axis=1), color='C2', linewidth=2, label='Emission')
ax2.set_xlabel('Phase [days]')
ax2.set_ylabel('Normalized flux [-]')
ax2.legend(loc='best')
ax2.grid(True, alpha=0.3, linestyle=':')

plt.tight_layout()
plt.show()
```

<p align="center">
  <img src="https://raw.githubusercontent.com/seap-udea/pryngles/kiss/gallery/wasp43_polarization_lightcurve.png" alt="Transit light curve for WASP-43b" width="70%"/>
</p>

### 5) Generate an instrument-specific signal

Use the detector layer to generate a synthetic observed signal with realistic uncertainties. This assumes the transit light curve was already computed in step 3 with `signal=DETECTOR_PROPERTIES`.

```python
signal_times, signal_flux, signal_error = system.detector.generate_signal(
  times_system * system.ut,
  total_flux,
)
```

Typical visualization of model + synthetic measurements:

```python
fig, ax = plt.subplots(figsize=(10, 5))
ax.errorbar(
  signal_times / pr.Consts.days,
  signal_flux,
  yerr=signal_error,
  markersize=2,
  fmt='o',
  color='black',
  alpha=1,
  label='Simulated signal',
)
ax.plot(times_days, total_flux, color='C0', linewidth=2, label='Light curve')
ax.set_xlabel('Phase [days]')
ax.set_ylabel('Normalized flux [-]')
ax.legend(loc='best')
ax.grid(True, alpha=0.3, linestyle=':')
plt.tight_layout()
plt.show()
```

<p align="center">
  <img src="https://raw.githubusercontent.com/seap-udea/pryngles/kiss/gallery/wasp43_synthetic_lightcurve_full.png" alt="Full synthetic light curve" width="70%"/>
</p>

---

For worked examples, see:

- [tutorials/Quickstart.ipynb](https://github.com/seap-udea/pryngles/blob/kiss/tutorials/Quickstart.ipynb) (minimal entry point aligned with `README.md`)
- [tutorials/System-Interface-Tutorial.ipynb](https://github.com/seap-udea/pryngles/blob/kiss/tutorials/System-Interface-Tutorial.ipynb) (complete workflow and combined effects)

---

## Engineering constraints and expectations

### Runtime and performance

The primary computational cost is the **number of spangles** times the **number of time samples** times the **number of enabled effects**.

Engineering improvements typically focus on:

- vectorization over spangles (pandas/numpy),
- minimizing per-time Python overhead,
- caching geometry that is time-invariant,
- reducing redundant trigonometric work across effects.

### Numerical robustness

The package heavily relies on third-party scientific libraries. Recent changes in dependencies (e.g., pandas strict dtype behavior, SciPy deprecations) can materially affect execution. Modernization work should include:

- explicit dtype management for spangle tables,
- careful handling of array shapes for interpolators and quadrature,
- repeatable notebook-based regression checks for representative scenarios.

### Python version

The repository targets **Python 3.12+** (see `README.md` and packaging configuration).

---

## Links and references inside the repo

- **README**: [README.md](https://github.com/seap-udea/pryngles/blob/kiss/README.md) (badges, papers, installation, quickstart, citations)
- **Release notes**: [WHATSNEW.md](https://github.com/seap-udea/pryngles/blob/kiss/WHATSNEW.md)
- **Tutorials**: [Quickstart.ipynb](https://github.com/seap-udea/pryngles/blob/kiss/tutorials/Quickstart.ipynb), [System-Interface-Tutorial.ipynb](https://github.com/seap-udea/pryngles/blob/kiss/tutorials/System-Interface-Tutorial.ipynb)
- **Core package**: [src/pryngles/](https://github.com/seap-udea/pryngles/tree/kiss/src/pryngles)

## External resources

- **GitHub repository**: https://github.com/seap-udea/pryngles
- **PyPI package**: https://pypi.org/project/pryngles/
- **API documentation (Read the Docs)**: https://pryngles.readthedocs.io
- **ASCL entry**: https://ascl.net/2205.016

