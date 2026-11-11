# Thermodynamics of Solution

A C++17 experimental aqueous chemical-equilibrium calculator, with modular ideal and nonideal activity models, a variable-temperature equilibrium engine, automated tests, and command-line inputs/outputs.

**Status:** MVP+1 implementation for *restricted aqueous chemistry*, **not** validated for engineering design, nuclear safety, regulatory use, or experimental prediction without reference data and further validation.

## Features

- Arbitrary declared aqueous basis components and secondary complexes with **signed integer formation stoichiometries**.
- Fixed aqueous component totals (mol/L), mass-action constraints, and optional electroneutrality.
- The balance method adjusts exactly **one** basis component (e.g., H+) as `free`, replacing its specified total with charge balance. If no free component exists, charge imbalance is only reported.
- Ideal activities or an approximate temperature-dependent Davies activity model for dilute mixtures.
- Model temperature from **273.15 to 323.15 K (0–50 °C)** and constant-enthalpy van't Hoff corrections from a specified reference temperature.
- pH calculated as `-log10(a_H+)` when `H+` exists in the model.
- Damped Newton iteration with finite-difference Jacobian and pivoted linear solves.
- Parseable JSON and CSV outputs, residual diagnostics, deterministic command-line behavior, and GitHub Actions CI.

## Build and run

Requires CMake 3.16+ and a C++17 compiler.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
ctest --test-dir build --output-on-failure

# New MVP+1 aqueous solver:
./build/aqueous_cli examples/carbonate_ideal.txt
./build/aqueous_cli examples/carbonate_davies.txt results.json species.csv

# Phase 2 ideal association solver remains available:
./build/thermo ideal examples/association.txt
```

Files are relative to the current working directory. `aqueous_cli` returns nonzero on a parse, convergence or model validation error.

## Input model schema

One record per line; whitespace-separated; lines beginning `#` are comments.

```text
temperature_k 298.15
reference_temperature_k 298.15
activity_model davies

# component NAME INTEGER_CHARGE TOTAL_MOL_PER_L|free
component H+ 1 free
component CO2* 0 0.001
component Na+ 1 0.001

# complex NAME INTEGER_CHARGE LOG10_K_REF DELTA_H_KJ_MOL
#         NUMBER_OF_BASIS_TERMS [BASIS_NAME SIGNED_INTEGER_STOICH]...
complex OH- -1 -14.0 55.8 1 H+ -1
complex HCO3- -1 -6.35 7.7 2 CO2* 1 H+ -1
complex CO3-- -2 -16.68 22.5 2 CO2* 1 H+ -2
```

For the `OH-` record, water is treated as a solvent of constant activity, and `K_w` is an effective equilibrium constant. This is **not** a full explicit solvent mass balance.

The example `CO2*` species represents lumped dissolved carbon dioxide and carbonic acid under a chosen convention. The reaction constants and enthalpies above are **illustrative and unverified**, rather than an internally validated thermodynamic database.

Each complex charge must equal the charge-weighted sum of its basis stoichiometry. All component totals must be strictly positive. A `free` component is determined by electroneutrality, rather than an independent total. Up to one free component is supported.

## Mathematical formulation

For each complex `j`, using standard concentration `c° = 1 mol/L`:

```text
a_j = 10^(log10 K_j(T)) * product_i (a_i ^ nu_ij)
a_s = gamma_s * (c_s/c°)

T_i = sum_s (nu_is * c_s), for each fixed-total component i
0 = sum_s z_s*c_s, when one component is marked free
I = 0.5 * sum_s (z_s^2 * c_s) [approximate molarity basis]
pH = -log10(a_H+)
```

Temperature adjustments use constant formation enthalpy:

```text
log10 K(T) = log10 K(T_ref)
    - DeltaH/(R ln(10)) * (1/T - 1/T_ref)
```

The nonideal Davies approximation uses:

```text
log10 gamma_i = -A(T) * z_i^2 * (sqrt(I)/(1+sqrt(I)) - 0.3 I)
```

Here A(T) is estimated from empirical pure-water dielectric and density approximations. The code uses mol/L as an approximation to the conventional molal ionic strength; this makes the nonideal activity coefficients **approximate**, and it is not a rigorous molal activity model. Neutral-species activity coefficients are taken as unity. The Davies option rejects solutions above **0.1 mol/L ionic strength** to keep the scope conservative; this is a software-domain restriction, not a universal statement of Davies model validity.

The nonlinear solve uses base-10 logarithmic free concentrations and, for Davies, an additional log ionic-strength unknown. Residuals are scaled by component totals and a charge scale; finite-difference Jacobian, partial-pivot elimination, and a backtracking line search seek a solution. Failure is signaled, not masked.

## Outputs

`aqueous_cli MODEL [JSON_FILE] [CSV_FILE]` prints JSON or saves JSON and CSV. JSON fields include:

- `temperature_K`, `activity_model`, `iterations`, `normalized_residual`
- `ionic_strength_mol_L`, `charge_residual_mol_L` and optional `pH`
- `concentrations_mol_L`, `activities`, `component_residuals_mol_L`

Only numerical success and mass-action/conservation consistency are established by the current tests. These are **verification** checks, not independent chemical accuracy validation.

## Automated test coverage

- Synthetic analytically solvable ideal association with benchmark free concentrations.
- Ideal and Davies numerical behavior; expected direction of charged activity coefficients.
- Carbonate species mass conservation, water mass action and charge-balance residual.
- pH sensitivity to nonideal activities and temperature.
- Formation-constant temperature dependence.
- Rejection of inconsistent reaction charges, unphysical temperature, and excessive ionic strength.

## Scientific limitations and outstanding validation

1. The included carbonate formation constants and enthalpies are not proven against a single curated literature database; **temperature-dependent pH predictions are not experimentally validated**.
2. No robust benchmark comparison with PHREEQC or other independent reference code has been completed. That comparison is the next release gate.
3. No explicit solvent chemical potential, ionic-size specific models, ion pairing parameterization, pressure dependence, redox, or phase change.
4. The legacy Pitzer-inspired functions are incomplete and **are not coupled** to the aqueous solver.
5. The approximate Davies concentration convention should be replaced by consistent molality conversion and temperature-dependent water activity before higher-accuracy work.
6. Initial guesses and numerical conditioning may fail for extreme or strongly coupled systems. The numerical implementation is a specialized solver, not a production geochemical engine.

## Architecture

- `include/aqueous_engine.h` — model data, results, public interface
- `src/aqueous_engine.cpp` — model validation, thermodynamics and nonlinear solver
- `src/aqueous_cli.cpp` — command-line interface and serialization
- `examples/carbonate_*.txt` — illustrative dilute carbonate/water input files
- `tests/aqueous_tests.cpp` — numerical verification and safeguards
- `speciation.cpp` — earlier ideal-only Phase 2 solver, retained for compatibility

The high-priority next steps are: independent PHREEQC comparison using identical constants, more chemical test systems, automated validation reports, a consistent thermodynamic data provenance scheme, and stricter dimensional/standard-state treatment.
