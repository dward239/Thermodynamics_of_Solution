# Thermodynamics of Solution

Exploratory C++ numerical components for chemical solution calculations. This repository is **a work in progress**, not a validated thermodynamic equilibrium package.

## Current scope

- Load aqueous species and solid-phase records from a text file.
- Load binary interaction parameters.
- Calculate ionic strength from species concentrations and charges.
- Estimate activity coefficients using an **illustrative, incomplete** Pitzer-inspired expression.
- Calculate a saturation index from **supplied** activities and log10(Ksp).

**Not yet implemented:** A convergent, validated chemical speciation or multiphase equilibrium solver, electroneutrality and mass-balance coupling, temperature-dependent parameters, and a full published Pitzer formalism.

The original prototype contained duplicate headers, undeclared solver routines, incomplete input data, and a CI configuration referencing nonexistent commands. This refactor establishes a compilable, testable foundation without presenting those incomplete calculations as scientifically validated.

## Build

Requires CMake >= 3.16 and a C++17 compiler. No Eigen dependency is currently required.

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
./build/thermo
```

To load a species database:

```sh
./build/thermo path/to/species.txt path/to/pitzer_params.txt
```

The checked-in `database.txt` is empty and is not a usable chemistry dataset. Records in `pitzer_params.txt` and `solids.txt` are examples; their scientific provenance and validity have not been established. Do not use them for design, safety, or regulatory decisions.

## Next engineering steps

1. Define a consistent aqueous reaction/species database schema (primary species, secondary complexes, reactions, log K, units, standard states).
2. Supply traceable literature data and validation cases.
3. Formulate independent elemental balances, electroneutrality, and mass-action residuals.
4. Add a safeguarded nonlinear solver with meaningful convergence reporting.
5. Implement and validate a complete activity-coefficient model before coupling it to the solver.
