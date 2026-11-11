# Thermodynamics of Solution

Experimental C++ software for exploring solution thermodynamics and chemical equilibria.

## Phase 2: ideal aqueous association solver

A small **ideal-solution mass-action solver** is now provided. It is deliberately limited to finite, positive component totals and aqueous complexes expressed as positive-integer combinations of independent primary components. It enforces each primary-component material balance through damped Newton iterations in logarithmic concentration coordinates.

For each complex j:

```text
c_j / c° = 10^(log10 K_j) * product_i (c_i / c°)^(nu_ij)
T_i = c_i + sum_j nu_ij * c_j
c° = 1 mol/L
```

Here K is an assumed, dimensionless formation constant at fixed temperature. With ideal activities, activity coefficients are unity. Formation constants supplied in the example files are **synthetic**, not chemical literature data.

The solver uses a scaled material-balance residual, analytic Jacobian, pivoted Gaussian elimination, and backtracking line search. Nonconvergence throws an exception rather than returning an apparently valid result.

**Important limitations:** This formulation does not independently enforce electroneutrality, determine pH, support redox equilibria, precipitate solids, use a thermodynamically complete Pitzer model, or accept arbitrary reaction networks. Charge residuals are reported for diagnostic purposes. For physically meaningful electrolyte cases, the supplied component totals must already be charge-consistent.

## Build and test

C++17 compiler and CMake 3.16 or newer:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build
ctest --test-dir build --output-on-failure
./build/thermo ideal examples/association.txt
```

The tests compare results to two closed-form synthetic equilibrium cases: A + B ⇌ AB and 2A ⇌ A2. They also cover invalid inputs, charge inconsistency, deliberate charge imbalance, and nonconvergence behavior.

## Ideal model file

```text
# component NAME TOTAL_MOL_PER_L ELECTRIC_CHARGE
component A+ 0.01 1
component B- 0.01 -1
# complex NAME ELECTRIC_CHARGE LOG10_FORMATION_CONSTANT COUNT PRIMARY STOICH ...
complex AB 0 2.0 2 A+ 1 B- 1
```

All complex stoichiometries must reference declared primary components, use positive integer coefficients, and conserve electric charge. Each species name must be unique. All totals must be strictly positive and represent mol/L. No zero-total components are supported by this logarithmic implementation.

## Other experimental components

- Species and solid record parsers (`species.cpp`)
- Ionic-strength calculation and an **illustrative incomplete** Pitzer-inspired approximation (`pitzer.cpp`)
- Saturation-index calculation from specified activities (`solver.cpp`)

`database.txt` is currently empty. The supplied `pitzer_params.txt` and `solids.txt` contain **unverified example data**; do not use these values for a quantitative chemical engineering calculation without establishing provenance and physical validity. Neither activity-coefficient nor solid-equilibrium functions are coupled to the ideal solver.

## Roadmap

1. Add benchmark equilibria with traceable literature thermodynamic data and standard-state conventions.
2. Support explicit reactions with a validated independent-component basis, including water equilibrium.
3. Introduce a defensible electroneutrality/pH option and robust zero-component handling.
4. Add analytic or automatic differentiation and conditioning diagnostics for general networks.
5. Implement and validate activity models separately before coupling nonideal activities to speciation.
6. Add solid-phase stability and phase appearance/disappearance only after aqueous models are validated.

**Not validated for design, regulatory, nuclear safety, or process engineering decisions.**
