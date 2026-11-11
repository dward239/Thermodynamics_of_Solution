# Restricted radiolytic hydrogen accumulation and gas–liquid kinetics

This engineering approximation adds one **hydrogen-specific** kinetic model
to the zero-dimensional radioactive vessel foundation. It represents
dissolved H2, gas-phase H2, and an explicit first-order dissolved H2
consumption sink. It is not a full radical radiolysis mechanism.

## Equations and units

Let L and G be dissolved and gas-phase H2 in moles, respectively.
Gross radiolytic production is computed from an independently specified
G(H2), deposited radiation energy and the existing gross radiolysis module.
An idealized first-order model uses

  dL/dt = q_gross - k_c * L - k_t * (L - r*G)
  dG/dt = k_t * (L - r*G)

where rates k_c and k_t are in s^-1. The dimensionless user-input
r is the equilibrium dissolved-to-headspace **mole ratio**. For a
linear concentration partition model, r=(c_liquid/c_gas)*(V_liquid/V_gas).
Consumption removes modeled H2 moles into an **unmodeled product pool**.
Thus it is a loss from the H2 inventory, not destruction of hydrogen atoms.

Backward Euler with user-selected substeps (1..100000) solves the
coupled 2x2 reservoir system, retaining nonnegative inventories and
auditing the identity

  initial L + initial G + gross H2
    = final L + final G + H2 consumed to unmodeled products.

In the integrated vessel step, transfer of the gas-phase H2 quantity
modifies the rigid-headspace ideal-gas pressure. Original inert gas
inventory is audited separately. Existing radioactive decay, prescribed
radiation-energy absorption and lumped thermal calculations are reused.

## Interpretation and exclusions

- Gross G(H2) is user-provided and treated as a constant yield over the
  step, not measured or computed from the radiation energy spectrum.
- An interval's integrated gross generation is assumed uniformly
  distributed across the numerical substeps; large intervals may have
  significant source-averaging errors.
- Recombination is a user-supplied pseudo-first-order sink, not an
  explicit H/O radical reaction network. Reaction products and reaction
  enthalpy are NOT added to the aqueous material or energy balances.
- Liquid/gas partition ratio and rates are fixed during each step.
  Henry-law data, temperature and pressure dependence, multiple liquids,
  varying headspace volume and nonideal gases are not supported.
- Pressure is obtained from a single, uniform vessel temperature and
  ideal gas law. Other radiolytic products are explicitly rejected.
- No feedback from this gas chemistry to nuclide speciation, activity
  coefficients, water composition or solid equilibria has been implemented.
- No experimental radiolysis or nuclear-waste benchmark data have been
  supplied or used. All tests use invented nuclear parameters.

**Do not use this prototype for safety or regulatory predictions.**
