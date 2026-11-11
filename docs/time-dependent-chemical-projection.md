# Time-dependent neutral equilibrium projection: scientific boundaries

This is the **first coupling increment** between the accepted adaptive
radionuclide-vessel time grid and the generalized chemical-equilibrium
solvers. It is NOT a fully coupled multiphase thermochemical DAE.

## Execution contract

The new `TimeDependentChemistryModel` accepts:

- an ideal-dilute, neutral species/phase template with **zero placeholder
  component totals**, plus optional <=8 initially absent pure solids;
- independent time-invariant background chemical-component moles;
- the restricted radioactive vessel's nuclide network and accepted steps.

At t0 and before **each accepted candidate** commits, a callback performs:

1. Convert each nuclide's atom count into element moles (Avogadro's
   constant); sum isotope inventories of the same element.
2. Add fixed background component totals.
3. Remove zero-total components and their dependent species/solids from
   the *current* equilibrium solve. Return explicit zero results for
   absent species; do not invent numerical concentration seeds.
4. Solve the surviving ideal-dilute, neutral multi-fluid mass-action
   equilibrium, optionally enumerating pure-solid candidates.
5. Independently reconstruct component inventories from all output
   fluid species and solid moles, checking scaled mass closure.
6. Only after successful chemical projection, allow the adaptive
   controller to commit its accepted vessel candidate. If chemistry
   fails, throw without returning a partially qualified trajectory.

A rejected adaptive *numerical* trial never invokes chemical projection.
Projection does **not** modify nuclide decay, deposited energy, temperature,
pressure, gas/H2 kinetics or adaptive time-step selection. The accepted
time grid is still set exclusively by the vessel/H2 dynamics.

## Verified synthetic examples

- A half-life-10-second parent element starts at one mole of atoms, with
  no daughter seed. It produces an initially absent daughter element,
  which appears in the chemical species distribution after decay.
- A neutral solute partitions between distinct aqueous and organic
  liquids according to specified standard chemical potentials.
- Pure XY(s) precipitates as daughter Y becomes available, with fluid
  X/Y inventories and solid XY(s) collectively conserving components.
- Time-invariant background components persist independently of decay.
- All-zero chemical inventories remain exactly absent, not 1e-30 seeds.
- Duplicate/missing phases/species and unsupported charged chemistry
  reject explicitly. An unrepresentable daughter element appearing
  during the run aborts rather than committing an invalid sample.

## Scientific and numerical limits

- **One-way operator splitting:** projected equilibrium does not influence
  temperature, radiolysis, pressure, mass-transfer coefficients, or decay.
  It is not a converged iterative coupling, nonlinear DAE, or global
  multiphase chemical-reactor material and energy balance.
- Thermodynamics are **neutral and ideal-dilute**. Charged electrolyte
  Davies/Pitzer variants and charge-balanced solids are not coupled here.
- Phase definitions are fixed volumes. Phase appearance/disappearance
  applies only to <=8 ideal pure-solid candidates, not fluid phase stability.
- The dimensionless supplied species standard chemical potentials
  `mu0/(RT)` are **held constant over the trajectory**. Changes in
  vessel temperature do not adjust equilibrium constants or properties.
  This is thermodynamically inconsistent for general nonisothermal
  systems; all synthetic tests are an illustrative restricted regime.
  A real implementation needs validated `mu0(T,P)`/EOS functions and
  a temperature-consistent equilibrium callback.
- The chemical component inventory is independent of the separate H2
  radiolysis reservoir. No atom balance maps the H2 consumption sink
  to explicit reaction products; user-defined chemical species are not
  allowed to imply that they physically share those H2 moles.
- Nuclide chemical identity is aggregated to element totals, so
  isotope-specific partitioning/fractionation is NOT represented.
- No actual thermodynamic, radiolysis, or nuclide-reference database
  and no experimental comparison are provided by this change.

The conservation audit covers only the declared chemical-component
inventory, fluid species and pure solids. It does NOT close all nuclear,
hydrogen, vapor, oxygen or energy balances across the full vessel.

See [SysML v1-style architecture views](../architecture/sysml-v2/views/README.md)
and [requirements traceability](../architecture/sysml-v2/requirements-traceability.md).
