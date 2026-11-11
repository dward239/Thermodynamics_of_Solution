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
4. When temperature-dependent thermodynamics is configured, update the
   surviving standard-state inputs using the candidate vessel temperature.
5. Solve the surviving ideal-dilute, neutral multi-fluid mass-action
   equilibrium, optionally enumerating pure-solid candidates.
6. Independently reconstruct component inventories from all output
   fluid species and solid moles, checking scaled mass closure.
7. Only after successful chemical projection, allow the adaptive
   controller to commit its accepted vessel candidate. If chemistry
   fails, throw without returning a partially qualified trajectory.

A rejected adaptive *numerical* trial never invokes chemical projection.
Projection does **not** modify nuclide decay, deposited energy, temperature,
pressure, gas/H2 kinetics or adaptive time-step selection. The accepted
time grid is still set exclusively by the vessel/H2 dynamics.

## Temperature-dependent standard states

Set `TimeDependentChemistryModel::thermodynamics` to opt in. Set the
reference temperature and the inclusive minimum and maximum temperatures
over which the supplied data can be used. The defaults are 298.15 K and
273.15–323.15 K. Each species and candidate solid must have a finite
`standardEnthalpyJPerMolById` entry, including species absent at t0.
Use a common component reference basis for the enthalpies and the
template values of `standardMuOverRT`. Specify a zero enthalpy explicitly
when that is the intended value on that basis.

With constant standard enthalpy, the integrated Gibbs–Helmholtz relation
gives `mu°(T)/(RT) = mu°(Tref)/(RTref) + (H°/R)(1/T - 1/Tref)`.
The corresponding reaction constants follow the constant-enthalpy van't
Hoff relation used, for example, in
[USGS PHREEQC](https://water.usgs.gov/water-resources/software/PHREEQC/documentation/phreeqc3-html/phreeqc3-36.htm).
Liquid standard concentrations stay fixed. Each ideal gas phase uses
`c°(T) = c°(Tref) Tref/T`, so its activity `c/c°(T)` equals `p/p°` with
`p° = c°(Tref) R Tref`. These are chemical phase partial pressures on the
equilibrium model's own fixed volumes.

The same correction applies at t0 and at each candidate that passes
the vessel's numerical checks. If the candidate is outside the supplied
temperature range, the run aborts before committing it. With no
thermodynamics configuration, the template's dimensionless potentials
and standard concentrations apply unchanged at every temperature.

## Verified synthetic examples

- A half-life-10-second parent element starts at one mole of atoms, with
  no daughter seed. It produces an initially absent daughter element,
  which appears in the chemical species distribution after decay.
- A neutral solute partitions between distinct aqueous and organic
  liquids according to specified standard chemical potentials.
- Pure XY(s) precipitates as daughter Y becomes available, with fluid
  X/Y inventories and solid XY(s) collectively conserving components.
- Time-invariant background components persist independently of decay.
- Decay heating shifts liquid partitioning, ideal gas partitioning and
  pure-solid saturation in agreement with analytical solutions.
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
- The optional temperature correction assumes constant standard
  enthalpies and fixed liquid standard concentrations. Pressure
  corrections, heat capacity effects on standard potentials, and
  nonideal gas fugacity are outside this model. Physical predictions
  require thermodynamic data appropriate to the phases and range.
- The chemical component inventory is independent of the separate H2
  radiolysis reservoir. No atom balance maps the H2 consumption sink
  to explicit reaction products; user-defined chemical species are not
  allowed to imply that they physically share those H2 moles.
  Projected gas species also do not contribute to the vessel headspace
  pressure calculation.
- Nuclide chemical identity is aggregated to element totals, so
  isotope-specific partitioning/fractionation is NOT represented.
- No actual thermodynamic, radiolysis, or nuclide-reference database
  and no experimental comparison are provided by this change.

The conservation audit covers only the declared chemical-component
inventory, fluid species and pure solids. It does NOT close all nuclear,
hydrogen, vapor, oxygen or energy balances across the full vessel.

See [SysML v1-style architecture views](../architecture/sysml-v2/views/README.md)
and [requirements traceability](../architecture/sysml-v2/requirements-traceability.md).
