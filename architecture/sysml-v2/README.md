# SysML v2 software architecture baseline

This directory holds the **intended** architecture. Its components distinguish
implemented mathematical kernels from planned coupling and phase behavior.
SysML v2 textual notation has not yet been independently parser-validated
in CI. Architectural agreement is a design objective, not a verified claim.

## Allocation and status

| SysML v2 definition | C++ implementation | Status |
| --- | --- | --- |
| ChemicalSystemModel | `include/multiphase_system.h` | Phase/species/reaction structural validation |
| NuclideInventoryModel | `include/decay_network.h` | Inventory per nuclide, not chemically partitioned |
| NuclearDecayModel | `src/decay_network.cpp` | First-order branching; fission unsupported |
| RadiationDepositionModel | `NuclideDefinition::depositionFractions` | User-prescribed absorption; no shielding solver |
| RadiolysisSourceModel | `src/radiolysis_sources.cpp` | Gross G-value production, not net gas |
| ThermalBalanceModel | `src/lumped_thermal.cpp` | Lumped constant-capacity energy balance |
| ChemicalEquilibriumModel | `include/aqueous_engine.h` | Existing restricted aqueous solver |
| SimulationController | Not implemented | Future coupling/integrator |

## Engineering requirement traceability

| ID | Required behavior | Verification |
| --- | --- | --- |
| RT-001 | Define arbitrary finite aqueous, organic liquid, gas and solid phase populations. | `tests/foundation_tests.cpp` |
| RT-002 | Reject species/reaction networks violating charge or component conservation. | `tests/foundation_tests.cpp` |
| RT-003 | Track radioactive half-lives, branches, daughters and cumulative decays. | `tests/foundation_tests.cpp`, Bateman benchmark |
| RT-004 | Differentiate radiation emission from energy deposited in each phase. | `tests/foundation_tests.cpp` |
| RT-005 | Convert user-supplied G-values from 100 eV events to mol per phase. | `tests/foundation_tests.cpp` |
| RT-006 | Integrate lumped sensible-heat balance with linear heat loss. | `tests/foundation_tests.cpp` |
| RT-007 | Solve general multiple-liquid, gas and solid chemical equilibria. | **Not implemented** |
| RT-008 | Couple decay, absorbed dose, chemistry, temperature, gas and pressure in time. | **Not implemented** |
| RT-009 | Validate predictions against independent nuclear/chemical reference data. | **Not implemented** |

## Critical separation-of-concerns decisions

- A **phase** is distinct from a physical **compartment**. Two immiscible
  liquids can share one vessel; two water tanks can be separate compartments.
- An isotope is a nuclear identity; a chemical species is a molecular or ionic
  form assigned to a phase. Isotope-to-species inventories require future coupling.
- Nuclear decay energy is not automatically equal to phase-deposited energy;
  deposited fractions are prescribed pending radiation transport.
- Gross radiolytic G-value output does not account for radical recombination,
  chemical reactions, gas partitioning or pressure. It is not net H2 generation.
- The heat model requires an explicit sensible-heat fraction because chemical
  enthalpy storage and phase change are not solved.
- The existing aqueous solver **remains limited to 273.15–323.15 K**, even
  if a separate energy balance predicts hotter conditions.

## Next architecture milestones

1. Isotopic inventories with chemical-element and phase allocations.
2. Multiphase thermodynamic state and active-phase equilibrium.
3. Defined vessel compartments, interfaces, mass-transfer and gas pressure.
4. A time integrator for chemical equilibrium constraints plus decay and kinetics.
5. Phase absorbed dose, spectrum/LET and detailed radiolytic chemistry.
6. Independent reference datasets and benchmark acceptance criteria.
