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
| RadiolysisSourceModel | `src/radiolysis_sources.cpp` and `src/hydrogen_kinetics.cpp` | Gross G-value and restricted net H2 accumulation model |
| ThermalBalanceModel | `src/lumped_thermal.cpp` | Lumped constant-capacity energy balance |
| ChemicalEquilibriumModel | `include/aqueous_engine.h`, `include/ideal_multiphase_equilibrium.h` | Original aqueous solver; independent ideal-dilute fluid solver with restricted single-phase charge balance and bounded neutral multiple-pure-solid extension |
| SimulationController | `src/vessel_simulation.cpp`, `src/vessel_cli.cpp`, `src/radiolytic_hydrogen_vessel.cpp`, `src/adaptive_vessel_integration.cpp` | Adaptive rollback-enabled restricted decay/H2/heat/pressure integration; not full chemistry/DAE |
| Nuclide-to-chemical-component bridge | `src/inventory_chemistry_bridge.cpp` | Per-isotope atoms to elemental moles, no phase allocation |
| Neutral linear phase partition | `src/phase_partition.cpp` | Arbitrary fixed-volume liquid/gas phases; neutral dilute species only |
| Restricted decay-to-neutral partition bridge | `src/decay_phase_partition.cpp` | Element totals from daughter nuclei; each element treated as one neutral pseudo-solute |

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
| RT-010 | Convert isotopic inventories into chemical elemental mole totals. | `tests/inventory_bridge_tests.cpp` |
| RT-011 | Distribute neutral dilute solutes among arbitrary liquid and gas phases with mass conservation. | `tests/phase_partition_tests.cpp` |
| RT-012 | Demonstrate source-separated decay heating and headspace pressure in a closed vessel. | `tests/vessel_tests.cpp` |
| RT-013 | Accept user-defined nuclide networks and time grids without recompiling code. | `examples/synthetic_vessel.txt`, `synthetic_vessel_cli` CTest |
| RT-014 | Map newly generated daughter elements to neutral multi-liquid partition inventories and reject missing parameters. | `tests/decay_phase_partition_tests.cpp` |
| RT-015 | Solve simultaneous ideal-dilute fluid mass action and component conservation with arbitrary neutral species across liquids and gases. | `tests/ideal_multiphase_tests.cpp` |
| RT-016 | Select up to eight initially absent pure solids by equilibrium saturation and complementarity. | `tests/ideal_multisolid_tests.cpp` |
| RT-017 | Support a designated electroneutral fluid phase by inferring one free counterion component total. | `tests/ideal_charge_balance_tests.cpp` |
| RT-018 | Compute a bounded, self-consistent approximate Davies ionic strength and charged-species activities in the generalized fluid model. | `tests/davies_multiphase_tests.cpp` |
| RT-019 | Track radiolytic H2 gross generation, dissolved/headspace transfer and first-order modeled loss while conserving H2 molecular inventory plus loss. | `tests/hydrogen_kinetics_tests.cpp` |
| RT-020 | Reject locally inaccurate time-integration steps and record attempts. | `tests/adaptive_vessel_tests.cpp` |
| RT-021 | Preserve accepted state across rejected trials and failure. | `tests/adaptive_vessel_tests.cpp` |
| RT-022 | Reach requested absolute end time with increasing accepted timestamps. | `tests/adaptive_vessel_tests.cpp` |
| RT-023 | Audit total nucleus count per accepted interval. | `tests/adaptive_vessel_tests.cpp` |
| RT-024 | Audit H2 inventory plus consumed sink per accepted interval. | `tests/adaptive_vessel_tests.cpp` |
| RT-025 | Enforce accepted/rejected step budgets and minimum step. | `tests/adaptive_vessel_tests.cpp` |
| RT-026 | Report interval energy, thermal response and pressure. | `tests/adaptive_vessel_tests.cpp` |

## Critical separation-of-concerns decisions

- A **phase** is distinct from a physical **compartment**. Two immiscible
  liquids can share one vessel; two water tanks can be separate compartments.
- An isotope is a nuclear identity; a chemical species is a molecular or ionic
  form assigned to a phase. Isotope-to-species inventories require future coupling.
- Nuclear decay energy is not automatically equal to phase-deposited energy;
  deposited fractions are prescribed pending radiation transport.
- Gross radiolytic G-value output does not account for radical recombination,
  chemical reactions or gas partitioning; the vessel pressure calculation
  uses *external net gas transfer*, not gross H2 production.
- The heat model requires an explicit sensible-heat fraction because chemical
  enthalpy storage and phase change are not solved.
- The existing aqueous solver **remains limited to 273.15–323.15 K**, even
  if a separate energy balance predicts hotter conditions.

## Next architecture milestones

1. Isotope-to-phase and isotope-to-chemical-species allocations (elemental totals are now available).
2. Multiphase thermodynamic state and active-phase equilibrium.
3. Defined vessel compartments, interfaces, mass-transfer and gas pressure.
4. A time integrator for chemical equilibrium constraints plus decay and kinetics.
5. Phase absorbed dose, spectrum/LET and detailed radiolytic chemistry.
6. Independent reference datasets and benchmark acceptance criteria.

## New ideal-dilute chemical-potential solver boundary

The ideal-fluid solver has arbitrary dynamic species, component and phase
counts, but currently supports **neutral species only** with supplied
dimensionless standard chemical potentials and fixed-volume phases.
It does not supersede the original ionic aqueous solver or provide
general nonideal multiphase equilibrium. See
`docs/general-ideal-fluid-equilibrium.md` for equations and assumptions.

### Candidate pure-solid phase stability

The optional ideal-fluid extension enumerates active subsets of at most eight
pure solids, requiring positive inventories, component conservation and
nonpositive affinity for absent solids. It is an experimental first step and
is not a production solid-activity model. See
`docs/ideal-multisolid-equilibrium.md`.

### Restricted charged-species conservation

`ideal_multiphase_equilibrium` now supports exactly one designated charged
fluid phase whose electroneutrality replaces one component mass balance;
its free component total is inferred. This does not implement multi-phase
electrochemical potentials, nonideal electrolytes or charged-solid selection.
See `docs/general-ideal-charge-balance.md`.

### Approximate nonideal extension

An optional `davies_multiphase_equilibrium` adapter now iterates ionic
strength and charged-aqueous-species activity corrections with the original
molarity-based Davies function. It rejects unsupported temperature and high
ionic strength; molality-correct Pitzer and multi-liquid electrochemical
potentials remain unimplemented. See `docs/davies-multiphase-equilibrium.md`.

### Hydrogen kinetic coupling limit

The new `radiolytic_hydrogen_vessel` tracks two H2 reservoirs and updates
pressure using net headspace transfer. The consumption sink feeds an
**unmodeled product pool** and no explicit atom or chemical enthalpy balance;
the model is not closed thermochemical reactor simulation. See
`docs/radiolytic-hydrogen-kinetics.md`.

## SysML v1-style diagram views of SysML v2 concepts

GitHub-renderable [diagram views](views/README.md) are organized by BDD,
IBD, ACT, STM, SD, REQ, and PAR concerns. Textual SysML v2 files include
connected `part` and `port` definitions, `action` and `state` definitions,
and software requirement definitions. See `TimeIntegrationConnections.sysml`,
`TimeIntegrationBehavior.sysml`, and `TimeIntegrationRequirements.sysml`.

**Validation status:** no SysML v2 parser or Mermaid renderer runs in CI.
These are engineering model views, *not* proven executable model links.
Connections express software data exchanges, not pipes or liquid interfaces.
See [requirement traceability](requirements-traceability.md) and the
[adaptive time-integration scope](../../docs/adaptive-time-integration.md).
