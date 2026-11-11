# Adaptive time-integration requirements-to-code traceability

**Scope:** restricted zero-dimensional radiolytic H2 vessel. General
chemical equilibrium and a coupled thermochemical DAE remain unimplemented.

| Requirement | Meaning | Allocation | Evidence |
| --- | --- | --- | --- |
| RT-020 | Reject locally inaccurate full-step trial using two-half-step error | `src/adaptive_vessel_integration.cpp::estimateError` | `testAdaptiveReductionAndRollback` |
| RT-021 | No accepted-state mutation on rejected trials or failures | `integrateRadiolyticHydrogenVessel` | caller-input immutability test |
| RT-022 | Advance monotonically to requested **absolute** end time | `AdaptiveVesselTrajectory` | `testStoppedAtNonzeroStart` |
| RT-023 | Audit conserved nucleus counts for one-parent-to-one-daughter nuclear transitions | `totalAtoms` + interval closure check | `testIntegratedRadiolysisAndHeating` |
| RT-024 | Audit gross H2, liquid H2, gas H2, and explicit consumed pool | `HydrogenKineticsStep` + interval closure check | `testIntegratedRadiolysisAndHeating` |
| RT-025 | Bound accepted/rejected attempts and minimum dt | `AdaptiveVesselOptions` | `testGuardsAndFailureDoesNotMutate` |
| RT-026 | Record interval radiation energy, temperature and gas pressure | `AdaptiveVesselSample` | `testIntegratedRadiolysisAndHeating` |
| RT-027 | Convert accepted nuclide inventories to elemental moles and add background inventories | `TimeDependentChemistryModel` and `elementalMolesFromNuclides` | `phaseAndDaughter` |
| RT-028 | Suppress absent daughter species/solids without injecting positive seed concentrations | `src/time_dependent_chemistry.cpp::project` | `phaseAndDaughter`, `allZero`, `dynamicSolidAppearance` |
| RT-029 | Execute chemistry only for accepted candidates before commit; abort on chemical failure | `AdaptiveVesselSample` acceptance callback | `testAcceptedObserverAndFailure`, `abortWhenNewDaughterHasNoSpecies` |
| RT-030 | Independently reconstruct declared conserved components across fluids and solids | `TimeDependentChemistrySample.maximumScaledComponentClosure` | `dynamicSolidAppearance` |
| RT-031 | Reject charged and invalid/unsupported model templates | `src/time_dependent_chemistry.cpp::checkModel` | `invalidModels` |

## Numerical acceptance and conservation

For each independently compared state variable `x_j`, the controller computes
`e_j = |x_full,j - x_half,j| / (atol_j + rtol*max(|x_full,j|, |x_half,j|))`
and accepts only if `max_j e_j <= 1`. The two-half-step state is
retained, not the one-full-step candidate. H2 kinetics use backward Euler,
so this is a *local discretization* estimate, not an error guarantee for
the physical H2 yield or gas-transfer model.

Within the accepted interval, C++ checks:
- `sum(N_after) - sum(N_before)` for the **restricted nuclear decay network**
- `H2_initial + H2_gross - H2_final - H2_consumed` in moles H2.

The bookkeeping checks have finite-precision thresholds: the C++ code uses
`1e-9 * max(1, reference inventory)`. These are not exact-equality
requirements and do not establish full atomic conservation when radiolytic
products are omitted. The SysML v2 requirements are conceptual and are
not interpreted by an automated verification engine yet.

**Unfulfilled:** RT-007 general nonideal multiphase phase equilibrium;
RT-008 general equilibrium-kinetics-energy DAE coupling (only a restricted
subset is now implemented); RT-009 physical and experimental validation.

## Restricted chemical projection (not DAE)

RT-027 through RT-031 are implemented only for neutral ideal-dilute
fluid chemistry with <=8 initially absent pure solids. The accepted-sample
observer is called **before commit**, but projected species do not affect
the adaptive local-error criterion or any vessel derivatives. The input
nuclear/thermal/radiolysis state is one-way; chemistry is re-equilibrated
as an algebraic diagnostic. At each accepted point, the chemical component
balance is reconstructed from species amounts plus solid amounts.

See [model limits](../../docs/time-dependent-chemical-projection.md).
