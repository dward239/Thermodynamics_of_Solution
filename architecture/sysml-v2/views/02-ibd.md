# IBD-like view — internal parts, ports and exchanges

Corresponds to the SysML v1 **Internal Block Diagram** (IBD) concerns:
part usages, interfaces, and item flows. Textual source:
[TimeIntegrationConnections.sysml](../TimeIntegrationConnections.sysml).

```mermaid
flowchart LR
  subgraph System["«block» AdaptiveVesselSoftwareIBD"]
    direction LR
    A["controller: AdaptiveController"]
    S["vessel: VesselStepKernel"]
    D["decay: DecayKernel"]
    H["hydrogen: RadiolysisHydrogenKernel"]
    T["thermal: ThermalBalanceKernel"]
    A -->|"StepRequest (t, dt, accepted state)"| S
    S -->|"StepResult (new state, pressure)"| A
    S -->|"NuclideState"| D
    D -->|"DaughterInventory"| S
    D -->|"DepositedRadiationEnergy"| S
    S -->|"AbsorbedEnergy, G(H2)"| H
    H -->|"Dissolved + headspace H2"| S
    S -->|"DepositedEnergy"| T
    T -->|"Temperature"| S
  end
  U["Caller: configuration and initial state"] --> A
  A --> V["Trajectory: accepted history + rejection log"]
```

| Exchange | Units / semantics | Source API |
| --- | --- | --- |
| State and step duration | seconds, K, moles, atoms | `adaptive_vessel_integration.h` |
| Nuclide state | atoms by nuclide ID | `decay_network.h` |
| Radiation deposition | J by phase | `decay_network.h` |
| Gross H2 | mol over interval, *not* net pressure gas | `radiolysis_sources.h` |
| Dissolved / gas H2 | mol per reservoir | `hydrogen_kinetics.h` |
| Thermal result | temperature K | `lumped_thermal.h` |
| Pressure | Pa, ideal gas | `vessel_simulation.h` |

These are **numerical software interfaces**, not fluid ports, pipes, or mass
transfer geometries. Connector payloads in the SysML source are intentionally
coarse **conceptual item types**, not complete generated C++ bindings.
