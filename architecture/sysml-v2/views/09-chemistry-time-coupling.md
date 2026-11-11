# Accepted-state chemical coupling: IBD / ACT / STM / sequence

This optional variant extends the adaptive vessel controller. It is a
**one-way chemical equilibrium projection** of accepted radionuclide states,
not a fully implicit chemical/thermal/radiolysis DAE. The source architecture
is [ChemistryTimeCoupling.sysml](../ChemistryTimeCoupling.sysml).

## IBD-like connected parts and item flows

```mermaid
flowchart LR
  subgraph V["«block» Adaptive chemistry-coupled simulation"]
    C["controller: AdaptiveVesselIntegrator"]
    N["mapper: NuclideElementMapper"]
    E["equilibrium: IdealNeutralMultiphaseKernel"]
    A["audit: ComponentConservation"]
    C -->|"accepted candidate nuclide atoms"| N
    C -->|"candidate vessel temperature"| E
    N -->|"element mol + background mol"| E
    E -->|"species mol / solid mol"| A
    A -->|"pass/fail, component residual"| C
  end
  D["Nuclide decay / thermal / H2 vessel model"] -->|"candidate states"| C
  C -->|"accepted history with chemical states"| O["Trajectory"]
```

No reverse arrow exists from `equilibrium` back to vessel heat, pressure,
decay, or radiolytic gas kinetics. No charged-species flow is supported.

## Activity-like decision flow

```mermaid
flowchart TD
  A["Candidate passes vessel step-doubling check"] --> B["Audit nucleus and H2 closure"]
  B --> C{"Accepted observer configured?"}
  C -->|No| D["Commit vessel candidate"]
  C -->|Yes| E["Nuclides → elements + background"]
  E --> F["Exclude zero-inventory species and solids"]
  F --> T["Evaluate configured standard states"]
  T -->|valid| G["Solve ideal fluid / pure-solid equilibrium"]
  T -->|out of range / invalid| J
  G --> H["Audit fluid and solid component totals"]
  H --> I{"Projection succeeds?"}
  I -->|Yes| D
  I -->|No| J([Throw: no partially qualified trajectory])
  D --> K["Record accepted state and continue"]
```

Rejected vessel **numerical** trial steps never enter this flow.

## State-machine-like failure handling

```mermaid
stateDiagram-v2
  [*] --> AwaitCandidate
  AwaitCandidate --> MapInventory: accepted candidate
  MapInventory --> EvaluateStandardStates: all components mapped
  EvaluateStandardStates --> SolveEquilibrium: standard states ready
  SolveEquilibrium --> AuditComponents: converged
  AuditComponents --> ReadyToCommit: scaled closure <= 1
  MapInventory --> Failed: unsupported / missing component
  EvaluateStandardStates --> Failed: out of range / invalid
  SolveEquilibrium --> Failed: invalid / nonconverged equilibrium
  AuditComponents --> Failed: conservation failure
  ReadyToCommit --> [*]
  Failed --> [*]
```

## Sequence-like callback and numerical library interactions

```mermaid
sequenceDiagram
  autonumber
  participant Adaptive as AdaptiveVesselIntegrator
  participant Callback as AcceptedSampleObserver
  participant Bridge as NuclideElementBridge
  participant Eq as IdealMultiphaseEquilibrium
  participant Audit as ComponentMassAudit
  Adaptive->>Callback: accepted candidate (pre-commit)
  Callback->>Bridge: elementalMolesFromNuclides(atoms)
  Bridge-->>Callback: element mole totals
  Callback->>Callback: add background and suppress zero species
  Callback->>Callback: evaluate standard states at candidate temperature
  alt valid standard states, equilibrium, and closure
    Callback->>Eq: solveIdealMultiphase / WithSolids
    Eq-->>Callback: fluid species and pure-solid inventories
    Callback->>Audit: reconstruct conserved component totals
    Audit-->>Callback: closure passed
    Callback-->>Adaptive: projection accepted
    Adaptive->>Adaptive: commit candidate and history
  else invalid temperature, chemistry, or conservation
    Callback-->>Adaptive: exception
    Adaptive->>Adaptive: abort without returning partial trajectory
  end
```

**Verification:** `tests/time_dependent_chemistry_tests.cpp` and
`tests/adaptive_vessel_tests.cpp`. These verify a numerical interface,
not scientifically validated radioisotope phase chemistry.
