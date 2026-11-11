# Allocation view — behavior, software parts and C++ evidence

This view is the **functional allocation** concern familiar from SysML v1
«allocate» relationships. Textual source:
[TimeIntegrationAllocations.sysml](../TimeIntegrationAllocations.sysml).

```mermaid
flowchart LR
  subgraph Functional["«activity» behavior"]
    Trial["EvaluateAdaptiveTrial"]
    Good["AcceptedBranch"]
    Bad["RejectedBranch"]
  end
  subgraph Components["«block» logical software parts"]
    Controller["controller: AdaptiveController"]
    Vessel["vessel: VesselStepKernel"]
    H2["hydrogen: RadiolysisHydrogenKernel"]
  end
  subgraph Code["C++ implementation and tests"]
    A["adaptive_vessel_integration.cpp"]
    B["radiolytic_hydrogen_vessel.cpp"]
    C["hydrogen_kinetics.cpp"]
    D["adaptive_vessel_tests.cpp"]
  end
  Trial -.->|«allocate»| Controller
  Good -.->|«allocate»| Controller
  Bad -.->|«allocate»| Controller
  Controller --> A
  Vessel --> B
  H2 --> C
  A -->|«verify»| D
```

The SysML v2 `allocation` usages allocate control-flow behavior to
the logical controller part. The dotted links are *architectural mappings*,
not measured execution times or verified generated bindings. The
implementation files and C++ test are the currently available evidence.
