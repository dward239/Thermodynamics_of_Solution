# SD-like view — one adaptive trial

This is the **SysML v1 Sequence Diagram** analogue. It shows calls and
results rather than implementing UML Interaction objects.

```mermaid
sequenceDiagram
  autonumber
  participant Client
  participant Adaptive as AdaptiveVesselIntegrator
  participant Step as RadiolyticHydrogenVessel
  participant Decay as DecayNetwork
  participant H2 as HydrogenKinetics
  participant Thermal as VesselThermalAndPressure
  Client->>Adaptive: integrate(initial, endTime, options)
  Adaptive->>Adaptive: validate and select dt
  loop While t < endTime
    Adaptive->>Step: advance(accepted, dt)
    Step->>Decay: advance(nuclides, dt)
    Decay-->>Step: daughters and deposited energy
    Step->>H2: advance(dissolved, gas, gross H2)
    H2-->>Step: next H2, consumed, net transfer
    Step->>Thermal: heat + pressure update
    Thermal-->>Step: T, P, gas moles
    Step-->>Adaptive: full-step trial
    Adaptive->>Step: advance(accepted, dt/2)
    Step-->>Adaptive: first half-step
    Adaptive->>Step: advance(first-half state, dt/2)
    Step-->>Adaptive: second half-step
    Adaptive->>Adaptive: estimate error and audit
    alt Accepted
      Adaptive->>Adaptive: commit half-step result and record history
    else Rejected
      Adaptive->>Adaptive: log attempt, reduce dt, keep accepted state
    end
  end
  Adaptive-->>Client: accepted history, diagnostics, final state
```

There is no GUI, asynchronous message queue, or separate parallel physics
workers. The calls are synchronous C++ calculations. At present the
underlying vessel step itself computes decay twice for its thermal path;
a future refactor may eliminate this duplicated internal call.

## Accepted-step chemistry call sequence

When the optional coupling wrapper is used, the adaptive controller calls
its accepted-state observer before committing a candidate. That observer
maps nuclide counts to elemental moles, removes absent components/species,
calls the ideal fluid/multiple-pure-solid equilibrium solver, and audits
mass conservation. A chemistry exception aborts. See the
[extended sequence](09-chemistry-time-coupling.md).
