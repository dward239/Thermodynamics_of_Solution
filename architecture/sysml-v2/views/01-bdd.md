# BDD-like view — software component definitions

This view uses the **SysML v1 Block Definition Diagram** visual vocabulary.
The formal source counterparts are the `part def` elements in
[RadioThermochemicalArchitecture.sysml](../RadioThermochemicalArchitecture.sysml)
and [TimeIntegrationConnections.sysml](../TimeIntegrationConnections.sysml).

```mermaid
classDiagram
direction TB
class AdaptiveVesselIntegrator {
  <<block>>
  +integrateRadiolyticHydrogenVessel()
  +AdaptiveVesselOptions
  +AdaptiveVesselTrajectory
}
class RadiolyticHydrogenVessel {
  <<block>>
  +advanceRadiolyticHydrogenVessel()
}
class NuclearDecayKernel {
  <<block>>
  +advance()
}
class RadiationSourceKernel {
  <<block>>
  +grossRadiolyticProductionMol()
}
class HydrogenKineticsKernel {
  <<block>>
  +advanceHydrogenKinetics()
}
class ThermalKernel {
  <<block>>
  +advanceLumpedTemperatureK()
}
class PressureKernel {
  <<block>>
  +idealGasPressure()
}
class IdealMultiphaseKernel {
  <<block>>
  +solveIdealMultiphase()
}
AdaptiveVesselIntegrator --> RadiolyticHydrogenVessel : evaluates full and half steps
RadiolyticHydrogenVessel --> NuclearDecayKernel : decays nuclides
RadiolyticHydrogenVessel --> RadiationSourceKernel : computes gross H2
RadiolyticHydrogenVessel --> HydrogenKineticsKernel : transfers and consumes H2
RadiolyticHydrogenVessel --> ThermalKernel : predicts temperature
RadiolyticHydrogenVessel --> PressureKernel : updates pressure
AdaptiveVesselIntegrator ..> IdealMultiphaseKernel : planned DAE chemistry coupling
```

**Implemented:** solid-line dependencies; dotted dependence to the ideal
multiphase kernel is a **future integration**, not a currently invoked solver.
This is a software architecture BDD, not a physical vessel layout.
