#pragma once
namespace thermo {
struct LumpedThermalProperties {
    double massKg=0.0;
    double heatCapacityJPerKgK=0.0;
    double conductanceWPerK=0.0;
    double ambientTemperatureK=298.15;
};
// Constant average heat deposition rate on the interval; solves C*dT/dt = Q-UA*dT.
// The user-specified heat fraction separates absorbed energy from chemical
// energy storage; this model cannot infer heat fractions from G-values.
double advanceLumpedTemperatureK(double initialTemperatureK,
                                  double intervalSeconds,
                                  double depositedEnergyJ,
                                  double sensibleHeatFraction,
                                  const LumpedThermalProperties& properties);
} // namespace thermo
