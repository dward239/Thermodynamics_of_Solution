#include "lumped_thermal.h"
#include <cmath>
#include <stdexcept>
namespace thermo {
double advanceLumpedTemperatureK(double initialTemperatureK,
                                  double intervalSeconds,
                                  double depositedEnergyJ,
                                  double sensibleHeatFraction,
                                  const LumpedThermalProperties& p) {
    if(!std::isfinite(initialTemperatureK) || initialTemperatureK<=0.0 ||
       !std::isfinite(intervalSeconds) || intervalSeconds<0.0 ||
       !std::isfinite(depositedEnergyJ) || depositedEnergyJ<0.0 ||
       !std::isfinite(sensibleHeatFraction) || sensibleHeatFraction<0.0 ||
       sensibleHeatFraction>1.0 ||
       !std::isfinite(p.massKg) || p.massKg<=0.0 ||
       !std::isfinite(p.heatCapacityJPerKgK) || p.heatCapacityJPerKgK<=0.0 ||
       !std::isfinite(p.conductanceWPerK) || p.conductanceWPerK<0.0 ||
       !std::isfinite(p.ambientTemperatureK) || p.ambientTemperatureK<=0.0)
        throw std::invalid_argument("Invalid lumped thermal state or properties");
    if(intervalSeconds==0.0) {
        if(depositedEnergyJ!=0.0) throw std::invalid_argument("Energy deposited during zero dt");
        return initialTemperatureK;
    }
    const double capacity=p.massKg*p.heatCapacityJPerKgK;
    const double heatJ=depositedEnergyJ*sensibleHeatFraction;
    double result;
    if(p.conductanceWPerK==0.0)
        result=initialTemperatureK+heatJ/capacity;
    else {
        const double factor=std::exp(-p.conductanceWPerK*intervalSeconds/capacity);
        const double steady=p.ambientTemperatureK+
                            heatJ/(intervalSeconds*p.conductanceWPerK);
        result=steady+(initialTemperatureK-steady)*factor;
    }
    if(!std::isfinite(result) || result<=0.0)
        throw std::runtime_error("Nonphysical lumped temperature after energy balance");
    return result;
}
} // namespace thermo
