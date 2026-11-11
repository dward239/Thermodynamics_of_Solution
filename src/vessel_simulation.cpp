#include "vessel_simulation.h"
#include <cmath>
#include <stdexcept>
namespace thermo {
VesselStepResult advanceClosedVessel(const VesselState& initial,
    double dtSeconds,const DecayNetwork& network,const VesselModel& model,
    double netGasAddedMol) {
    constexpr double R=8.31446261815324;
    if(!std::isfinite(initial.timeSeconds)||initial.timeSeconds<0 ||
       !std::isfinite(initial.temperatureK)||initial.temperatureK<=0 ||
       !std::isfinite(initial.headspaceGasMol)||initial.headspaceGasMol<0 ||
       !std::isfinite(model.headspaceVolumeM3)||model.headspaceVolumeM3<=0 ||
       !std::isfinite(model.sensibleHeatFraction)||
       model.sensibleHeatFraction<0||model.sensibleHeatFraction>1 ||
       !std::isfinite(netGasAddedMol)||!std::isfinite(dtSeconds)||dtSeconds<0 ||
       !std::isfinite(initial.timeSeconds+dtSeconds))
        throw std::invalid_argument("Invalid vessel input or time step");
    const double gas=initial.headspaceGasMol+netGasAddedMol;
    if(!std::isfinite(gas)||gas<0)
        throw std::invalid_argument("Invalid headspace inventory after net transfer");
    VesselStepResult out;
    out.decay=network.advance(initial.nuclideAtoms,dtSeconds);
    out.grossRadiolysisMolByPhaseSpecies=
        grossRadiolyticProductionMol(out.decay.depositedEnergyJByPhase,
                                     model.radiolysisYields);
    double deposited=0.0;
    for(const auto& [phase,energy]:out.decay.depositedEnergyJByPhase)
        deposited+=energy;
    if(!std::isfinite(deposited))
        throw std::runtime_error("Deposited energy overflow");
    out.state=initial;
    out.state.timeSeconds+=dtSeconds;
    out.state.nuclideAtoms=out.decay.finalAtoms;
    // Gross G-value output is diagnostic. Only explicit NET transfer changes gas.
    out.state.headspaceGasMol=gas;
    out.state.temperatureK=advanceLumpedTemperatureK(
        initial.temperatureK,dtSeconds,deposited,
        model.sensibleHeatFraction,model.thermal);
    out.headspacePressurePa=gas*R*out.state.temperatureK/model.headspaceVolumeM3;
    if(!std::isfinite(out.headspacePressurePa)||out.headspacePressurePa<0)
        throw std::runtime_error("Nonphysical ideal gas headspace pressure");
    return out;
}
} // namespace thermo
