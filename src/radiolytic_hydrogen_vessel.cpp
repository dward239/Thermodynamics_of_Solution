#include "radiolytic_hydrogen_vessel.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace thermo {
RadiolyticHydrogenVesselStep advanceRadiolyticHydrogenVessel(
    const RadiolyticHydrogenVesselState& initial,
    double dt,
    const DecayNetwork& network,
    const RadiolyticHydrogenVesselModel& model) {
    if(model.liquidPhaseId.empty())
        throw std::invalid_argument("Missing radiolytic H2 source liquid phase");
    if(!std::isfinite(initial.vessel.headspaceGasMol)||
       !std::isfinite(initial.hydrogen.headspaceMol)||
       initial.hydrogen.headspaceMol<0.0||
       initial.hydrogen.headspaceMol>
           initial.vessel.headspaceGasMol+1e-12*
           std::max(1.0,initial.vessel.headspaceGasMol))
        throw std::invalid_argument("Tracked headspace H2 exceeds all gas");
    unsigned h2YieldCount=0;
    for(const auto& yield:model.vessel.radiolysisYields) {
        if(yield.phaseId==model.liquidPhaseId&&yield.speciesId=="H2")
            ++h2YieldCount;
        else if(yield.moleculesPer100eV>0.0)
            throw std::invalid_argument(
                "Only designated liquid H2 yield is supported in coupled gas model");
    }
    if(h2YieldCount!=1)
        throw std::invalid_argument("Exactly one designated gross H2 yield is required");
    // Obtain gross H2 from phase-absorbed emitted radiation energy.
    const auto decay=network.advance(initial.vessel.nuclideAtoms,dt);
    const auto gross=grossRadiolyticProductionMol(
        decay.depositedEnergyJByPhase,model.vessel.radiolysisYields);
    const auto key=model.liquidPhaseId+"::H2";
    const auto it=gross.find(key);
    if(it==gross.end())
        throw std::runtime_error("Designated H2 gross production missing");
    RadiolyticHydrogenVesselStep out;
    out.hydrogenStep=advanceHydrogenKinetics(
        initial.hydrogen,dt,it->second,model.kinetics);
    // This original vessel kernel is deterministic; it independently computes
    // the same decay step for thermal accounting, but energy is not added twice.
    out.vesselStep=advanceClosedVessel(initial.vessel,dt,network,model.vessel,
                                      out.hydrogenStep.netTransferToHeadspaceMol);
    out.state.vessel=out.vesselStep.state;
    out.state.hydrogen=out.hydrogenStep.state;
    const double expectedInertGas=
        initial.vessel.headspaceGasMol-initial.hydrogen.headspaceMol;
    const double actualInertGas=
        out.state.vessel.headspaceGasMol-out.state.hydrogen.headspaceMol;
    if(!std::isfinite(expectedInertGas)||!std::isfinite(actualInertGas)||
       std::abs(expectedInertGas-actualInertGas)>
           1e-10*std::max(1.0,expectedInertGas))
        throw std::runtime_error("Headspace inert-gas inventory not conserved");
    return out;
}
} // namespace thermo
