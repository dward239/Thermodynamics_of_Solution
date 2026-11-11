#include "decay_phase_partition.h"
#include "inventory_chemistry_bridge.h"
#include <stdexcept>
#include <vector>
namespace thermo {
DecayPartitionResult advanceDecayAndLinearRepartition(
    const DecayNetwork& network,
    const std::map<std::string,double>& atoms,
    double dt,
    const std::vector<PartitionPhase>& phases,
    const std::string& referencePhase,
    const std::map<std::string,std::map<std::string,double>>& ratios) {
    DecayPartitionResult result;
    result.decay=network.advance(atoms,dt);
    result.elementTotalMol=elementalMolesFromNuclides(
        network,result.decay.finalAtoms);
    std::vector<PartitionSolute> solutes;
    for(const auto& [element,amount]:result.elementTotalMol) {
        if(amount==0.0)continue;
        const auto it=ratios.find(element);
        if(it==ratios.end())
            throw std::invalid_argument(
                "Missing neutral phase-partition coefficient set for daughter/element: "+element);
        solutes.push_back({element,amount,it->second});
    }
    result.phaseDistribution=solveLinearPhasePartition(
        phases,referencePhase,solutes);
    return result;
}
} // namespace thermo
