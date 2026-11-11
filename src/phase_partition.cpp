#include "phase_partition.h"
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
namespace thermo {
PartitionSolution solveLinearPhasePartition(
    const std::vector<PartitionPhase>& phases,
    const std::string& referencePhase,
    const std::vector<PartitionSolute>& solutes) {
    if(referencePhase.empty()||phases.empty())
        throw std::invalid_argument("Partition model requires phases and reference");
    std::map<std::string,double> volumes;
    for(const auto& phase:phases) {
        if(phase.id.empty()||!volumes.emplace(phase.id,phase.volumeM3).second||
           !std::isfinite(phase.volumeM3)||phase.volumeM3<=0)
            throw std::invalid_argument("Invalid partition phase");
    }
    if(!volumes.count(referencePhase))
        throw std::invalid_argument("Reference phase does not exist");
    PartitionSolution result;
    std::set<std::string> seen;
    for(const auto& solute:solutes) {
        if(solute.id.empty()||!seen.insert(solute.id).second||
           !std::isfinite(solute.totalMoles)||solute.totalMoles<0)
            throw std::invalid_argument("Invalid partition solute");
        if(solute.concentrationRatioToReference.size()!=phases.size())
            throw std::invalid_argument("Every phase needs a partition ratio");
        double denominator=0.0;
        for(const auto& phase:phases) {
            const auto it=solute.concentrationRatioToReference.find(phase.id);
            if(it==solute.concentrationRatioToReference.end()||
               !std::isfinite(it->second)||it->second<=0)
                throw std::invalid_argument("Invalid/missing partition coefficient");
            if(phase.id==referencePhase && std::abs(it->second-1.0)>1e-12)
                throw std::invalid_argument("Reference partition coefficient must equal one");
            denominator+=it->second*phase.volumeM3;
        }
        if(!std::isfinite(denominator)||denominator<=0)
            throw std::invalid_argument("Invalid partition normalization");
        const double cRef=solute.totalMoles/denominator;
        if(!std::isfinite(cRef))throw std::overflow_error("Partition concentration overflow");
        for(const auto& phase:phases) {
            const double concentration=cRef*solute.concentrationRatioToReference.at(phase.id);
            const double amount=concentration*phase.volumeM3;
            if(!std::isfinite(concentration)||!std::isfinite(amount))
                throw std::overflow_error("Partition result overflow");
            result.concentrationMolM3[solute.id][phase.id]=concentration;
            result.moles[solute.id][phase.id]=amount;
        }
    }
    return result;
}
} // namespace thermo
