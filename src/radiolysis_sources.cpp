#include "radiolysis_sources.h"
#include <cmath>
#include <stdexcept>
namespace thermo {
std::map<std::string,double> grossRadiolyticProductionMol(
    const std::map<std::string,double>& depositedEnergyJByPhase,
    const std::vector<RadiolysisYield>& yields) {
    constexpr double elementaryChargeJPerEV=1.602176634e-19;
    constexpr double avogadroPerMol=6.02214076e23;
    constexpr double joulesPerMolOf100EV=100.0*elementaryChargeJPerEV*avogadroPerMol;
    std::map<std::string,double> output;
    for(const auto& [phase,energy]:depositedEnergyJByPhase)
        if(phase.empty() || !std::isfinite(energy) || energy<0.0)
            throw std::invalid_argument("Invalid deposited energy phase or value");
    for(const auto& y:yields) {
        if(y.phaseId.empty() || y.speciesId.empty() ||
           !std::isfinite(y.moleculesPer100eV) || y.moleculesPer100eV<0)
            throw std::invalid_argument("Invalid radiolytic G-value definition");
        auto it=depositedEnergyJByPhase.find(y.phaseId);
        if(it==depositedEnergyJByPhase.end())
            throw std::invalid_argument("Radiolysis yield refers to phase without deposited energy");
        const double moles=it->second*y.moleculesPer100eV/joulesPerMolOf100EV;
        if(!std::isfinite(moles)) throw std::overflow_error("Radiolytic production overflow");
        output[y.phaseId+"::"+y.speciesId]+=moles;
    }
    return output;
}
} // namespace thermo
