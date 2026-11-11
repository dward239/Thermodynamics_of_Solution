#include "multiphase_system.h"
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
namespace thermo {
void MultiphaseSystem::validate() const {
    if (phases.empty()) throw std::invalid_argument("At least one phase is required");
    std::set<std::string> phaseIds;
    for (const auto& p: phases) {
        if (p.id.empty() || !phaseIds.insert(p.id).second)
            throw std::invalid_argument("Duplicate/empty phase id: " + p.id);
        if (!std::isfinite(p.volumeM3) || p.volumeM3 <= 0.0 ||
            !std::isfinite(p.temperatureK) || p.temperatureK <= 0.0)
            throw std::invalid_argument("Invalid phase volume or temperature: " + p.id);
    }
    std::map<std::string, const ChemicalSpecies*> speciesById;
    for (const auto& s: species) {
        if (s.id.empty() || !speciesById.emplace(s.id, &s).second)
            throw std::invalid_argument("Duplicate/empty species id: " + s.id);
        if (!phaseIds.count(s.phaseId))
            throw std::invalid_argument("Species references unknown phase: " + s.id);
        if (s.composition.empty())
            throw std::invalid_argument("Species missing composition: " + s.id);
        for (const auto& [component,count]: s.composition)
            if (component.empty() || count <= 0)
                throw std::invalid_argument("Invalid species composition: " + s.id);
    }
    std::set<std::string> reactionIds;
    for (const auto& r: reactions) {
        if (r.id.empty() || !reactionIds.insert(r.id).second ||
            !std::isfinite(r.logKReference) || r.coefficients.size() < 2)
            throw std::invalid_argument("Invalid equilibrium reaction: " + r.id);
        long long chargeBalance=0;
        bool hasReactant=false, hasProduct=false;
        std::map<std::string, long long> componentBalance;
        for (const auto& [speciesName,stoich]: r.coefficients) {
            auto it=speciesById.find(speciesName);
            if (it==speciesById.end() || stoich==0)
                throw std::invalid_argument("Missing/zero reaction species: " + speciesName);
            hasReactant |= stoich < 0;
            hasProduct |= stoich > 0;
            chargeBalance += static_cast<long long>(stoich)*it->second->charge;
            for (const auto& [component,count]:it->second->composition)
                componentBalance[component]+=static_cast<long long>(stoich)*count;
        }
        if (!hasReactant || !hasProduct || chargeBalance != 0)
            throw std::invalid_argument("Reaction missing side or violates charge: " + r.id);
        for (const auto& [component,balance]:componentBalance)
            if (balance != 0)
                throw std::invalid_argument("Reaction violates conservation of " + component +
                                            ": " + r.id);
    }
}
} // namespace thermo
