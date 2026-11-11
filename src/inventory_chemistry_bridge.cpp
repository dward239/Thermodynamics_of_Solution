#include "inventory_chemistry_bridge.h"
#include <cmath>
#include <set>
#include <stdexcept>
namespace thermo {
std::map<std::string,double> elementalMolesFromNuclides(
    const DecayNetwork& network,
    const std::map<std::string,double>& nuclideAtoms) {
    constexpr double avogadro=6.02214076e23;
    std::map<std::string,std::string> identity;
    for(const auto& n:network.nuclides())identity[n.id]=n.element;
    std::map<std::string,double> elemental;
    for(const auto& [id,atoms]:nuclideAtoms) {
        if(!identity.count(id) || !std::isfinite(atoms) || atoms<0)
            throw std::invalid_argument("Unknown nuclide or invalid atom count: "+id);
        elemental[identity.at(id)]+=atoms/avogadro;
        if(!std::isfinite(elemental.at(identity.at(id))))
            throw std::overflow_error("Elemental mole total overflow");
    }
    return elemental;
}
void validateElementCoverage(const MultiphaseSystem& system,
    const std::map<std::string,double>& elementalMoles) {
    system.validate();
    std::set<std::string> present;
    for(const auto& sp:system.species)
        for(const auto& [element,count]:sp.composition)present.insert(element);
    for(const auto& [element,moles]:elementalMoles) {
        if(element.empty() || !std::isfinite(moles) || moles<0)
            throw std::invalid_argument("Invalid elemental inventory: "+element);
        if(moles>0 && !present.count(element))
            throw std::invalid_argument("No chemical species defined for element: "+element);
    }
}
} // namespace thermo
