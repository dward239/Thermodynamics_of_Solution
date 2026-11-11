#include "time_dependent_chemistry.h"
#include "inventory_chemistry_bridge.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

namespace thermo {
namespace {
void requireNonnegativeFinite(double value, const std::string& label) {
    if(!std::isfinite(value) || value<0.0)
        throw std::invalid_argument("Invalid nonnegative chemical amount: "+label);
}

std::set<std::string> checkModel(
    const DecayNetwork& network, const TimeDependentChemistryModel& model) {
    if(model.equilibriumTemplate.components.empty() ||
       model.equilibriumTemplate.phases.empty() ||
       model.equilibriumTemplate.species.empty())
        throw std::invalid_argument("Chemistry template must contain components, phases and species");
    if(!model.equilibriumTemplate.electroneutralPhaseId.empty())
        throw std::invalid_argument("Charged chemical projection not supported");
    if(model.maximumEquilibriumIterations<=0 ||
       model.maximumEquilibriumIterations>100000 ||
       !std::isfinite(model.equilibriumRelativeTolerance) ||
       model.equilibriumRelativeTolerance<=0.0 ||
       model.equilibriumRelativeTolerance>=0.01 ||
       !std::isfinite(model.conservationRelativeTolerance) ||
       model.conservationRelativeTolerance<=0 ||
       !std::isfinite(model.conservationAbsoluteToleranceMol) ||
       model.conservationAbsoluteToleranceMol<=0)
        throw std::invalid_argument("Invalid chemical projection convergence settings");
    if(model.candidateSolids.size()>8)
        throw std::invalid_argument("More than eight candidate pure solids");
    std::set<std::string> ids;
    for(const auto& c:model.equilibriumTemplate.components) {
        if(c.id.empty() || !ids.insert(c.id).second || !c.fixedTotal ||
           c.totalMol!=0.0)
            throw std::invalid_argument(
                "Template components must be unique fixed-total placeholders with zero total");
    }
    for(const auto& nuclide:network.nuclides())
        if(!ids.count(nuclide.element))
            throw std::invalid_argument("Missing chemical component for nuclide element: "+nuclide.element);
    for(const auto& [name,amount]:model.backgroundComponentMol) {
        if(!ids.count(name))
            throw std::invalid_argument("Unknown background component: "+name);
        requireNonnegativeFinite(amount,name);
    }
    for(const auto& species:model.equilibriumTemplate.species) {
        if(species.charge!=0)
            throw std::invalid_argument("Charged species not supported by time projection");
        for(const auto& [component,count]:species.componentStoichiometry)
            if(!ids.count(component) || count<=0)
                throw std::invalid_argument("Unknown/invalid species component in "+species.id);
    }
    for(const auto& solid:model.candidateSolids)
        for(const auto& [component,count]:solid.componentStoichiometry)
            if(!ids.count(component) || count<=0)
                throw std::invalid_argument("Unknown/invalid solid component: "+solid.id);
    return ids;
}

TimeDependentChemistrySample project(
    const AdaptiveVesselSample& vesselSample,
    const DecayNetwork& network,
    const TimeDependentChemistryModel& config,
    const std::set<std::string>& ids) {
    TimeDependentChemistrySample out;
    out.vessel=vesselSample;
    out.inputComponentMol=config.backgroundComponentMol;
    const auto fromNuclides=elementalMolesFromNuclides(
        network,vesselSample.state.vessel.nuclideAtoms);
    for(const auto& [name,amount]:fromNuclides)
        out.inputComponentMol[name]+=amount;
    for(const auto& name:ids)
        if(!out.inputComponentMol.count(name))
            out.inputComponentMol[name]=0.0;

    IdealMultiphaseModel reduced=config.equilibriumTemplate;
    reduced.components.clear();
    std::set<std::string> active;
    for(const auto& c:config.equilibriumTemplate.components) {
        const double total=out.inputComponentMol.at(c.id);
        requireNonnegativeFinite(total,c.id);
        if(total>0.0) {
            active.insert(c.id);
            reduced.components.push_back({c.id,total,true});
        }
    }
    reduced.species.clear();
    for(const auto& sp:config.equilibriumTemplate.species) {
        bool allowed=true;
        for(const auto& [name,count]:sp.componentStoichiometry)
            if(!active.count(name))allowed=false;
        if(allowed)reduced.species.push_back(sp);
        out.fluid.speciesMoles[sp.id]=0.0;
        out.fluid.concentrationsMolM3[sp.id]=0.0;
    }

    std::vector<IdealPureSolid> activeSolids;
    for(const auto& solid:config.candidateSolids) {
        bool allowed=true;
        for(const auto& [name,count]:solid.componentStoichiometry)
            if(!active.count(name))allowed=false;
        if(allowed)activeSolids.push_back(solid);
        out.solidMoles[solid.id]=0.0;
    }
    if(!reduced.components.empty()) {
        if(reduced.species.empty())
            throw std::invalid_argument("Positive component inventory has no permissible fluid species");
        if(activeSolids.empty()) {
            out.fluid=solveIdealMultiphase(reduced,config.maximumEquilibriumIterations,
                                           config.equilibriumRelativeTolerance);
        } else {
            auto solids=solveIdealMultiphaseWithSolids(
                reduced,activeSolids,config.maximumEquilibriumIterations,
                config.equilibriumRelativeTolerance);
            out.fluid=std::move(solids.fluidSolution);
            for(const auto& [name,amount]:solids.solidMoles)
                out.solidMoles[name]=amount;
        }
    }
    // Populate absent components, species and solids with explicit zeros.
    for(const auto& c:config.equilibriumTemplate.components) {
        const auto name=c.id;
        if(!out.fluid.calculatedComponentMol.count(name))
            out.fluid.calculatedComponentMol[name]=0.0;
        if(!out.fluid.componentResidualMol.count(name))
            out.fluid.componentResidualMol[name]=0.0;
    }
    for(const auto& sp:config.equilibriumTemplate.species) {
        if(!out.fluid.speciesMoles.count(sp.id))
            out.fluid.speciesMoles[sp.id]=0.0;
        if(!out.fluid.concentrationsMolM3.count(sp.id))
            out.fluid.concentrationsMolM3[sp.id]=0.0;
    }
    // Independent audit from returned phase species + solid inventories.
    // Compare each conserved component with its background + nuclide input.
    for(const auto& name:ids) {
        double calculated=0.0;
        for(const auto& sp:config.equilibriumTemplate.species) {
            const auto it=sp.componentStoichiometry.find(name);
            if(it!=sp.componentStoichiometry.end())
                calculated+=it->second*out.fluid.speciesMoles.at(sp.id);
        }
        for(const auto& solid:config.candidateSolids) {
            const auto it=solid.componentStoichiometry.find(name);
            if(it!=solid.componentStoichiometry.end())
                calculated+=it->second*out.solidMoles.at(solid.id);
        }
        const double total=out.inputComponentMol.at(name);
        const double scale=config.conservationAbsoluteToleranceMol+
                           config.conservationRelativeTolerance*total;
        const double residual=std::abs(calculated-total)/scale;
        if(!std::isfinite(residual))
            throw std::runtime_error("Nonfinite chemical component closure: "+name);
        out.maximumScaledComponentClosure=
            std::max(out.maximumScaledComponentClosure,residual);
        // Report component totals including precipitated solids in this audit.
        if(residual>1.0)
            throw std::runtime_error("Chemical component conservation violated: "+name);
    }
    return out;
}
} // namespace

TimeDependentChemistryTrajectory integrateVesselWithChemicalProjection(
    const RadiolyticHydrogenVesselState& initial,
    double endTimeSeconds,
    const DecayNetwork& network,
    const RadiolyticHydrogenVesselModel& vesselModel,
    const TimeDependentChemistryModel& chemistry,
    const AdaptiveVesselOptions& integrationOptions) {
    const auto names=checkModel(network,chemistry);
    TimeDependentChemistryTrajectory out;
    out.vessel=integrateRadiolyticHydrogenVessel(
        initial,endTimeSeconds,network,vesselModel,integrationOptions);
    out.samples.reserve(out.vessel.history.size());
    for(const auto& accepted:out.vessel.history)
        out.samples.push_back(project(accepted,network,chemistry,names));
    if(out.samples.size()!=out.vessel.history.size())
        throw std::runtime_error("Chemical and vessel time grids are inconsistent");
    return out;
}

} // namespace thermo
