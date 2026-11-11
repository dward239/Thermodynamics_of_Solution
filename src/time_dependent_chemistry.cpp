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
constexpr double gasConstantJPerMolK=8.31446261815324;

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
    std::set<std::string> phaseIds;
    for(const auto& phase:model.equilibriumTemplate.phases) {
        if(phase.id.empty() || !phaseIds.insert(phase.id).second ||
           !std::isfinite(phase.volumeM3) || phase.volumeM3<=0 ||
           !std::isfinite(phase.standardConcentrationMolM3) ||
           phase.standardConcentrationMolM3<=0)
            throw std::invalid_argument("Invalid/duplicate fluid phase in chemistry template");
        if(phase.kind!=IdealFluidPhaseKind::AqueousLiquid &&
           phase.kind!=IdealFluidPhaseKind::OrganicLiquid &&
           phase.kind!=IdealFluidPhaseKind::Gas)
            throw std::invalid_argument("Unknown ideal fluid phase kind");
    }
    std::set<std::string> speciesNames;
    for(const auto& species:model.equilibriumTemplate.species) {
        if(species.id.empty() || !speciesNames.insert(species.id).second ||
           !phaseIds.count(species.phaseId) ||
           !std::isfinite(species.standardMuOverRT) ||
           species.componentStoichiometry.empty())
            throw std::invalid_argument("Invalid/duplicate fluid species in chemistry template");
    }
    std::set<std::string> equilibriumNames=speciesNames;
    for(const auto& solid:model.candidateSolids) {
        if(solid.id.empty() || !equilibriumNames.insert(solid.id).second ||
           !std::isfinite(solid.standardMuOverRT) ||
           solid.componentStoichiometry.empty())
            throw std::invalid_argument("Invalid/duplicate pure solid candidate");
    }
    if(model.thermodynamics) {
        const auto& thermo=*model.thermodynamics;
        if(!std::isfinite(thermo.minimumTemperatureK) ||
           !std::isfinite(thermo.maximumTemperatureK) ||
           !std::isfinite(thermo.referenceTemperatureK) ||
           thermo.minimumTemperatureK<=0.0 ||
           thermo.maximumTemperatureK<=thermo.minimumTemperatureK ||
           thermo.referenceTemperatureK<thermo.minimumTemperatureK ||
           thermo.referenceTemperatureK>thermo.maximumTemperatureK)
            throw std::invalid_argument("Invalid chemical thermodynamics temperature range");
        const auto& enthalpies=thermo.standardEnthalpyJPerMolById;
        for(const auto& [id,enthalpy]:enthalpies)
            if(!equilibriumNames.count(id) || !std::isfinite(enthalpy))
                throw std::invalid_argument("Invalid standard enthalpy for "+id);
        for(const auto& id:equilibriumNames)
            if(!enthalpies.count(id))
                throw std::invalid_argument("Missing standard enthalpy for "+id);
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

void updateThermodynamicInputs(
    IdealMultiphaseModel& fluid, std::vector<IdealPureSolid>& solids,
    const ChemicalProjectionThermodynamics& thermo, double temperatureK) {
    if(!std::isfinite(temperatureK) || temperatureK<thermo.minimumTemperatureK ||
       temperatureK>thermo.maximumTemperatureK)
        throw std::domain_error("Vessel temperature outside chemical thermodynamics range: "+
                                std::to_string(temperatureK)+" K");
    const double referenceK=thermo.referenceTemperatureK;
    if(temperatureK==referenceK)
        return;

    // Integrated Gibbs-Helmholtz relation at constant standard enthalpy.
    const double inverseTemperatureDifference=(referenceK-temperatureK)/temperatureK/referenceK;
    const auto adjustPotential=[&](auto& species) {
        const double enthalpy=thermo.standardEnthalpyJPerMolById.at(species.id);
        species.standardMuOverRT+=enthalpy/gasConstantJPerMolK*inverseTemperatureDifference;
        if(!std::isfinite(species.standardMuOverRT))
            throw std::runtime_error("Nonfinite standard chemical potential for "+species.id);
    };
    for(auto& species:fluid.species)
        adjustPotential(species);
    for(auto& solid:solids)
        adjustPotential(solid);
    for(auto& phase:fluid.phases) {
        if(phase.kind!=IdealFluidPhaseKind::Gas)
            continue;
        // At fixed standard pressure, c/c°(T) equals p/p°.
        phase.standardConcentrationMolM3*=referenceK/temperatureK;
        if(!std::isfinite(phase.standardConcentrationMolM3) ||
           phase.standardConcentrationMolM3<=0.0)
            throw std::runtime_error("Invalid gas standard concentration for "+phase.id);
    }
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
    if(config.thermodynamics)
        updateThermodynamicInputs(reduced,activeSolids,*config.thermodynamics,
                                  vesselSample.state.vessel.temperatureK);
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
    // For the optional solid equilibrium path the core solver does not
    // populate calculatedComponentMol; compute it consistently below from
    // aqueous/organic/gas species before adding the solid inventory.
    for(const auto& sp:config.equilibriumTemplate.species) {
        if(!out.fluid.speciesMoles.count(sp.id))
            out.fluid.speciesMoles[sp.id]=0.0;
        if(!out.fluid.concentrationsMolM3.count(sp.id))
            out.fluid.concentrationsMolM3[sp.id]=0.0;
    }
    // Independent audit from returned phase species + solid inventories.
    // Compare each conserved component with its background + nuclide input.
    for(const auto& name:ids) {
        double fluidAmount=0.0;
        for(const auto& sp:config.equilibriumTemplate.species) {
            const auto it=sp.componentStoichiometry.find(name);
            if(it!=sp.componentStoichiometry.end())
                fluidAmount+=it->second*out.fluid.speciesMoles.at(sp.id);
        }
        double solidAmount=0.0;
        for(const auto& solid:config.candidateSolids) {
            const auto it=solid.componentStoichiometry.find(name);
            if(it!=solid.componentStoichiometry.end())
                solidAmount+=it->second*out.solidMoles.at(solid.id);
        }
        const double calculated=fluidAmount+solidAmount;
        out.fluid.calculatedComponentMol[name]=fluidAmount;
        const double total=out.inputComponentMol.at(name);
        out.fluid.componentResidualMol[name]=calculated-total;
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
    // Chemistry is projected inline at each accepted candidate, including
    // t0, before that candidate is committed by the adaptive controller.
    // Trial rejection never invokes this callback. Chemistry remains one-way:
    // the projected species do NOT influence time-step error estimation.
    out.vessel=integrateRadiolyticHydrogenVessel(
        initial,endTimeSeconds,network,vesselModel,integrationOptions,
        [&](const AdaptiveVesselSample& accepted) {
            out.samples.push_back(project(accepted,network,chemistry,names));
        });
    if(out.samples.size()!=out.vessel.history.size())
        throw std::runtime_error("Chemical and vessel time grids are inconsistent");
    return out;
}

} // namespace thermo
