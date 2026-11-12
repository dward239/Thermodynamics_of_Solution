#include "time_dependent_chemistry.h"
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
using namespace thermo;
namespace {
template<class T> void field(std::istringstream& in,T& value,const char* name,int line){
    if(!(in>>value))throw std::runtime_error(
        std::string("Missing or invalid ")+name+" at line "+std::to_string(line));
}
void trailing(std::istringstream& in,int line){
    std::string extra;
    if(in>>extra)throw std::runtime_error("Extra data at line "+std::to_string(line));
}
unsigned positiveCount(std::istringstream& in,const char* name,int line){
    int value=0;
    field(in,value,name,line);
    if(value<1 || value>100000)
        throw std::runtime_error(std::string(name)+" must be 1 through 100000 at line "+std::to_string(line));
    return static_cast<unsigned>(value);
}
std::map<std::string,int> stoichiometry(std::istringstream& in,int line){
    std::map<std::string,int> result;
    std::string component;
    while(in>>component){
        int count=0;
        field(in,count,"component coefficient",line);
        if(count<=0 || !result.emplace(component,count).second)
            throw std::runtime_error("Invalid or duplicate component term at line "+std::to_string(line));
    }
    if(result.empty())
        throw std::runtime_error("Missing component terms at line "+std::to_string(line));
    return result;
}
struct Input {
    std::vector<NuclideDefinition> nuclides;
    std::map<std::string,std::vector<DaughterBranch>> daughters;
    std::map<std::string,std::map<std::string,double>> depositions;
    VesselState initial;
    VesselModel vessel;
    double dt=0;
    unsigned steps=0;
    bool seenThermal=false,seenHeadspace=false,seenSimulation=false;
    bool adaptive=false;
    double endTimeSeconds=0.0;
    AdaptiveVesselOptions options;
    RadiolyticHydrogenVesselModel hydrogenModel;
    HydrogenKineticsState hydrogenInitial;
    TimeDependentChemistryModel chemistry;
    bool projectChemistry=false;
};
Input read(const std::string& path){
    std::ifstream f(path);
    if(!f)throw std::runtime_error("Cannot read vessel input: "+path);
    Input input;
    const std::set<std::string> singletonRecords={"simulate","adaptive","adaptive_tolerances",
        "adaptive_limits","hydrogen","chem_temperature"};
    std::set<std::string> seen;
    std::map<std::string,double> enthalpies;
    std::string line;
    int number=0;
    while(std::getline(f,line)){
        ++number;
        const auto comment=line.find('#');
        if(comment!=std::string::npos)line.resize(comment);
        std::istringstream in(line);
        std::string key;
        if(!(in>>key))continue;
        if(singletonRecords.count(key) && !seen.insert(key).second)
            throw std::runtime_error("Duplicate "+key+" record at line "+std::to_string(number));
        if(key=="nuclide"){
            NuclideDefinition n;
            field(in,n.id,"nuclide id",number);
            field(in,n.element,"element name",number);
            field(in,n.halfLifeSeconds,"half-life seconds",number);
            field(in,n.emittedEnergyJPerDecay,"emitted J per decay",number);
            input.nuclides.push_back(n);
        }else if(key=="daughter"){
            std::string parent;DaughterBranch b;
            field(in,parent,"parent id",number);
            field(in,b.daughter,"daughter id",number);
            field(in,b.fraction,"daughter branch fraction",number);
            input.daughters[parent].push_back(b);
        }else if(key=="deposit"){
            std::string parent,phase;double fraction;
            field(in,parent,"source nuclide",number);
            field(in,phase,"deposition phase",number);
            field(in,fraction,"deposition fraction",number);
            if(!input.depositions[parent].emplace(phase,fraction).second)
                throw std::runtime_error("Duplicate deposition key at line "+std::to_string(number));
        }else if(key=="atoms"){
            std::string id;double count;
            field(in,id,"inventory nuclide",number);
            field(in,count,"atom count",number);
            if(!input.initial.nuclideAtoms.emplace(id,count).second)
                throw std::runtime_error("Duplicate inventory at line "+std::to_string(number));
        }else if(key=="gvalue"){
            RadiolysisYield g;
            field(in,g.phaseId,"yield phase",number);
            field(in,g.speciesId,"yield species",number);
            field(in,g.moleculesPer100eV,"G molecules per 100eV",number);
            input.vessel.radiolysisYields.push_back(g);
        }else if(key=="thermal"){
            if(input.seenThermal)throw std::runtime_error("Duplicate thermal record");
            input.seenThermal=true;
            field(in,input.vessel.thermal.massKg,"mass kg",number);
            field(in,input.vessel.thermal.heatCapacityJPerKgK,"specific heat",number);
            field(in,input.vessel.thermal.conductanceWPerK,"conductance",number);
            field(in,input.vessel.thermal.ambientTemperatureK,"ambient T",number);
            field(in,input.vessel.sensibleHeatFraction,"sensible heat fraction",number);
        }else if(key=="headspace"){
            if(input.seenHeadspace)throw std::runtime_error("Duplicate headspace record");
            input.seenHeadspace=true;
            field(in,input.vessel.headspaceVolumeM3,"headspace volume m3",number);
            field(in,input.initial.headspaceGasMol,"headspace gas moles",number);
            field(in,input.initial.temperatureK,"initial T K",number);
        }else if(key=="simulate"){
            input.seenSimulation=true;
            field(in,input.dt,"time step seconds",number);
            input.steps=positiveCount(in,"step count",number);
        }else if(key=="adaptive"){
            input.adaptive=true;
            field(in,input.endTimeSeconds,"end time in seconds",number);
            field(in,input.options.initialStepSeconds,"initial step",number);
            field(in,input.options.minimumStepSeconds,"minimum step",number);
            field(in,input.options.maximumStepSeconds,"maximum step",number);
            field(in,input.options.relativeTolerance,"relative tolerance",number);
        }else if(key=="adaptive_tolerances"){
            field(in,input.options.absoluteTemperatureToleranceK,"temperature tolerance",number);
            field(in,input.options.absolutePressureTolerancePa,"pressure tolerance",number);
            field(in,input.options.absoluteMolesTolerance,"mole tolerance",number);
            field(in,input.options.absoluteAtomsTolerance,"atom tolerance",number);
        }else if(key=="adaptive_limits"){
            input.options.maximumAcceptedSteps=positiveCount(in,"accepted step limit",number);
            input.options.maximumRejectedSteps=positiveCount(in,"rejected step limit",number);
        }else if(key=="hydrogen"){
            auto& model=input.hydrogenModel;
            field(in,model.liquidPhaseId,"H2 source phase",number);
            field(in,input.hydrogenInitial.dissolvedMol,"initial dissolved H2 moles",number);
            field(in,input.hydrogenInitial.headspaceMol,"initial headspace H2 moles",number);
            field(in,model.kinetics.transferRatePerSecond,"H2 transfer rate",number);
            field(in,model.kinetics.dissolvedConsumptionRatePerSecond,"H2 consumption rate",number);
            field(in,model.kinetics.equilibriumDissolvedToGasMoleRatio,"equilibrium H2 mole ratio",number);
            model.kinetics.substeps=positiveCount(in,"H2 substeps",number);
        }else if(key=="chem_component"){
            input.projectChemistry=true;
            std::string id;
            double background=0.0;
            field(in,id,"component ID",number);
            field(in,background,"background component moles",number);
            input.chemistry.equilibriumTemplate.components.push_back({id,0.0,true});
            input.chemistry.backgroundComponentMol.emplace(id,background);
        }else if(key=="chem_phase"){
            input.projectChemistry=true;
            IdealPhase phase;
            std::string kind;
            field(in,phase.id,"phase ID",number);
            field(in,kind,"phase kind",number);
            field(in,phase.volumeM3,"phase volume",number);
            field(in,phase.standardConcentrationMolM3,"standard concentration",number);
            if(kind=="aqueous") phase.kind=IdealFluidPhaseKind::AqueousLiquid;
            else if(kind=="organic") phase.kind=IdealFluidPhaseKind::OrganicLiquid;
            else if(kind=="gas") phase.kind=IdealFluidPhaseKind::Gas;
            else throw std::runtime_error("Unknown chemistry phase kind at line "+std::to_string(number));
            input.chemistry.equilibriumTemplate.phases.push_back(phase);
        }else if(key=="chem_species"){
            input.projectChemistry=true;
            IdealSpecies species;
            field(in,species.id,"species ID",number);
            field(in,species.phaseId,"phase ID",number);
            field(in,species.standardMuOverRT,"standard chemical potential / RT",number);
            species.componentStoichiometry=stoichiometry(in,number);
            input.chemistry.equilibriumTemplate.species.push_back(std::move(species));
        }else if(key=="chem_solid"){
            input.projectChemistry=true;
            IdealPureSolid solid;
            field(in,solid.id,"solid ID",number);
            field(in,solid.standardMuOverRT,"standard chemical potential / RT",number);
            solid.componentStoichiometry=stoichiometry(in,number);
            input.chemistry.candidateSolids.push_back(std::move(solid));
        }else if(key=="chem_temperature"){
            input.projectChemistry=true;
            auto& data=input.chemistry.thermodynamics.emplace();
            field(in,data.referenceTemperatureK,"reference temperature",number);
            field(in,data.minimumTemperatureK,"minimum temperature",number);
            field(in,data.maximumTemperatureK,"maximum temperature",number);
        }else if(key=="chem_enthalpy"){
            input.projectChemistry=true;
            std::string id;
            double enthalpy=0.0;
            field(in,id,"species or solid ID",number);
            field(in,enthalpy,"standard enthalpy in J/mol",number);
            if(!enthalpies.emplace(id,enthalpy).second)
                throw std::runtime_error("Duplicate enthalpy at line "+std::to_string(number));
        }else{
            throw std::runtime_error("Unknown input keyword at line "+std::to_string(number)+": "+key);
        }
        trailing(in,number);
    }
    if(f.bad()) throw std::runtime_error("Error reading vessel input: "+path);
    if(!input.seenHeadspace||!input.seenThermal||input.nuclides.empty()||input.initial.nuclideAtoms.empty())
        throw std::runtime_error("Missing headspace, thermal, nuclide or inventory records");
    if(input.seenSimulation==input.adaptive)
        throw std::runtime_error("Specify exactly one simulate or adaptive record");
    if(input.seenSimulation && (!std::isfinite(input.dt) || input.dt<=0))
        throw std::runtime_error("Fixed step must be finite and positive");
    if(input.adaptive && !seen.count("hydrogen"))
        throw std::runtime_error("An adaptive run requires a hydrogen record");
    if(!input.adaptive && (seen.count("hydrogen") || seen.count("adaptive_tolerances") ||
                           seen.count("adaptive_limits") || input.projectChemistry))
        throw std::runtime_error("Hydrogen, adaptive controls and chemistry records require an adaptive run");
    if(!enthalpies.empty() && !input.chemistry.thermodynamics)
        throw std::runtime_error("chem_enthalpy requires chem_temperature");
    if(input.chemistry.thermodynamics)
        input.chemistry.thermodynamics->standardEnthalpyJPerMolById=std::move(enthalpies);
    if(input.adaptive) input.hydrogenModel.vessel=input.vessel;
    std::map<std::string,bool> valid;
    for(const auto& n:input.nuclides){
        if(!valid.emplace(n.id,true).second)throw std::runtime_error("Duplicate nuclide "+n.id);
    }
    for(const auto& [parent,children]:input.daughters)
        if(!valid.count(parent))throw std::runtime_error("Unknown decay parent "+parent);
    for(const auto& [parent,phases]:input.depositions)
        if(!valid.count(parent))throw std::runtime_error("Unknown deposition source "+parent);
    for(auto& n:input.nuclides){
        n.daughters=input.daughters[n.id];
        n.depositionFractions=input.depositions[n.id];
    }
    return input;
}
void print(std::ostream& out,double time,double temp,double pressure,double gas,
           double emitted,double absorbed,double gross){
    out<<std::setprecision(17)<<time<<','<<temp<<','<<pressure<<','
       <<gas<<','<<emitted<<','<<absorbed<<','<<gross<<'\n';
}

std::ostream& output(std::ofstream& file,const char* path){
    if(!path) return std::cout;
    file.open(path);
    if(!file) throw std::runtime_error("Cannot open output CSV");
    return file;
}

void csvColumn(std::ostream& out,const std::string& name){
    out<<",\"";
    for(char c:name){
        if(c=='"') out<<'"';
        out<<c;
    }
    out<<'"';
}

void writeAdaptive(std::ostream& out,const Input& input,const AdaptiveVesselTrajectory& trajectory,
                   const std::vector<TimeDependentChemistrySample>* chemistry=nullptr){
    out<<"time_s,step_s,normalized_step_error,temperature_K,headspace_pressure_Pa,headspace_gas_mol,"
         "dissolved_H2_mol,headspace_H2_mol,cumulative_consumed_H2_mol,interval_emitted_J,"
         "interval_absorbed_J,interval_gross_H2_mol,interval_consumed_H2_mol,"
         "interval_H2_closure_mol,interval_nuclide_closure_atoms";
    for(const auto& nuclide:input.nuclides) csvColumn(out,"atoms:"+nuclide.id);
    if(chemistry){
        out<<",maximum_scaled_component_closure";
        for(const auto& component:input.chemistry.equilibriumTemplate.components)
            csvColumn(out,"component_mol:"+component.id);
        for(const auto& species:input.chemistry.equilibriumTemplate.species)
            csvColumn(out,"species_mol:"+species.id);
        for(const auto& solid:input.chemistry.candidateSolids)
            csvColumn(out,"solid_mol:"+solid.id);
    }
    out<<'\n'<<std::setprecision(17);
    for(std::size_t i=0;i<trajectory.history.size();++i){
        const auto& s=trajectory.history[i];
        const auto& vessel=s.state.vessel;
        const auto& h2=s.state.hydrogen;
        out<<vessel.timeSeconds<<','<<s.acceptedStepSeconds<<','<<s.estimatedNormalizedError<<','
           <<vessel.temperatureK<<','<<s.pressurePa<<','<<vessel.headspaceGasMol<<','
           <<h2.dissolvedMol<<','<<h2.headspaceMol<<','<<h2.cumulativeConsumedMol<<','
           <<s.intervalEmittedEnergyJ<<','<<s.intervalDepositedEnergyJ<<','<<s.intervalGrossH2Mol<<','
           <<s.intervalConsumedH2Mol<<','<<s.intervalH2ClosureResidualMol<<','
           <<s.intervalNuclideAtomClosureResidual;
        for(const auto& nuclide:input.nuclides){
            const auto it=vessel.nuclideAtoms.find(nuclide.id);
            out<<','<<(it==vessel.nuclideAtoms.end()?0.0:it->second);
        }
        if(chemistry){
            const auto& projected=chemistry->at(i);
            out<<','<<projected.maximumScaledComponentClosure;
            for(const auto& component:input.chemistry.equilibriumTemplate.components)
                out<<','<<projected.inputComponentMol.at(component.id);
            for(const auto& species:input.chemistry.equilibriumTemplate.species)
                out<<','<<projected.fluid.speciesMoles.at(species.id);
            for(const auto& solid:input.chemistry.candidateSolids)
                out<<','<<projected.solidMoles.at(solid.id);
        }
        out<<'\n';
    }
    out.flush();
    if(!out) throw std::runtime_error("Could not write adaptive trajectory CSV");
}
}
int main(int argc,char** argv){
    if(argc<2||argc>3){
        std::cerr<<"Usage: vessel_cli MODEL.txt [OUTPUT.csv]\n";
        return 2;
    }
    try{
        const auto input=read(argv[1]);
        const DecayNetwork network(input.nuclides);
        std::ofstream file;
        const char* outputPath=argc==3?argv[2]:nullptr;
        if(input.adaptive){
            const RadiolyticHydrogenVesselState initial{input.initial,input.hydrogenInitial};
            if(input.projectChemistry){
                const auto trajectory=integrateVesselWithChemicalProjection(
                    initial,input.endTimeSeconds,network,input.hydrogenModel,input.chemistry,input.options);
                writeAdaptive(output(file,outputPath),input,trajectory.vessel,&trajectory.samples);
            }else{
                const auto trajectory=integrateRadiolyticHydrogenVessel(
                    initial,input.endTimeSeconds,network,input.hydrogenModel,input.options);
                writeAdaptive(output(file,outputPath),input,trajectory);
            }
            return 0;
        }
        auto& out=output(file,outputPath);
        out<<"time_s,temperature_K,headspace_pressure_Pa,headspace_gas_mol,"
               "interval_emitted_J,interval_absorbed_J,interval_gross_products_mol\n";
        auto state=input.initial;
        const double initialPressure=state.headspaceGasMol*8.31446261815324*
                                     state.temperatureK/input.vessel.headspaceVolumeM3;
        print(out,0,state.temperatureK,initialPressure,state.headspaceGasMol,0,0,0);
        for(unsigned k=0;k<input.steps;++k){
            // This demonstration deliberately sets NET gas transfer to zero.
            // Gross radiolytic products are not automatically accumulated.
            const auto next=advanceClosedVessel(state,input.dt,network,input.vessel);
            double absorbed=0,gross=0;
            for(const auto& [phase,e]:next.decay.depositedEnergyJByPhase)absorbed+=e;
            for(const auto& [name,n]:next.grossRadiolysisMolByPhaseSpecies)gross+=n;
            state=next.state;
            print(out,state.timeSeconds,state.temperatureK,next.headspacePressurePa,
                  state.headspaceGasMol,next.decay.emittedEnergyJ,absorbed,gross);
        }
        return 0;
    }catch(const std::exception& e){
        std::cerr<<"Vessel run failed: "<<e.what()<<'\n';
        return 1;
    }
}
