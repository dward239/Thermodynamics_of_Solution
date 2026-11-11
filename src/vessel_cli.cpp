#include "vessel_simulation.h"
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
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
struct Input {
    std::vector<NuclideDefinition> nuclides;
    std::map<std::string,std::vector<DaughterBranch>> daughters;
    std::map<std::string,std::map<std::string,double>> depositions;
    VesselState initial;
    VesselModel vessel;
    double dt=0;
    unsigned steps=0;
    bool seenThermal=false,seenHeadspace=false,seenSimulation=false;
};
Input read(const std::string& path){
    std::ifstream f(path);
    if(!f)throw std::runtime_error("Cannot read vessel input: "+path);
    Input input;
    std::string line;
    int number=0;
    while(std::getline(f,line)){
        ++number;
        const auto comment=line.find('#');
        if(comment!=std::string::npos)line.resize(comment);
        std::istringstream in(line);
        std::string key;
        if(!(in>>key))continue;
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
            if(input.seenSimulation)throw std::runtime_error("Duplicate simulate record");
            input.seenSimulation=true;
            field(in,input.dt,"time step seconds",number);
            field(in,input.steps,"number of steps",number);
        }else{
            throw std::runtime_error("Unknown input keyword at line "+std::to_string(number)+": "+key);
        }
        trailing(in,number);
    }
    if(!input.seenHeadspace||!input.seenThermal||!input.seenSimulation||
       input.nuclides.empty()||input.initial.nuclideAtoms.empty()||
       !std::isfinite(input.dt)||input.dt<=0||input.steps==0||input.steps>100000)
        throw std::runtime_error("Missing/invalid vessel, thermal, simulate or inventory records");
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
        std::ostream* out=&std::cout;
        if(argc==3){
            file.open(argv[2]);
            if(!file)throw std::runtime_error("Cannot open output CSV");
            out=&file;
        }
        *out<<"time_s,temperature_K,headspace_pressure_Pa,headspace_gas_mol,"
               "interval_emitted_J,interval_absorbed_J,interval_gross_products_mol\n";
        auto state=input.initial;
        const double initialPressure=state.headspaceGasMol*8.31446261815324*
                                     state.temperatureK/input.vessel.headspaceVolumeM3;
        print(*out,0,state.temperatureK,initialPressure,state.headspaceGasMol,0,0,0);
        for(unsigned k=0;k<input.steps;++k){
            // This demonstration deliberately sets NET gas transfer to zero.
            // Gross radiolytic products are not automatically accumulated.
            const auto next=advanceClosedVessel(state,input.dt,network,input.vessel);
            double absorbed=0,gross=0;
            for(const auto& [phase,e]:next.decay.depositedEnergyJByPhase)absorbed+=e;
            for(const auto& [name,n]:next.grossRadiolysisMolByPhaseSpecies)gross+=n;
            state=next.state;
            print(*out,state.timeSeconds,state.temperatureK,next.headspacePressurePa,
                  state.headspaceGasMol,next.decay.emittedEnergyJ,absorbed,gross);
        }
        return 0;
    }catch(const std::exception& e){
        std::cerr<<"Vessel demonstration failed: "<<e.what()<<'\n';
        return 1;
    }
}
