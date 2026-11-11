#include "decay_phase_partition.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace thermo;
static void closeTo(double a,double b,double t,const char* msg){
    if(!std::isfinite(a)||std::abs(a-b)>t*std::max(1.0,std::abs(b)))
        throw std::runtime_error(msg);
}
int main(){
    constexpr double NA=6.02214076e23;
    const DecayNetwork network({{"X","ElementX",10,0,{{"Y",1.0}},{}},
                                {"Y","ElementY",0,0,{},{}}});
    const std::vector<PartitionPhase> phases={{"aqueous",1.0},{"organic",1.0}};
    const std::map<std::string,std::map<std::string,double>> ratios={
        {"ElementX",{{"aqueous",1.0},{"organic",3.0}}},
        {"ElementY",{{"aqueous",1.0},{"organic",1.0}}}
    };
    auto r=advanceDecayAndLinearRepartition(
        network,{{"X",NA}},10.0,phases,"aqueous",ratios);
    closeTo(r.elementTotalMol.at("ElementX"),0.5,1e-11,"Surviving parent element");
    closeTo(r.elementTotalMol.at("ElementY"),0.5,1e-11,"New daughter element");
    closeTo(r.phaseDistribution.moles.at("ElementX").at("aqueous"),
            0.125,1e-11,"X aqueous partition");
    closeTo(r.phaseDistribution.moles.at("ElementX").at("organic"),
            0.375,1e-11,"X organic partition");
    closeTo(r.phaseDistribution.moles.at("ElementY").at("aqueous"),
            0.25,1e-11,"Y aqueous partition");
    closeTo(r.phaseDistribution.moles.at("ElementY").at("organic"),
            0.25,1e-11,"Y organic partition");
    double total=0;
    for(const auto& [element,p]:r.phaseDistribution.moles)
        for(const auto& [phase,n]:p)total+=n;
    closeTo(total,1.0,1e-11,"Coupled nuclear and chemical mole conservation");
    bool rejected=false;
    try{advanceDecayAndLinearRepartition(
        network,{{"X",NA}},10,phases,"aqueous",
        {{"ElementX",{{"aqueous",1},{"organic",3}}}});}
    catch(const std::invalid_argument&){rejected=true;}
    if(!rejected)throw std::runtime_error("Missing daughter partition data accepted");
    std::cout<<"Restricted nuclear-to-neutral-phase coupling checks passed\n";
}
