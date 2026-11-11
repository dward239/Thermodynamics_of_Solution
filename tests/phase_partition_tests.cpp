#include "phase_partition.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace thermo;
static void near(double a,double b,double tol,const char* message){
    if(!std::isfinite(a)||std::abs(a-b)>tol*std::max(1.0,std::abs(b)))
        throw std::runtime_error(message);
}
int main(){
    const std::vector<PartitionPhase> phases={{"water",1.0},{"organic_1",2.0},
                                              {"organic_2",0.5},{"headspace",3.0}};
    const std::vector<PartitionSolute> solutes={
        {"NeutralA",10.0,{{"water",1.0},{"organic_1",2.0},
                           {"organic_2",4.0},{"headspace",0.5}}},
        {"NeutralB",20.0,{{"water",1.0},{"organic_1",0.5},
                           {"organic_2",0.25},{"headspace",1.0}}}};
    const auto result=solveLinearPhasePartition(phases,"water",solutes);
    for(const auto& s:solutes){
        double total=0;
        for(const auto& p:phases)
            total+=result.moles.at(s.id).at(p.id);
        near(total,s.totalMoles,1e-12,"Partition mass conservation");
        const double cref=result.concentrationMolM3.at(s.id).at("water");
        for(const auto& p:phases)
            near(result.concentrationMolM3.at(s.id).at(p.id)/cref,
                 s.concentrationRatioToReference.at(p.id),1e-12,
                 "Partition ratio equilibrium");
    }
    near(result.moles.at("NeutralA").at("water"),
         10.0/(1.0+4.0+2.0+1.5),1e-12,"Closed-form distribution");
    bool failed=false;
    try{solveLinearPhasePartition(phases,"water",
        {{"Bad",1,{{"water",1.0},{"organic_1",2.0}}}});}
    catch(const std::invalid_argument&){failed=true;}
    if(!failed)throw std::runtime_error("Missing phase coefficient accepted");
    failed=false;
    try{solveLinearPhasePartition(phases,"water",
        {{"Bad",1,{{"water",1.0},{"organic_1",2.0},
                   {"organic_2",4.0},{"headspace",-1.0}}}});}
    catch(const std::invalid_argument&){failed=true;}
    if(!failed)throw std::runtime_error("Negative partition coefficient accepted");
    std::cout<<"Multi-liquid and gas ideal-linear partition tests passed\n";
}
