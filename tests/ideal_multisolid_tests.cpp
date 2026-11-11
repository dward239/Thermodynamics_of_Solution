#include "ideal_multiphase_equilibrium.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace thermo;
static void closeTo(double a,double b,double tol,const char* msg){
    if(!std::isfinite(a)||
       std::abs(a-b)>tol*std::max(1.0,std::abs(b)))
        throw std::runtime_error(msg);
}
static IdealMultiphaseModel fluid(const std::vector<IdealComponent>& totals){
    IdealMultiphaseModel m;
    m.components=totals;
    m.phases={{"water",0.001,1000.0}};
    for(const auto& c:totals)
        m.species.push_back({c.id+"(aq)","water",{{c.id,1}},0.0});
    return m;
}
static IdealPureSolid solid(const char* name,const std::map<std::string,int>& nu,
                            double ksp){
    return {name,nu,std::log(ksp)};
}
static void testSingleSolid(){
    auto m=fluid({{"A",0.01},{"B",0.01}});
    auto r=solveIdealMultiphaseWithSolids(
        m,{solid("AB(s)",{{"A",1},{"B",1}},1e-5)});
    const double free=std::sqrt(1e-5);
    closeTo(r.fluidSolution.speciesMoles.at("A(aq)"),free,1e-9,"Saturated A");
    closeTo(r.fluidSolution.speciesMoles.at("B(aq)"),free,1e-9,"Saturated B");
    closeTo(r.solidMoles.at("AB(s)"),0.01-free,1e-9,"Precipitated AB");
    closeTo(r.solidAffinityOverRT.at("AB(s)"),0,1e-9,"AB affinity");
    closeTo(r.fluidSolution.componentResidualMol.at("A"),0,1e-9,"Full A conservation");
}
static void testNoSolid(){
    auto m=fluid({{"A",0.01},{"B",0.01}});
    auto r=solveIdealMultiphaseWithSolids(
        m,{solid("AB(s)",{{"A",1},{"B",1}},1e-2)});
    closeTo(r.solidMoles.at("AB(s)"),0,1e-12,"Solid absent when undersaturated");
    if(r.solidAffinityOverRT.at("AB(s)")>=0)
        throw std::runtime_error("Expected negative absent-solid affinity");
}
static void testCompetingPolymorphs(){
    auto m=fluid({{"A",0.01},{"B",0.01}});
    auto r=solveIdealMultiphaseWithSolids(m,{
        solid("LessStable",{{"A",1},{"B",1}},1e-5),
        solid("MoreStable",{{"A",1},{"B",1}},1e-6)});
    closeTo(r.fluidSolution.speciesMoles.at("A(aq)"),0.001,1e-9,
            "Stable polymorph saturation");
    closeTo(r.solidMoles.at("LessStable"),0,1e-9,"Unstable polymorph absent");
    closeTo(r.solidMoles.at("MoreStable"),0.009,1e-9,"Stable polymorph present");
    if(r.solidAffinityOverRT.at("LessStable")>=0)
        throw std::runtime_error("Competing metastable phase not undersaturated");
}
static void testTwoSimultaneousSolids(){
    auto m=fluid({{"A",0.02},{"B",0.01},{"C",0.01}});
    auto r=solveIdealMultiphaseWithSolids(m,{
        solid("AB(s)",{{"A",1},{"B",1}},1e-6),
        solid("AC(s)",{{"A",1},{"C",1}},1e-6)});
    const double a=std::sqrt(2e-6);
    const double b=1e-6/a;
    closeTo(r.fluidSolution.speciesMoles.at("A(aq)"),a,1e-8,"Shared A");
    closeTo(r.fluidSolution.speciesMoles.at("B(aq)"),b,1e-8,"Shared B");
    closeTo(r.fluidSolution.speciesMoles.at("C(aq)"),b,1e-8,"Shared C");
    closeTo(r.solidMoles.at("AB(s)"),0.01-b,1e-8,"AB solid inventory");
    closeTo(r.solidMoles.at("AC(s)"),0.01-b,1e-8,"AC solid inventory");
    closeTo(r.fluidSolution.componentResidualMol.at("A"),0,1e-8,"A conservation");
    closeTo(r.fluidSolution.componentResidualMol.at("B"),0,1e-8,"B conservation");
    closeTo(r.fluidSolution.componentResidualMol.at("C"),0,1e-8,"C conservation");
}
static void testInvalid(){
    auto m=fluid({{"A",0.01},{"B",0.01}});
    bool invalid=false;
    try{solveIdealMultiphaseWithSolids(m,{
        solid("A(aq)",{{"A",1},{"B",1}},1e-5)});}
    catch(const std::invalid_argument&){invalid=true;}
    if(!invalid)throw std::runtime_error("Did not reject duplicate phase/species id");
    invalid=false;
    try{solveIdealMultiphaseWithSolids(m,{
        solid("Bad",{{"Unknown",1}},1e-5)});}
    catch(const std::invalid_argument&){invalid=true;}
    if(!invalid)throw std::runtime_error("Did not reject unknown solid component");
}
int main(){
    testSingleSolid();
    testNoSolid();
    testCompetingPolymorphs();
    testTwoSimultaneousSolids();
    testInvalid();
    std::cout<<"Ideal pure-solid active phase equilibrium tests passed\n";
}
