#include "ideal_multiphase_equilibrium.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace thermo;
static void closeTo(double actual,double expected,double tol,const char* msg){
    if(!std::isfinite(actual)||
       std::abs(actual-expected)>tol*std::max(1.0,std::abs(expected)))
        throw std::runtime_error(msg);
}
static void shouldThrow(const IdealMultiphaseModel& m,const char* msg){
    bool failed=false;
    try{(void)solveIdealMultiphase(m);}
    catch(const std::invalid_argument&){failed=true;}
    catch(const std::runtime_error&){failed=true;}
    if(!failed)throw std::runtime_error(msg);
}
static void testTwoLiquidPhases(){
    IdealMultiphaseModel m;
    m.components={{"X",0.01}};
    m.phases={{"water",0.001,1000.0},{"organic",0.002,1000.0,IdealFluidPhaseKind::OrganicLiquid}};
    m.species={{"X(aq)","water",{{"X",1}},0.0},
               {"X(org)","organic",{{"X",1}},-std::log(2.0)}};
    auto r=solveIdealMultiphase(m);
    closeTo(r.speciesMoles.at("X(aq)"),0.002,1e-11,"Aqueous solute allocation");
    closeTo(r.speciesMoles.at("X(org)"),0.008,1e-11,"Organic solute allocation");
    closeTo(r.concentrationsMolM3.at("X(org)")/
            r.concentrationsMolM3.at("X(aq)"),2,1e-11,
            "Chemical-potential partition coefficient");
    closeTo(r.componentResidualMol.at("X"),0,1e-11,"X material balance");
    if(r.maxNormalizedResidual>1e-10)
        throw std::runtime_error("Unexpected phase solver residual");
}
static void testFormationEquilibrium(){
    IdealMultiphaseModel m;
    m.components={{"A",0.01},{"B",0.01}};
    m.phases={{"water",0.001,1000.0}};
    m.species={{"A(aq)","water",{{"A",1}},0},
               {"B(aq)","water",{{"B",1}},0},
               {"AB(aq)","water",{{"A",1},{"B",1}},-std::log(100.0)}};
    auto r=solveIdealMultiphase(m);
    const double expectedFree=(std::sqrt(5.0)-1.0)/200.0;
    closeTo(r.speciesMoles.at("A(aq)"),expectedFree,1e-11,"Closed form free A");
    closeTo(r.speciesMoles.at("B(aq)"),expectedFree,1e-11,"Closed form free B");
    closeTo(r.speciesMoles.at("AB(aq)"),0.01-expectedFree,1e-11,
            "Closed form complex AB");
    closeTo(r.speciesMoles.at("AB(aq)")/
            (r.speciesMoles.at("A(aq)")*r.speciesMoles.at("B(aq)")),
            100.0,1e-10,"Mass action constant at 1 mol/L reference");
}
static void testThreePhaseDistribution(){
    IdealMultiphaseModel m;
    m.components={{"GaslikeNeutral",2}};
    m.phases={{"aqueous",1,1},{"organic",1,1,IdealFluidPhaseKind::OrganicLiquid},{"headspace",1,1,IdealFluidPhaseKind::Gas}};
    m.species={{"S(aq)","aqueous",{{"GaslikeNeutral",1}},0},
               {"S(org)","organic",{{"GaslikeNeutral",1}},-std::log(2.0)},
               {"S(gas)","headspace",{{"GaslikeNeutral",1}},-std::log(3.0)}};
    auto r=solveIdealMultiphase(m);
    closeTo(r.speciesMoles.at("S(aq)"),2.0/6,1e-11,"Aqueous fraction");
    closeTo(r.speciesMoles.at("S(org)"),4.0/6,1e-11,"Organic fraction");
    closeTo(r.speciesMoles.at("S(gas)"),1.0,1e-11,"Gas fraction");
    closeTo(r.componentResidualMol.at("GaslikeNeutral"),0,1e-11,
            "Three phase mass balance");
}
static void testInvalidModels(){
    IdealMultiphaseModel m;
    m.components={{"A",0.01},{"B",0.01}};
    m.phases={{"water",0.001,1000.0}};
    m.species={{"AB","water",{{"A",1},{"B",1}},0}};
    shouldThrow(m,"Rank-deficient stoichiometric Jacobian should reject");
    m.species.push_back({"A","water",{{"A",1}},0});
    m.species.push_back({"B","water",{{"B",1}},0});
    auto r=solveIdealMultiphase(m);
    closeTo(r.componentResidualMol.at("A"),0,1e-11,"Valid repaired model");
    auto invalid=m;
    invalid.species[0].phaseId="nowhere";
    shouldThrow(invalid,"Unknown phase should reject");
    invalid=m;invalid.species[0].componentStoichiometry["Unknown"]=1;
    shouldThrow(invalid,"Unknown conserved component should reject");
    invalid=m;invalid.components[0].totalMol=0;
    shouldThrow(invalid,"Zero total should reject");
}
int main(){
    testTwoLiquidPhases();
    testFormationEquilibrium();
    testThreePhaseDistribution();
    testInvalidModels();
    std::cout<<"General ideal-dilute multi-fluid mass action tests passed\n";
}
