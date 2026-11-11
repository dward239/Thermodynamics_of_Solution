#include "multiphase_system.h"
#include "decay_network.h"
#include "lumped_thermal.h"
#include "radiolysis_sources.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace thermo;
static void expect(bool b,const char* msg){if(!b)throw std::runtime_error(msg);}
static void closeTo(double v,double x,double tol,const char* msg){
    if(!std::isfinite(v)||std::abs(v-x)>tol*std::max(1.0,std::abs(x)))
        throw std::runtime_error(msg);
}
template<class F> static void expectInvalid(F f,const char* msg){
    bool threw=false;try{f();}catch(const std::invalid_argument&){threw=true;}
    expect(threw,msg);
}
static void testPhaseStructure(){
    MultiphaseSystem m;
    m.phases={{"aq",PhaseKind::AqueousLiquid,0.001,298.15},
              {"org",PhaseKind::OrganicLiquid,0.001,298.15},
              {"gas",PhaseKind::Gas,0.0005,298.15}};
    m.species={{"Na+","aq",1,{{"Na",1}}},
               {"Cl-","aq",-1,{{"Cl",1}}},
               {"NaCl(aq)","aq",0,{{"Na",1},{"Cl",1}}},
               {"NaCl(org)","org",0,{{"Na",1},{"Cl",1}}}};
    m.reactions={{"association",{{"Na+",-1},{"Cl-",-1},{"NaCl(aq)",1}},0.0},
                 {"partition",{{"NaCl(aq)",-1},{"NaCl(org)",1}},0.0}};
    m.validate();
    auto bad=m;bad.reactions[0].coefficients["Cl-"]=1;
    expectInvalid([&]{bad.validate();},"Reject invalid reaction");
    bad=m;bad.species[0].phaseId="missing";
    expectInvalid([&]{bad.validate();},"Reject unknown phase");
    bad=m;bad.phases[0].volumeM3=0;
    expectInvalid([&]{bad.validate();},"Reject zero phase volume");
    bad=m;bad.species[0].id=bad.species[1].id;
    expectInvalid([&]{bad.validate();},"Reject duplicate species");
    bad=m;bad.reactions[1].coefficients["NaCl(org)"]=2;
    expectInvalid([&]{bad.validate();},"Reject nonconserving partition reaction");
}
static void testDecay(){
    const DecayNetwork chain({{"X","X",10.0,2.0,{{"Y",1.0}},{{"aq",0.25},{"gas",0.1}}},
                              {"Y","Y",0,0,{},{}}});
    const auto x=chain.advance({{"X",8.0}},10.0);
    closeTo(x.finalAtoms.at("X"),4,1e-12,"Parent half life");
    closeTo(x.finalAtoms.at("Y"),4,1e-12,"Daughter ingrowth");
    closeTo(x.decaysByNuclide.at("X"),4,1e-12,"Cumulative decay events");
    closeTo(x.emittedEnergyJ,8,1e-12,"Emitted energy");
    closeTo(x.depositedEnergyJByPhase.at("aq"),2,1e-12,"Absorbed aqueous energy");
    closeTo(x.depositedEnergyJByPhase.at("gas"),0.8,1e-12,"Absorbed gas energy");
    auto zero=chain.advance({{"X",8.0}},0.0);
    closeTo(zero.finalAtoms.at("X"),8,1e-12,"Zero time");
    closeTo(zero.emittedEnergyJ,0,1e-12,"Zero time energy");
    expectInvalid([&]{chain.advance({{"X",-1}},1);},"Reject negative atoms");
    expectInvalid([&]{chain.advance({{"NoSuchNuclide",1}},1);},"Reject unknown nuclide");
    expectInvalid([&]{chain.advance({{"X",1}},1e9);},"Reject numerically stiff step");
    const DecayNetwork two({{"A","E1",10,1,{{"B",1}},{}},
                            {"B","E2",20,1,{{"C",1}},{}},
                            {"C","E3",0,0,{},{}}});
    auto a=two.advance({{"A",10.0}},10.0);
    const double b=20.0*(std::sqrt(0.5)-0.5);
    closeTo(a.finalAtoms.at("A"),5,1e-11,"Bateman A");
    closeTo(a.finalAtoms.at("B"),b,1e-11,"Bateman B");
    closeTo(a.finalAtoms.at("C"),5-b,1e-11,"Bateman C");
    closeTo(a.decaysByNuclide.at("B"),5-b,1e-11,"Cumulative secondary decay");
    const DecayNetwork branching({{"A","E1",10,0,{{"B",0.3},{"C",0.7}},{}},
                                  {"B","E2",0,0,{},{}},
                                  {"C","E3",0,0,{},{}}});
    auto branches=branching.advance({{"A",10}},10);
    closeTo(branches.finalAtoms.at("B"),1.5,1e-11,"30 percent branch");
    closeTo(branches.finalAtoms.at("C"),3.5,1e-11,"70 percent branch");
    expectInvalid([]{DecayNetwork fail({{"A","E",10,1,{{"B",0.5}},{}},
                                       {"B","E",0,0,{},{}}});},"Reject incomplete branches");
    expectInvalid([]{DecayNetwork fail({{"A","E",10,1,{{"B",1}},{{"aq",1.2}}},
                                       {"B","E",0,0,{},{}}});},"Reject over-unity deposition");
}
static void testRadiolysis(){
    const double energy=100.0*1.602176634e-19*6.02214076e23;
    const auto products=grossRadiolyticProductionMol({{"aq",energy}},
                                                     {{"aq","H2",1.0},{"aq","O2",0.5}});
    closeTo(products.at("aq::H2"),1,1e-12,"G value conversion");
    closeTo(products.at("aq::O2"),0.5,1e-12,"Fractional G value");
    expectInvalid([]{grossRadiolyticProductionMol({{"aq",1}},{{"gas","H2",1}});},
                  "Reject missing deposition");
    expectInvalid([]{grossRadiolyticProductionMol({{"aq",1}},{{"aq","H2",-1}});},
                  "Reject invalid G");
}
static void testThermal(){
    LumpedThermalProperties noLoss{2,500,0,300};
    closeTo(advanceLumpedTemperatureK(300,10,1000,1,noLoss),301,1e-12,"Adiabatic");
    closeTo(advanceLumpedTemperatureK(300,10,1000,0.5,noLoss),300.5,1e-12,"Heat fraction");
    LumpedThermalProperties cooling{1,100,10,300};
    closeTo(advanceLumpedTemperatureK(350,10,0,1,cooling),
            300+50*std::exp(-1.0),1e-12,"Newton cooling");
    expectInvalid([&]{advanceLumpedTemperatureK(300,0,100,1,noLoss);},"Zero dt with energy");
    expectInvalid([&]{advanceLumpedTemperatureK(300,10,100,1.1,noLoss);},
                  "Reject over-unity heat fraction");
}
int main(){
    testPhaseStructure();testDecay();testRadiolysis();testThermal();
    std::cout<<"Radio-thermochemical foundation tests passed\n";
}
