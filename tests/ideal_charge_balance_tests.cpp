#include "ideal_multiphase_equilibrium.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace thermo;
static void require(bool value,const char* message){
    if(!value)throw std::runtime_error(message);
}
static void near(double got,double want,double tolerance,const char* message){
    if(!std::isfinite(got)||std::abs(got-want)>tolerance*std::max(1.0,std::abs(want)))
        throw std::runtime_error(message);
}
template<class F> static void expectInvalid(F f,const char* message){
    bool threw=false;
    try{f();}catch(const std::invalid_argument&){threw=true;}
    require(threw,message);
}
static IdealMultiphaseModel electrolyte(){
    IdealMultiphaseModel m;
    m.components={{"Na",0.01,true},{"Cl",0.0,false}};
    m.phases={{"water",0.001,1000.0}};
    m.species={
        {"Na+","water",{{"Na",1}},0.0,+1},
        {"Cl-","water",{{"Cl",1}},0.0,-1}
    };
    m.electroneutralPhaseId="water";
    return m;
}
static void testSimpleCounterion(){
    const auto r=solveIdealMultiphase(electrolyte());
    near(r.speciesMoles.at("Na+"),0.01,1e-11,"Na mass balance");
    near(r.speciesMoles.at("Cl-"),0.01,1e-11,"Electroneutral free counterion");
    near(r.calculatedComponentMol.at("Cl"),0.01,1e-11,"Inferred free Cl total");
    near(r.componentResidualMol.at("Na"),0,1e-11,"Fixed Na residual");
    require(r.componentResidualMol.count("Cl")==0,"Free total is not a fixed residual");
    near(r.chargeResidualMol,0,1e-11,"Charge balance residual");
}
static void testIonicAssociation(){
    auto m=electrolyte();
    m.species.push_back({"NaCl(aq)","water",{{"Na",1},{"Cl",1}},
                         -std::log(100.0),0});
    const auto r=solveIdealMultiphase(m);
    const double free=(std::sqrt(5.0)-1.0)/200.0;
    near(r.speciesMoles.at("Na+"),free,1e-10,"Free Na at equilibrium");
    near(r.speciesMoles.at("Cl-"),free,1e-10,"Free Cl from charge balance");
    near(r.speciesMoles.at("NaCl(aq)"),0.01-free,1e-10,"Neutral ion pair");
    near(r.speciesMoles.at("NaCl(aq)")/
         (r.speciesMoles.at("Na+")*r.speciesMoles.at("Cl-")),
         100,1e-9,"Ideal association mass action");
    near(r.chargeResidualMol,0,1e-11,"Charge balance with ion pairing");
    near(r.calculatedComponentMol.at("Cl"),0.01,1e-10,"Inferred Cl inventory");
}
static void testNeutralOtherPhases(){
    auto m=electrolyte();
    m.components.push_back({"X",0.004,true});
    m.phases.push_back({"organic",0.001,1000.0});
    m.phases.push_back({"headspace",0.001,1000.0});
    m.species.push_back({"X(aq)","water",{{"X",1}},0,0});
    m.species.push_back({"X(org)","organic",{{"X",1}},-std::log(2.0),0});
    m.species.push_back({"X(g)","headspace",{{"X",1}},-std::log(1.0),0});
    auto r=solveIdealMultiphase(m);
    near(r.speciesMoles.at("X(aq)"),0.001,1e-11,"Neutral solute aqueous");
    near(r.speciesMoles.at("X(org)"),0.002,1e-11,"Neutral solute organic");
    near(r.speciesMoles.at("X(g)"),0.001,1e-11,"Neutral solute gas");
    near(r.chargeResidualMol,0,1e-11,"Charged water remains electroneutral");
}
static void testInvalidChargeModels(){
    auto m=electrolyte();
    auto invalid=m;
    invalid.electroneutralPhaseId.clear();
    expectInvalid([&]{(void)solveIdealMultiphase(invalid);},
                  "Charged species without phase constraint must reject");
    invalid=m;
    invalid.components[1].fixedTotal=true;
    invalid.components[1].totalMol=0.01;
    expectInvalid([&]{(void)solveIdealMultiphase(invalid);},
                  "Fully specified inventory cannot also impose independent charge equation");
    invalid=m;
    invalid.components[0].fixedTotal=false;
    invalid.components[0].totalMol=0;
    expectInvalid([&]{(void)solveIdealMultiphase(invalid);},
                  "More than one free component must reject");
    invalid=m;
    invalid.components[1].totalMol=0.001;
    expectInvalid([&]{(void)solveIdealMultiphase(invalid);},
                  "Free component cannot also specify a total");
    invalid=m;
    invalid.species[1].charge=0;
    invalid.species[0].charge=0;
    expectInvalid([&]{(void)solveIdealMultiphase(invalid);},
                  "Unnecessary charge constraint must reject");
    invalid=m;
    invalid.phases.push_back({"organic",0.001,1000.0});
    invalid.species[1].phaseId="organic";
    expectInvalid([&]{(void)solveIdealMultiphase(invalid);},
                  "Charged interphase transfer unsupported");
    invalid=m;
    invalid.electroneutralPhaseId="missing";
    expectInvalid([&]{(void)solveIdealMultiphase(invalid);},
                  "Unknown electroneutral phase must reject");
    expectInvalid([&]{(void)solveIdealMultiphaseWithSolids(
        m,{{"NaCl(s)",{{"Na",1},{"Cl",1}},-10.0}});},
        "Charged solution with active solids unsupported");
}
int main(){
    testSimpleCounterion();
    testIonicAssociation();
    testNeutralOtherPhases();
    testInvalidChargeModels();
    std::cout<<"Charge-balanced ideal-dilute multiphase cases passed\n";
}
