#include "aqueous_engine.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <algorithm>
using namespace aqueous;
void require(bool v,const char* message) {if(!v) throw std::runtime_error(message);}
void close(double a,double b,double tolerance,const char* message) {
    require(std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<=tolerance,message);
}
Model association() {
    Model m;
    m.components={{"A+",1,true,0.01},{"B-",-1,true,0.01}};
    Complex ab;ab.name="AB";ab.charge=0;ab.logKReference=2.0;
    ab.basisStoichiometry={{"A+",1},{"B-",1}};
    m.complexes={ab};
    return m;
}
int main() {
    auto ideal=association();
    const auto r=solve(ideal);
    const double expected=2*0.01/(1+std::sqrt(1+4*100*0.01));
    close(r.concentrationMolL.at("A+"),expected,1e-10,"Analytical association");
    close(r.chargeResidualMolL,0,1e-12,"Neutrality diagnostic");
    require(r.normalizedResidual<1e-9,"Mass residual");
    require(toJson(r).find("\"activities\"")!=std::string::npos,"JSON generated");

    // Temperature dependence of the dimensionless formation constant:
    // d(log10K)/dT = DeltaH/(R ln10 T^2), constant DeltaH assumed.
    Complex reaction;reaction.name="AB";reaction.logKReference=2.0;reaction.deltaHkJMol=20.0;
    reaction.charge=0;reaction.basisStoichiometry={{"A+",1},{"B-",1}};
    close(logKAtTemperature(reaction,298.15,298.15),2.0,1e-13,"Reference K");
    require(logKAtTemperature(reaction,313.15,298.15)>2.0,"Endothermic formation K increases with T");

    ideal.activity=ActivityModel::Davies;
    const auto nonideal=solve(ideal);
    require(nonideal.ionicStrengthMolL>0 && nonideal.ionicStrengthMolL<0.1,"Davies I range");
    require(nonideal.activity.at("A+")<nonideal.concentrationMolL.at("A+"),"Davies gamma<1");
    close(nonideal.chargeResidualMolL,0,1e-11,"Davies charge balance diagnostic");

    // Carbonate with free H+ explicitly imposes electroneutrality.
    auto carbonate=readModel("examples/carbonate_ideal.txt");
    const auto water=solve(carbonate);
    require(water.hasPH && water.pH>4 && water.pH<12,"Calculated carbonate pH plausible");
    require(std::abs(water.chargeResidualMolL)<1e-10,"Carbonate charge balance");
    require(water.normalizedResidual<1e-9,"Carbonate mass balance");
    const double car=water.concentrationMolL.at("CO2*")
        +water.concentrationMolL.at("HCO3-")+water.concentrationMolL.at("CO3--");
    close(car,0.001,1e-10,"Carbon conservation");
    close(water.activity.at("H+")*water.activity.at("OH-"),1e-14,1e-19,"Water mass action");
    auto nonidealCarbonate=readModel("examples/carbonate_davies.txt");
    auto other=solve(nonidealCarbonate);
    require(other.hasPH && std::abs(other.chargeResidualMolL)<1e-10,"Davies carbonate electroneutral");
    require(other.ionicStrengthMolL>0 && other.ionicStrengthMolL<0.1,"Davies carbonate I");
    require(std::abs(other.pH-water.pH)>1e-5,"Nonideal correction changes pH");

    carbonate.temperatureK=313.15;
    const auto warm=solve(carbonate);
    require(std::abs(warm.pH-water.pH)>1e-4,"Temperature changes speciation/pH");

    bool caught=false;
    try { auto bad=carbonate;bad.temperatureK=350;solve(bad); }
    catch(const std::invalid_argument&) {caught=true;}
    require(caught,"Temperature limit is enforced");

    caught=false;
    try {auto bad=ideal;bad.components[0].totalMolL=2;bad.components[1].totalMolL=2;solve(bad);}
    catch(const std::exception&) {caught=true;}
    require(caught,"High ionic strength rejected");

    caught=false;
    try {auto bad=association();bad.complexes[0].charge=1;solve(bad);}
    catch(const std::invalid_argument&) {caught=true;}
    require(caught,"Charge-inconsistent reaction rejected");
    std::cout<<"MVP+1 numerical checks passed\n";
}
