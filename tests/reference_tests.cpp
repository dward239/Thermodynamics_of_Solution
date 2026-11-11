#include "aqueous_engine.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace aqueous;
static void close(double actual,double expected,double tolerance,const char* msg) {
    if(!std::isfinite(actual) || std::abs(actual-expected)>tolerance)
        throw std::runtime_error(msg);
}
int main() {
    // Independent ideal carbonate reference derived from a closed-form species
    // distribution at fixed total inorganic carbon, plus scalar charge bisection.
    // This intentionally does NOT reuse the multivariable Newton implementation.
    const double carbonTotal=1e-3;
    const double sodiumTotal=1e-3;
    const double K1=std::pow(10.,-6.35);
    const double K2=std::pow(10.,-10.33);
    const double Kw=1e-14;
    auto charge=[&](double H){
        const double denom=H*H+K1*H+K1*K2;
        const double bicarbonate=carbonTotal*K1*H/denom;
        const double carbonate=carbonTotal*K1*K2/denom;
        return H+sodiumTotal-bicarbonate-2*carbonate-Kw/H;
    };
    double low=1e-14, high=1e-1;
    if(!(charge(low)<0 && charge(high)>0))
        throw std::runtime_error("Analytic charge bracket was not valid");
    for(int i=0;i<180;++i) {
        double mid=std::sqrt(low*high);
        if(charge(mid)>0) high=mid; else low=mid;
    }
    const double H=std::sqrt(low*high);
    const double denom=H*H+K1*H+K1*K2;
    const double CO2=carbonTotal*H*H/denom;
    const double HCO3=carbonTotal*K1*H/denom;
    const double CO3=carbonTotal*K1*K2/denom;
    const auto result=solve(readModel("examples/carbonate_ideal.txt"));
    close(result.pH,-std::log10(H),2e-7,"Independent carbonate reference pH");
    close(result.concentrationMolL.at("CO2*"),CO2,1e-9,"Independent CO2");
    close(result.concentrationMolL.at("HCO3-"),HCO3,1e-9,"Independent bicarbonate");
    close(result.concentrationMolL.at("CO3--"),CO3,1e-9,"Independent carbonate");
    close(result.concentrationMolL.at("OH-"),Kw/H,1e-9,"Independent hydroxide");
    close(result.chargeResidualMolL,0,1e-9,"Electroneutrality");
    const double Kref=10.0;
    Complex rxn;rxn.logKReference=std::log10(Kref);rxn.deltaHkJMol=20.0;
    const double T=313.15,Tref=298.15,R=8.31446261815324;
    const double expectedLogK=std::log10(Kref)+rxn.deltaHkJMol*1000.0/
        (R*std::log(10.0))*(1.0/Tref-1.0/T);
    close(logKAtTemperature(rxn,T,Tref),expectedLogK,1e-12,
          "Independent van't Hoff calculation");
    std::cout<<"Independent analytic ideal carbonate reference passed\n";
}
