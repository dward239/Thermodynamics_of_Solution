#include "aqueous_engine.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace aqueous;
static void expect(bool ok, const char* msg) { if(!ok) throw std::runtime_error(msg); }
static void close(double actual,double expected,double tol,const char* msg) {
    expect(std::isfinite(actual)&&std::isfinite(expected)&&std::abs(actual-expected)<=tol,msg);
}
static Model model(double total,double logKsp) {
    Model m;
    m.components={{"A+",1,true,total},{"B-",-1,true,total}};
    Solid s;s.name="AB_s";s.logKspReference=logKsp;
    s.dissolutionStoichiometry={{"A+",1},{"B-",1}};
    m.solids.push_back(s);
    return m;
}
int main() {
    auto m=model(0.01,-5.0);
    auto result=solve(m);
    const double expectedFree=std::sqrt(1e-5);
    const double expectedPrecipitate=0.01-expectedFree;
    close(result.concentrationMolL.at("A+"),expectedFree,1e-9,"A at saturation");
    close(result.concentrationMolL.at("B-"),expectedFree,1e-9,"B at saturation");
    close(result.solidPrecipitatedMolL.at("AB_s"),expectedPrecipitate,1e-9,"Solid amount");
    close(result.saturationIndex.at("AB_s"),0,1e-7,"Saturation SI");
    close(result.concentrationMolL.at("A+")+result.solidPrecipitatedMolL.at("AB_s"),
          0.01,1e-9,"A atom conservation");
    close(result.chargeResidualMolL,0,1e-11,"Charge conservation");
    expect(toJson(result).find("\"solids_precipitated_mol_L\"")!=std::string::npos,
           "JSON contains precipitated solid");

    auto dilute=solve(model(0.001,-5));
    close(dilute.solidPrecipitatedMolL.at("AB_s"),0,1e-13,"No solid below Ksp");
    close(dilute.saturationIndex.at("AB_s"),-1.0,1e-7,"Below saturation SI");

    bool failed=false;
    try {auto bad=m;bad.activity=ActivityModel::Davies;solve(bad);}
    catch(const std::invalid_argument&) {failed=true;}
    expect(failed,"No nonideal precipitation support");
    failed=false;
    try {auto bad=m;bad.solids.push_back(bad.solids.front());solve(bad);}
    catch(const std::invalid_argument&) {failed=true;}
    expect(failed,"No competing solids support");
    failed=false;
    try {auto bad=m;bad.solids[0].dissolutionStoichiometry["B-"]=2;solve(bad);}
    catch(const std::invalid_argument&) {failed=true;}
    expect(failed,"Reject charged dissolution stoichiometry");
    failed=false;
    try {auto bad=m;bad.components[0].fixedTotal=false;solve(bad);}
    catch(const std::invalid_argument&) {failed=true;}
    expect(failed,"Reject free charge-balance component with solid");

    // Verify parser supports the restricted solid syntax.
    const auto parsed=readModel("examples/precipitation_ideal.txt");
    close(solve(parsed).solidPrecipitatedMolL.at("AB_s"),expectedPrecipitate,1e-9,
          "Parsed precipitation model");
    std::cout<<"Restricted ideal precipitation tests passed\n";
}
