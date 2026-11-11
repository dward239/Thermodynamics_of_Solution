#include "aqueous_sweep.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
using namespace aqueous;
static void expect(bool value,const std::string& message) {
    if(!value) throw std::runtime_error(message);
}
int main() {
    Model m;
    m.components={{"A+",1,true,0.01},{"B-",-1,true,0.01}};
    Complex ab;ab.name="AB";ab.charge=0;ab.logKReference=2;
    ab.basisStoichiometry={{"A+",1},{"B-",1}};
    m.complexes.push_back(ab);
    auto t=runSweep(m,"temperature","",298.15,313.15,4);
    expect(t.size()==4,"Sweep count");
    expect(std::abs(t[0].parameterValue-298.15)<1e-12,"Initial T");
    expect(std::abs(t[3].parameterValue-313.15)<1e-12,"Final T");
    expect(t[1].result.normalizedResidual<1e-9,"Sweep residual");
    expect(sweepCsv(t,"temperature").find("concentration_A+")!=std::string::npos,"Sweep CSV species");
    auto c=runSweep(m,"total","A+",0.001,0.01,3);
    expect(c[2].result.concentrationMolL.at("A+")>c[0].result.concentrationMolL.at("A+"),
           "Higher total gives higher free concentration");
    bool failed=false;
    try {runSweep(m,"temperature","",250,313,4);}
    catch(const std::invalid_argument&) {failed=true;}
    expect(failed,"Reject unsafe temperature sweep");
    failed=false;
    try {runSweep(m,"total","does_not_exist",0.001,0.01,2);}
    catch(const std::invalid_argument&) {failed=true;}
    expect(failed,"Reject unknown component");
    failed=false;
    try {runSweep(m,"total","A+",0,0.01,3);}
    catch(const std::invalid_argument&) {failed=true;}
    expect(failed,"Reject zero total");
    failed=false;
    try {runSweep(m,"temperature","",298.15,313.15,1);}
    catch(const std::invalid_argument&) {failed=true;}
    expect(failed,"Reject inadequate sweep count");
    std::cout<<"Sweep tests passed\n";
}
