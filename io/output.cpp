#include "io/output.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace thermo {
namespace {
void writeNumber(std::ostream& out,double value) {
    if(std::isfinite(value)) out<<value;
    else out<<"null";
}
}
void saveJson(const Result& r,const std::string& filename) {
    std::ofstream out(filename);
    if(!out) throw std::runtime_error("Cannot write JSON file: "+filename);
    out<<std::setprecision(17);
    out<<"{\n  \"temperature_C\": ";
    writeNumber(out,r.temperatureC);
    out<<",\n  \"ionic_strength_molal\": ";
    writeNumber(out,r.ionicStrength);
    out<<",\n  \"pH\": ";
    writeNumber(out,r.pH);
    out<<",\n  \"charge_residual_molal\": ";
    writeNumber(out,r.netCharge);
    out<<",\n  \"max_scaled_residual\": ";
    writeNumber(out,r.maxScaledResidual);
    out<<",\n  \"iterations\": "<<r.iterations<<",\n  \"species\": [\n";
    bool first=true;
    for(const auto& [name,m]:r.molalities) {
        if(!first) out<<",\n";
        first=false;
        // Species names are validated by input grammar in published examples.
        out<<"    {\"name\": \""<<name<<"\", \"molality\": ";
        writeNumber(out,m);
        out<<", \"activity\": ";
        writeNumber(out,r.activities.at(name));
        out<<", \"gamma\": ";
        writeNumber(out,r.gammas.at(name));
        out<<"}";
    }
    out<<"\n  ],\n  \"component_residuals_molal\": {";
    first=true;
    for(const auto& [name,v]:r.componentResiduals) {
        if(!first) out<<", ";
        first=false;
        out<<"\""<<name<<"\": ";
        writeNumber(out,v);
    }
    out<<"}\n}\n";
}
void saveCsv(const Result& r,const std::string& filename) {
    std::ofstream out(filename);
    if(!out) throw std::runtime_error("Cannot write CSV file: "+filename);
    out<<std::setprecision(17)<<"species,molality_mol_per_kg,activity,gamma\n";
    for(const auto& [name,m]:r.molalities)
        out<<name<<","<<m<<","<<r.activities.at(name)<<","<<r.gammas.at(name)<<"\n";
}
} // namespace thermo
