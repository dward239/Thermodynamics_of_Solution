#include "chemistry/model.h"
#include <cmath>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace thermo {
namespace {
void endOfLine(std::istringstream& in, int row) {
    std::string extra;
    if (in >> extra && extra[0] != '#')
        throw std::runtime_error("Unexpected text at line " + std::to_string(row));
}
}
double temperatureLogK(double logK25, double deltaHkJmol, double temperatureC) {
    if (!std::isfinite(logK25) || !std::isfinite(deltaHkJmol) ||
        !std::isfinite(temperatureC) || temperatureC < 0.0 || temperatureC > 80.0)
        throw std::invalid_argument("Invalid equilibrium data or temperature (valid: 0-80 C)");
    constexpr double R = 8.31446261815324; // J/(mol K)
    constexpr double Tref = 298.15;
    const double Tk = 273.15 + temperatureC;
    return logK25 - (1000.0 * deltaHkJmol / (R * std::log(10.0))) *
                            (1.0 / Tk - 1.0 / Tref);
}
Prepared prepare(Model model) {
    if (!std::isfinite(model.temperatureC) ||
        model.temperatureC < 0.0 || model.temperatureC > 80.0)
        throw std::invalid_argument("Temperature outside 0-80 Celsius model range");
    if (model.basis.empty()) throw std::invalid_argument("No primary basis species");
    std::map<std::string, std::size_t> index;
    Prepared out;
    out.model = std::move(model);
    const auto n = out.model.basis.size();
    for (std::size_t i = 0; i < n; ++i) {
        const auto& b = out.model.basis[i];
        if (b.name.empty() || !index.emplace(b.name,i).second)
            throw std::invalid_argument("Duplicate or empty primary species name: " + b.name);
        if ((!b.adjustable && (!std::isfinite(b.totalMolal) || b.totalMolal <= 0.0)) ||
            (b.adjustable && b.totalMolal != 0.0))
            throw std::invalid_argument("Primary totals must be positive; adjustable total must be omitted: " + b.name);
        Species s{b.name,b.charge,std::vector<int>(n,0),0.0};
        s.powers[i]=1;
        out.species.push_back(std::move(s));
        if (b.name == "H+") out.hydrogenIndex = static_cast<int>(i);
    }
    int adjustableCount=0;
    for (std::size_t i=0; i<n; ++i) if (out.model.basis[i].adjustable) {
        ++adjustableCount; out.adjustedIndex=static_cast<int>(i);
        if (out.model.adjustedBasis != out.model.basis[i].name)
            throw std::invalid_argument("Closure basis does not match adjustable species");
    }
    if (out.model.closure == Closure::Totals) {
        if (adjustableCount != 0 || !out.model.adjustedBasis.empty())
            throw std::invalid_argument("Totals closure does not allow adjustable primary components");
    } else {
        if (adjustableCount != 1)
            throw std::invalid_argument("Charge / fixed-pH closure requires exactly one adjustable component");
        if (out.model.closure == Closure::FixedPH) {
            if (out.adjustedIndex != out.hydrogenIndex)
                throw std::invalid_argument("Fixed-pH closure requires H+ as adjustable primary");
            if (!std::isfinite(out.model.fixedPH) || out.model.fixedPH < 0.0 ||
                out.model.fixedPH > 14.0)
                throw std::invalid_argument("Fixed pH must be between 0 and 14");
        }
    }
    for (const auto& f:out.model.formations) {
        if (f.name.empty() || !index.emplace(f.name,out.species.size()).second)
            throw std::invalid_argument("Duplicate or empty species name: " + f.name);
        if (!std::isfinite(f.logK25) || !std::isfinite(f.deltaHkJmol) || f.powers.empty())
            throw std::invalid_argument("Incomplete formation definition: " + f.name);
        Species s{f.name,f.charge,std::vector<int>(n,0),
                  temperatureLogK(f.logK25,f.deltaHkJmol,out.model.temperatureC)};
        long long charge=0;
        for (const auto& [name,coef]:f.powers) {
            const auto it=index.find(name);
            if (it==index.end() || it->second >= n || coef == 0)
                throw std::invalid_argument("Invalid formation basis or coefficient in " + f.name);
            s.powers[it->second]=coef;
            charge += static_cast<long long>(coef)*out.model.basis[it->second].charge;
        }
        if (charge != f.charge)
            throw std::invalid_argument("Formation electric charge not conserved: " + f.name);
        out.species.push_back(std::move(s));
    }
    return out;
}
Model readModel(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Cannot open model file: " + filename);
    Model model;
    bool haveT=false, haveActivity=false, haveClosure=false;
    std::string line;
    int row=0;
    while (std::getline(in,line)) {
        ++row;
        const auto first=line.find_first_not_of(" \t\r");
        if (first == std::string::npos || line[first]=='#') continue;
        std::istringstream is(line);
        std::string kind; is>>kind;
        if (kind=="temperature_C") {
            if (haveT || !(is>>model.temperatureC))
                throw std::runtime_error("Invalid temperature on line " + std::to_string(row));
            haveT=true;
        } else if (kind=="activity_model") {
            std::string value;
            if (haveActivity || !(is>>value))
                throw std::runtime_error("Invalid activity model on line " + std::to_string(row));
            if(value=="ideal") model.activity=ActivityModel::Ideal;
            else if(value=="davies") model.activity=ActivityModel::Davies;
            else throw std::runtime_error("Expected ideal or davies on line " + std::to_string(row));
            haveActivity=true;
        } else if (kind=="closure") {
            std::string value;
            if (haveClosure || !(is>>value))
                throw std::runtime_error("Invalid closure on line " + std::to_string(row));
            if(value=="totals") model.closure=Closure::Totals;
            else if(value=="charge") {
                model.closure=Closure::Charge;
                if(!(is>>model.adjustedBasis))
                    throw std::runtime_error("Charge closure needs adjustable species");
            } else if(value=="fixed_ph") {
                model.closure=Closure::FixedPH;
                if(!(is>>model.adjustedBasis>>model.fixedPH))
                    throw std::runtime_error("fixed_ph requires H+ and pH target");
            } else throw std::runtime_error("Unknown closure at line " + std::to_string(row));
            haveClosure=true;
        } else if (kind=="basis") {
            Basis b;
            std::string total;
            if(!(is>>b.name>>b.charge>>total))
                throw std::runtime_error("Invalid basis row " + std::to_string(row));
            if(total=="adjust") b.adjustable=true;
            else {
                std::size_t used=0;
                try {b.totalMolal=std::stod(total,&used);}
                catch(...) {throw std::runtime_error("Invalid molal total on line " + std::to_string(row));}
                if(used!=total.size())
                    throw std::runtime_error("Invalid molal total on line " + std::to_string(row));
            }
            model.basis.push_back(std::move(b));
        } else if (kind=="species") {
            Formation s;
            int count=0;
            if(!(is>>s.name>>s.charge>>s.logK25>>s.deltaHkJmol>>count) ||
                count<=0 || count>64)
                throw std::runtime_error("Invalid species row " + std::to_string(row));
            for(int i=0;i<count;++i) {
                std::string name; int coefficient;
                if(!(is>>name>>coefficient) || !s.powers.emplace(name,coefficient).second)
                    throw std::runtime_error("Invalid species term on line " + std::to_string(row));
            }
            model.formations.push_back(std::move(s));
        } else {
            throw std::runtime_error("Unknown keyword " + kind + " at line " + std::to_string(row));
        }
        endOfLine(is,row);
    }
    if(!haveT || !haveActivity || !haveClosure)
        throw std::runtime_error("Model must declare temperature_C, activity_model and closure");
    prepare(model);
    return model;
}
} // namespace thermo
