#include "aqueous_engine.h"
#include "aqueous_solid.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace aqueous {
namespace {
constexpr double kLn10 = 2.302585092994046;
constexpr double kGasConstant = 8.31446261815324; // J mol^-1 K^-1
struct Species {
    std::string name;
    int charge = 0;
    double logK = 0;
    std::vector<int> nu;
};
struct Context {
    std::vector<Species> species;
    std::vector<double> totals;
    std::vector<bool> fixed;
    double chargeScale = 1e-6;
    int freeIndex = -1;
};
struct Calculation {
    std::vector<double> concentration;
    std::vector<double> residual;
    double ionicStrength = 0.0;
    double charge = 0.0;
    double norm = std::numeric_limits<double>::infinity();
};
void checkTemperature(double T) {
    // Activity coefficients use a pure-water correlation within this range.
    if (!std::isfinite(T) || T < 273.15 || T > 323.15)
        throw std::invalid_argument("Temperature must be 273.15 to 323.15 K (0-50 C)");
}
double dielectric(double T) {
    const double c = T - 273.15;
    return 87.740 - 0.40008*c + 9.398e-4*c*c - 1.410e-6*c*c*c;
}
double density(double T) {
    // Pure-water correlation, approximate g/cm3 near atmospheric pressure.
    const double c = T - 273.15;
    return 1.0 - ((c + 288.9414)/(508929.2*(c+68.12963)))*(c-3.9863)*(c-3.9863);
}
double daviesA(double T) {
    const double epsilon = dielectric(T);
    return 1.82483e6 * std::sqrt(density(T)) / std::pow(epsilon*T, 1.5);
}
Context prepare(const Model& m) {
    checkTemperature(m.temperatureK);
    checkTemperature(m.referenceTemperatureK);
    if (m.components.empty()) throw std::invalid_argument("At least one component is required");
    Context ctx;
    const std::size_t n = m.components.size();
    std::map<std::string, std::size_t> indices;
    ctx.totals.resize(n);
    ctx.fixed.resize(n);
    for (std::size_t i=0; i<n; ++i) {
        const auto& c = m.components[i];
        if (c.name.empty() || !indices.emplace(c.name,i).second)
            throw std::invalid_argument("Duplicate or empty component name");
        if (c.fixedTotal) {
            if (!std::isfinite(c.totalMolL) || c.totalMolL<=0)
                throw std::invalid_argument("Fixed total must be positive and finite: "+c.name);
            ctx.chargeScale += std::abs(c.charge)*c.totalMolL;
        } else {
            if (ctx.freeIndex != -1) throw std::invalid_argument("Exactly zero or one free component is supported");
            ctx.freeIndex = static_cast<int>(i);
        }
        ctx.fixed[i]=c.fixedTotal;
        ctx.totals[i]=c.totalMolL;
        Species p;
        p.name=c.name;
        p.charge=c.charge;
        p.nu.assign(n,0);
        p.nu[i]=1;
        ctx.species.push_back(std::move(p));
    }
    for(const auto& c:m.complexes) {
        if(c.name.empty() || !indices.emplace(c.name,ctx.species.size()).second ||
            !std::isfinite(c.logKReference) || !std::isfinite(c.deltaHkJMol) ||
            c.basisStoichiometry.empty())
            throw std::invalid_argument("Invalid or duplicate complex: "+c.name);
        Species s;
        s.name=c.name;
        s.charge=c.charge;
        s.logK=logKAtTemperature(c,m.temperatureK,m.referenceTemperatureK);
        s.nu.assign(n,0);
        long long charge=0;
        for(const auto& [name,stoich]:c.basisStoichiometry) {
            auto it=indices.find(name);
            if(it==indices.end() || it->second>=n || stoich==0)
                throw std::invalid_argument("Invalid basis component in "+c.name+": "+name);
            s.nu[it->second]=stoich;
            charge += static_cast<long long>(stoich)*m.components[it->second].charge;
        }
        if(charge!=c.charge)
            throw std::invalid_argument("Charge is not conserved in formation of "+c.name);
        ctx.species.push_back(std::move(s));
    }
    return ctx;
}
bool calculate(const Model& m,const Context& ctx,const std::vector<double>& x,
               Calculation& out) {
    const std::size_t n=m.components.size();
    const bool isDavies=m.activity==ActivityModel::Davies;
    const double I=isDavies?std::pow(10.0,x[n]):0.0;
    if(isDavies && (!std::isfinite(I) || I<=0 || I>0.100000001)) return false;
    std::vector<double> lg(ctx.species.size(),0.0);
    if(isDavies) {
        for(std::size_t s=0;s<ctx.species.size();++s)
            lg[s]=daviesLogGamma(ctx.species[s].charge,I,m.temperatureK);
    }
    out.concentration.assign(ctx.species.size(),0.0);
    out.charge=0;
    out.ionicStrength=0;
    for(std::size_t s=0;s<ctx.species.size();++s) {
        const auto& sp=ctx.species[s];
        double logC=sp.logK-lg[s];
        for(std::size_t j=0;j<n;++j)
            logC += sp.nu[j]*(x[j]+lg[j]);
        // Prevent overflow and loss of meaningful mass during iteration.
        if(!std::isfinite(logC) || logC>200.0 || logC< -290.0) return false;
        const double c=std::pow(10.0,logC);
        if(!std::isfinite(c)) return false;
        out.concentration[s]=c;
        out.charge+=sp.charge*c;
        out.ionicStrength+=0.5*c*sp.charge*sp.charge;
    }
    out.residual.clear();
    for(std::size_t i=0;i<n;++i) {
        if(!ctx.fixed[i]) continue;
        double calculated=0;
        for(std::size_t s=0;s<ctx.species.size();++s)
            calculated += ctx.species[s].nu[i]*out.concentration[s];
        out.residual.push_back((calculated-ctx.totals[i])/ctx.totals[i]);
    }
    if(ctx.freeIndex!=-1) out.residual.push_back(out.charge/ctx.chargeScale);
    if(isDavies) out.residual.push_back(out.ionicStrength/I - 1.0);
    out.norm=0;
    for(double v:out.residual) {
        if(!std::isfinite(v)) return false;
        out.norm=std::max(out.norm,std::abs(v));
    }
    return true;
}
bool linearSolve(std::vector<std::vector<double>> a,std::vector<double> b,
                 std::vector<double>& answer) {
    const std::size_t n=b.size();
    for(std::size_t k=0;k<n;++k) {
        std::size_t p=k;
        for(std::size_t i=k+1;i<n;++i)
            if(std::abs(a[i][k])>std::abs(a[p][k])) p=i;
        if(!std::isfinite(a[p][k]) || std::abs(a[p][k])<1e-14) return false;
        std::swap(a[p],a[k]);std::swap(b[p],b[k]);
        for(std::size_t i=k+1;i<n;++i) {
            double f=a[i][k]/a[k][k];
            for(std::size_t j=k+1;j<n;++j) a[i][j]-=f*a[k][j];
            b[i]-=f*b[k];
        }
    }
    answer.assign(n,0);
    for(std::size_t i=n;i-- >0;) {
        double v=b[i];
        for(std::size_t j=i+1;j<n;++j) v-=a[i][j]*answer[j];
        answer[i]=v/a[i][i];
        if(!std::isfinite(answer[i])) return false;
    }
    return true;
}
std::string escapeJson(const std::string& s) {
    std::string out;
    for(char c:s) {
        if(c=='"' || c=='\\') out+='\\';
        if(c=='\n') out+="\\n";
        else out+=c;
    }
    return out;
}
}
double logKAtTemperature(const Complex& r,double T,double Tref) {
    checkTemperature(T);checkTemperature(Tref);
    return r.logKReference - (r.deltaHkJMol*1000.0/(kGasConstant*kLn10))*(1.0/T-1.0/Tref);
}
double daviesLogGamma(int charge,double I,double T) {
    checkTemperature(T);
    if(!std::isfinite(I) || I<0 || I>0.100000001)
        throw std::invalid_argument("Davies approximation restricted to ionic strength <= 0.1 mol/L");
    if(charge==0) return 0.0;
    const double root=std::sqrt(I);
    return -daviesA(T)*charge*charge*(root/(1.0+root)-0.3*I);
}
Model readModel(const std::string& filename) {
    std::ifstream file(filename);
    if(!file) throw std::runtime_error("Cannot open input: "+filename);
    Model m;
    std::string line;
    int number=0;
    while(std::getline(file,line)) {
        ++number;
        auto pos=line.find_first_not_of(" \t\r");
        if(pos==std::string::npos || line[pos]=='#') continue;
        std::istringstream in(line);
        std::string key;in>>key;
        if(key=="temperature_k") {
            if(!(in>>m.temperatureK)) throw std::runtime_error("Bad temperature");
        } else if(key=="reference_temperature_k") {
            if(!(in>>m.referenceTemperatureK)) throw std::runtime_error("Bad reference temperature");
        } else if(key=="activity_model") {
            std::string model;
            if(!(in>>model)) throw std::runtime_error("Bad activity model");
            if(model=="ideal") m.activity=ActivityModel::Ideal;
            else if(model=="davies") m.activity=ActivityModel::Davies;
            else throw std::runtime_error("Unknown activity model: "+model);
        } else if(key=="component") {
            Component c;
            std::string total;
            if(!(in>>c.name>>c.charge>>total)) throw std::runtime_error("Bad component row");
            if(total=="free") c.fixedTotal=false;
            else {
                try {
                    std::size_t length=0;
                    c.totalMolL=std::stod(total,&length);
                    if(length!=total.size()) throw std::invalid_argument("trailing");
                } catch(...) {throw std::runtime_error("Invalid total at line "+std::to_string(number));}
            }
            m.components.push_back(std::move(c));
        } else if(key=="solid") {
            Solid s; int count=0;
            if(!(in>>s.name>>s.logKspReference>>s.deltaHkJMol>>count) || count<1)
                throw std::runtime_error("Bad solid at line "+std::to_string(number));
            for(int i=0;i<count;++i) {
                std::string name; int nu=0;
                if(!(in>>name>>nu) || nu<=0 || !s.dissolutionStoichiometry.emplace(name,nu).second)
                    throw std::runtime_error("Bad dissolution stoichiometry at line "+std::to_string(number));
            }
            m.solids.push_back(std::move(s));
        } else if(key=="complex") {
            Complex c;int count=0;
            if(!(in>>c.name>>c.charge>>c.logKReference>>c.deltaHkJMol>>count) || count<1)
                throw std::runtime_error("Bad complex at line "+std::to_string(number));
            for(int i=0;i<count;++i) {
                std::string name;int v=0;
                if(!(in>>name>>v) || !c.basisStoichiometry.emplace(name,v).second)
                    throw std::runtime_error("Bad stoichiometry at line "+std::to_string(number));
            }
            m.complexes.push_back(std::move(c));
        } else throw std::runtime_error("Unknown keyword at line "+std::to_string(number)+": "+key);
        std::string extra;
        if(in>>extra && extra[0]!='#')
            throw std::runtime_error("Extra data at line "+std::to_string(number));
    }
    prepare(m);
    return m;
}
Result solve(const Model& model,const Options& options) {
    if(options.maxIterations<1 || !std::isfinite(options.relativeTolerance) ||
       options.relativeTolerance<=0 || options.relativeTolerance>=0.1)
        throw std::invalid_argument("Invalid solver options");
    if(!model.solids.empty()) return solveSingleIdealSolid(model, options);
    const Context ctx=prepare(model);
    const std::size_t n=model.components.size();
    const bool davies=model.activity==ActivityModel::Davies;
    std::vector<double> x(n+(davies?1:0),0.0);
    for(std::size_t i=0;i<n;++i) {
        const auto& c=model.components[i];
        x[i]=c.fixedTotal?std::log10(c.totalMolL)-0.2:
              (c.name=="H+"?-6.0:-6.0);
    }
    if(davies) {
        double initialI=1e-8;
        for(const auto& c:model.components)
            if(c.fixedTotal) initialI+=0.5*c.charge*c.charge*c.totalMolL;
        if(initialI>0.1) throw std::invalid_argument("Input likely exceeds Davies validity range");
        x[n]=std::log10(initialI);
    }
    Calculation current;
    if(!calculate(model,ctx,x,current))
        throw std::runtime_error("Initial equilibrium state invalid");
    int iterations=0;
    while(current.norm>options.relativeTolerance && iterations<options.maxIterations) {
        const std::size_t width=x.size();
        std::vector<std::vector<double>> J(width,std::vector<double>(width,0.0));
        for(std::size_t j=0;j<width;++j) {
            const double h=1e-5;
            auto test=x;test[j]+=h;
            Calculation next;
            if(!calculate(model,ctx,test,next)) {
                test[j]=x[j]-h;
                if(!calculate(model,ctx,test,next))
                    throw std::runtime_error("Jacobian evaluation failed");
                for(std::size_t i=0;i<width;++i) J[i][j]=(current.residual[i]-next.residual[i])/h;
            } else for(std::size_t i=0;i<width;++i)
                J[i][j]=(next.residual[i]-current.residual[i])/h;
        }
        auto rhs=current.residual;
        for(double& v:rhs) v=-v;
        std::vector<double> step;
        if(!linearSolve(J,rhs,step)) throw std::runtime_error("Singular equilibrium Jacobian");
        bool accepted=false;
        for(int attempt=0;attempt<45;++attempt) {
            const double alpha=std::ldexp(1.0,-attempt);
            auto trial=x;
            for(std::size_t j=0;j<width;++j) trial[j]+=alpha*step[j];
            Calculation candidate;
            if(calculate(model,ctx,trial,candidate) && candidate.norm<current.norm) {
                current=std::move(candidate);x=std::move(trial);accepted=true;break;
            }
        }
        if(!accepted) throw std::runtime_error("Equilibrium line search failed");
        ++iterations;
    }
    if(current.norm>options.relativeTolerance)
        throw std::runtime_error("Equilibrium failed to converge in specified iterations");
    Result result;
    result.temperatureK=model.temperatureK;result.model=model.activity;
    result.iterations=iterations;result.normalizedResidual=current.norm;
    result.ionicStrengthMolL=current.ionicStrength;result.chargeResidualMolL=current.charge;
    for(std::size_t s=0;s<ctx.species.size();++s) {
        const auto& species=ctx.species[s];
        const double c=current.concentration[s];
        result.concentrationMolL[species.name]=c;
        const double lg=davies?daviesLogGamma(species.charge,current.ionicStrength,model.temperatureK):0.0;
        result.activity[species.name]=c*std::pow(10.0,lg); // c/c°; c° = 1 mol/L
    }
    for(std::size_t i=0;i<n;++i) {
        if(!ctx.fixed[i]) continue;
        double actual=0;
        for(std::size_t s=0;s<ctx.species.size();++s)
            actual+=ctx.species[s].nu[i]*current.concentration[s];
        result.componentResidualMolL[model.components[i].name]=actual-ctx.totals[i];
    }
    auto h=result.activity.find("H+");
    if(h!=result.activity.end()) {
        result.hasPH=true;
        result.pH=-std::log10(h->second);
    }
    return result;
}
std::string toJson(const Result& r) {
    std::ostringstream out;out<<std::setprecision(15);
    out<<"{\n  \"temperature_K\": "<<r.temperatureK<<",\n"
       <<"  \"activity_model\": \""<<(r.model==ActivityModel::Ideal?"ideal":"davies")<<"\",\n"
       <<"  \"iterations\": "<<r.iterations<<",\n"
       <<"  \"normalized_residual\": "<<r.normalizedResidual<<",\n"
       <<"  \"ionic_strength_mol_L\": "<<r.ionicStrengthMolL<<",\n"
       <<"  \"charge_residual_mol_L\": "<<r.chargeResidualMolL<<",\n";
    if(r.hasPH) out<<"  \"pH\": "<<r.pH<<",\n";
    for(const auto* field : {"concentrations_mol_L","activities","component_residuals_mol_L",
                             "solids_precipitated_mol_L","saturation_indices"}) {
        const std::string key(field);
        const auto& data=key=="concentrations_mol_L"?r.concentrationMolL:
                         key=="activities"?r.activity:
                         key=="component_residuals_mol_L"?r.componentResidualMolL:
                         key=="solids_precipitated_mol_L"?r.solidPrecipitatedMolL:r.saturationIndex;
        out<<"  \""<<field<<"\": {";
        bool first=true;
        for(const auto& [name,value]:data) {
            if(!first) out<<",";
            first=false;
            out<<"\n    \""<<escapeJson(name)<<"\": "<<value;
        }
        if(!data.empty()) out<<"\n  ";
        out<<"}"<<(key=="saturation_indices"?"\n":",\n");
    }
    out<<"}\n";
    return out.str();
}
void writeCsv(const Result& r,const std::string& filename) {
    std::ofstream out(filename);
    if(!out) throw std::runtime_error("Cannot write CSV: "+filename);
    out<<"species,concentration_mol_L,dimensionless_activity\n";
    out<<std::setprecision(15);
    for(const auto& [name,c]:r.concentrationMolL)
        out<<name<<","<<c<<","<<r.activity.at(name)<<"\n";
}
}
