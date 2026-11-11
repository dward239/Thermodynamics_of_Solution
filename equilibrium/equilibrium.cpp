#include "equilibrium/equilibrium.h"
#include "activity/activity.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace thermo {
namespace {
constexpr double Ln10 = 2.302585092994045684;
using Vector = std::vector<double>;
using Matrix = std::vector<Vector>;
struct State {
    Vector concentrations;
    Vector gamma;
    Vector activities;
    Vector totals;
    Vector residuals;
    double I = 0.0;
    double calculatedI = 0.0;
    double charge = 0.0;
    double norm = std::numeric_limits<double>::infinity();
};
bool linearSolve(Matrix a, Vector b, Vector& dx) {
    const std::size_t n=b.size();
    for(std::size_t j=0;j<n;++j) {
        std::size_t pivot=j;
        for(std::size_t i=j+1;i<n;++i)
            if(std::abs(a[i][j])>std::abs(a[pivot][j])) pivot=i;
        if(!std::isfinite(a[pivot][j]) || std::abs(a[pivot][j])<1e-14)
            return false;
        std::swap(a[pivot],a[j]);
        std::swap(b[pivot],b[j]);
        for(std::size_t i=j+1;i<n;++i) {
            const double ratio=a[i][j]/a[j][j];
            for(std::size_t k=j+1;k<n;++k) a[i][k]-=ratio*a[j][k];
            b[i]-=ratio*b[j];
        }
    }
    dx.assign(n,0.0);
    for(std::size_t j=n;j-- >0;) {
        double value=b[j];
        for(std::size_t k=j+1;k<n;++k) value-=a[j][k]*dx[k];
        dx[j]=value/a[j][j];
        if(!std::isfinite(dx[j])) return false;
    }
    return true;
}
bool evaluate(const Prepared& prepared, const Vector& x, State& state) {
    const auto& model=prepared.model;
    const std::size_t n=model.basis.size();
    const bool nonideal=model.activity==ActivityModel::Davies;
    if(x.size()!=n+(nonideal?1u:0u)) return false;
    state=State{};
    state.I=nonideal?std::pow(10.0,x[n]):0.0;
    if(nonideal && (!std::isfinite(state.I) || state.I>0.1 || state.I<=0.0)) return false;
    const std::size_t ns=prepared.species.size();
    state.concentrations.assign(ns,0.0);
    state.gamma.assign(ns,1.0);
    state.activities.assign(ns,0.0);
    state.totals.assign(n,0.0);
    Vector logGamma(ns,0.0);
    try {
        for(std::size_t s=0;s<ns;++s) {
            logGamma[s]=activityLogGamma(model.activity,prepared.species[s].charge,
                                        state.I,model.temperatureC);
            state.gamma[s]=std::pow(10.0,logGamma[s]);
        }
    } catch (const std::exception&) {return false;}
    for(std::size_t s=0;s<ns;++s) {
        const auto& species=prepared.species[s];
        double logMolality=species.logK;
        if(s<n) logMolality=x[s];
        else {
            for(std::size_t j=0;j<n;++j)
                logMolality+=species.powers[j]*(x[j]+logGamma[j]);
            logMolality-=logGamma[s];
        }
        if(!std::isfinite(logMolality) || logMolality>70.0 || logMolality< -300.0) return false;
        const double m=std::pow(10.0,logMolality);
        if(!std::isfinite(m)) return false;
        state.concentrations[s]=m;
        state.activities[s]=m*state.gamma[s];
        state.charge += prepared.species[s].charge*m;
        state.calculatedI += 0.5*prepared.species[s].charge*
                              prepared.species[s].charge*m;
        for(std::size_t j=0;j<n;++j)
            state.totals[j]+=species.powers[j]*m;
    }
    // One closure equation per primary unknown: fixed totals except adjustable index.
    const double chargeScale=std::max(1e-8,[&](){
        double sum=0.0;
        for(std::size_t s=0;s<ns;++s)
            sum+=std::abs(prepared.species[s].charge)*state.concentrations[s];
        return sum;
    }());
    for(std::size_t j=0;j<n;++j) {
        if(static_cast<int>(j)==prepared.adjustedIndex) {
            if(model.closure==Closure::Charge)
                state.residuals.push_back(state.charge/chargeScale);
            else
                state.residuals.push_back(x[j]+logGamma[j]+model.fixedPH);
        } else {
            const double total=model.basis[j].totalMolal;
            state.residuals.push_back((state.totals[j]-total)/total);
        }
    }
    if(nonideal) {
        // Relative self-consistency of ionic strength.
        state.residuals.push_back((state.calculatedI-state.I)/
                                  std::max(state.calculatedI,state.I));
    }
    state.norm=0.0;
    for(double v:state.residuals) {
        if(!std::isfinite(v)) return false;
        state.norm=std::max(state.norm,std::abs(v));
    }
    return true;
}
}
Result solve(const Prepared& prepared, const SolverOptions& options) {
    // Repeat validation for callers who assembled a Prepared instance manually.
    const Prepared checked=prepare(prepared.model);
    const auto& model=checked.model;
    if(options.maxIterations<=0 || !std::isfinite(options.tolerance) ||
       options.tolerance<=0.0 || options.tolerance>=1.0)
        throw std::invalid_argument("Invalid solver iteration/tolerance settings");
    const std::size_t n=model.basis.size();
    const bool nonideal=model.activity==ActivityModel::Davies;
    Vector x(n+(nonideal?1u:0u));
    for(std::size_t j=0;j<n;++j) {
        double guess=model.basis[j].adjustable?1e-7:model.basis[j].totalMolal;
        x[j]=std::log10(guess);
    }
    if(checked.hydrogenIndex>=0) {
        if(model.closure==Closure::FixedPH)
            x[checked.hydrogenIndex]=-model.fixedPH;
        else if(model.closure==Closure::Charge)
            x[checked.hydrogenIndex]=-7.0;
    }
    if(nonideal) {
        double initialI=1e-6;
        for(std::size_t j=0;j<n;++j)
            initialI+=0.5*model.basis[j].charge*model.basis[j].charge*
                      (model.basis[j].adjustable?1e-7:model.basis[j].totalMolal);
        x[n]=std::log10(std::min(0.05,initialI));
    }
    State current;
    if(!evaluate(checked,x,current))
        throw std::runtime_error("Invalid initial state for equilibrium solve");
    int iterations=0;
    while(current.norm>options.tolerance && iterations<options.maxIterations) {
        const std::size_t dimensions=x.size();
        Matrix jac(dimensions,Vector(dimensions));
        for(std::size_t col=0;col<dimensions;++col) {
            // Centered differences in logarithmic coordinates, with a safeguarded fallback.
            constexpr double h=1e-5;
            Vector plus=x,minus=x;
            plus[col]+=h; minus[col]-=h;
            State p,m;
            const bool hasP=evaluate(checked,plus,p);
            const bool hasM=evaluate(checked,minus,m);
            if(!hasP && !hasM)
                throw std::runtime_error("Cannot construct finite-difference Jacobian");
            for(std::size_t row=0;row<dimensions;++row) {
                if(hasP && hasM) jac[row][col]=(p.residuals[row]-m.residuals[row])/(2*h);
                else if(hasP) jac[row][col]=(p.residuals[row]-current.residuals[row])/h;
                else jac[row][col]=(current.residuals[row]-m.residuals[row])/h;
            }
        }
        Vector rhs=current.residuals,dx;
        for(double& value:rhs) value=-value;
        if(!linearSolve(jac,rhs,dx))
            throw std::runtime_error("Equilibrium Jacobian singular or ill-conditioned");
        bool accepted=false;
        for(int attempt=0;attempt<50;++attempt) {
            const double alpha=std::ldexp(1.0,-attempt);
            Vector proposed=x;
            for(std::size_t j=0;j<dimensions;++j)
                proposed[j]+=alpha*dx[j];
            State trial;
            if(evaluate(checked,proposed,trial) &&
                trial.norm < current.norm*(1.0-1e-4*alpha)) {
                x=std::move(proposed);
                current=std::move(trial);
                accepted=true;
                break;
            }
        }
        if(!accepted)
            throw std::runtime_error("Equilibrium backtracking line search failed");
        ++iterations;
    }
    if(current.norm>options.tolerance)
        throw std::runtime_error("Equilibrium solver failed convergence criterion");
    Result result;
    result.iterations=iterations;
    result.temperatureC=model.temperatureC;
    result.ionicStrength=current.calculatedI;
    result.maxScaledResidual=current.norm;
    result.netCharge=current.charge;
    for(std::size_t s=0;s<checked.species.size();++s) {
        const auto& name=checked.species[s].name;
        result.molalities[name]=current.concentrations[s];
        result.activities[name]=current.activities[s];
        result.gammas[name]=current.gamma[s];
    }
    for(std::size_t j=0;j<n;++j) {
        const auto& b=model.basis[j];
        result.calculatedTotals[b.name]=current.totals[j];
        if(!b.adjustable)
            result.componentResiduals[b.name]=current.totals[j]-b.totalMolal;
    }
    if(checked.hydrogenIndex>=0)
        result.pH=-std::log10(current.activities[checked.hydrogenIndex]);
    else
        result.pH=std::numeric_limits<double>::quiet_NaN();
    return result;
}
} // namespace thermo
