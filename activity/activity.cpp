#include "activity/activity.h"
#include <cmath>
#include <stdexcept>

namespace thermo {
double waterDensityKgM3(double t) {
    if (!std::isfinite(t) || t<0.0 || t>80.0)
        throw std::invalid_argument("Water property temperature outside 0-80 C");
    // Pure liquid water near 1 atm; formula tabulated by ORNL/MEASUR.
    return 1000.0*(1.0 - ((t+288.9414)/(508929.2*(t+68.12963))) *
                                      (t-3.9863)*(t-3.9863));
}
double waterDielectric(double t) {
    if (!std::isfinite(t) || t<0.0 || t>80.0)
        throw std::invalid_argument("Water dielectric temperature outside 0-80 C");
    // Malmberg & Maryott (1956), NBS research paper 2641, 0-100 Celsius.
    return 87.740-0.40008*t+9.398e-4*t*t-1.410e-6*t*t*t;
}
double debyeHuckelA(double t) {
    const double densityGcm3=waterDensityKgM3(t)/1000.0;
    const double epsilon=waterDielectric(t);
    const double absoluteT=t+273.15;
    // Molality-based A; density in g/cm3 (Reaktoro DH model).
    return 1.824829238e6 * std::sqrt(densityGcm3) /
           std::pow(epsilon*absoluteT,1.5);
}
double daviesLogGamma(int charge, double ionicStrengthMolal, double t) {
    if (!std::isfinite(ionicStrengthMolal) || ionicStrengthMolal < 0.0 ||
        ionicStrengthMolal > 0.1)
        throw std::out_of_range("Davies ionic strength outside 0-0.1 mol/kg");
    if(charge==0) return 0.0; // Neutral species idealized, no Setschenow term.
    const double root=std::sqrt(ionicStrengthMolal);
    return -debyeHuckelA(t) * static_cast<double>(charge*charge) *
           (root/(1.0+root)-0.3*ionicStrengthMolal);
}
double activityLogGamma(ActivityModel model, int charge, double I, double t) {
    if(model==ActivityModel::Ideal) return 0.0;
    return daviesLogGamma(charge,I,t);
}
} // namespace thermo
