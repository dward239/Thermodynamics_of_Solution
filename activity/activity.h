#pragma once
#include "chemistry/model.h"

namespace thermo {
// Pure-water properties near atmospheric pressure, used solely in Davies A(T).
double waterDensityKgM3(double temperatureC);
double waterDielectric(double temperatureC);
double debyeHuckelA(double temperatureC);
double daviesLogGamma(int charge, double ionicStrengthMolal, double temperatureC);
double activityLogGamma(ActivityModel model, int charge, double ionicStrengthMolal,
                        double temperatureC);
} // namespace thermo
