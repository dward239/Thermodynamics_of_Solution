#pragma once
#include "equilibrium/equilibrium.h"
#include <string>

namespace thermo {
void saveJson(const Result& result,const std::string& filename);
void saveCsv(const Result& result,const std::string& filename);
} // namespace thermo
