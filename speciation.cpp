#include "speciation.h"
#include <stdexcept>

void solveSpeciation(std::vector<AqueousSpecies>&,
    Solution&, const TotalInput&) {
    throw std::logic_error(
        "Speciation is not implemented: requires explicit reactions, "
        "mass balances, activity iteration and numerical validation.");
}
