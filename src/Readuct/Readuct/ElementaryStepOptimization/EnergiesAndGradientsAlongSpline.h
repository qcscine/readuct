/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef READUCT_ELEMENTARYSTEPOPTIMIZATION_ENERGIESANDGRADIENTSALONGSPLINE_H
#define READUCT_ELEMENTARYSTEPOPTIMIZATION_ENERGIESANDGRADIENTSALONGSPLINE_H

#include "PointSequence.h"
#include <Utils/Typenames.h>
#include <vector>

// TODO: Improve documentation

namespace Scine {
namespace Readuct {

namespace ElementaryStepOptimization {

/*!
 * Struct for values possibly required by reaction path cost calculators.
 */
struct EnergiesAndGradientsAlongSpline {
  PointSequence uValues{};
  std::vector<double> energies{};
  std::vector<Utils::GradientCollection> gradients{};
};

} // namespace ElementaryStepOptimization

} // namespace Readuct
} // namespace Scine
#endif // ELEMENTARYSTEPOPTIMIZATION_ENERGIESANDGRADIENTSALONGSPLINE_H