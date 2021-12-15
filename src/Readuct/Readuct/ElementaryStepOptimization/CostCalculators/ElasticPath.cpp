/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ElasticPath.h"
#include "CommonTerms.h"
#include <Utils/Math/BSplines/BSpline.h>
#include <Utils/Math/BSplines/ControlPointDerivatives.h>
#include <Utils/Typenames.h>

namespace Scine {
namespace Readuct {

using namespace Utils::AutomaticDifferentiation;

namespace ElementaryStepOptimization {

namespace CostBasedOptimization {

FirstND ElasticPath::calculateCostContribution(const Utils::BSplines::BSpline& spline, double u, double /*energy*/,
                                               const Utils::GradientCollection& /*gradients*/) const {
  auto cost = elasticContribution(spline, u);
  return cost;
}

Utils::AutomaticDifferentiation::FirstND ElasticPath::elasticContribution(const Utils::BSplines::BSpline& spline, double u) const {
  return square(CommonTerms::derivativeOfSquaredTangentNorm(spline, u));
}

bool ElasticPath::energiesRequiredImpl() const {
  return false;
}

std::unique_ptr<ReactionPathCostCalculator> ElasticPath::cloneImpl() const {
  return std::make_unique<ElasticPath>(*this);
}

} // namespace CostBasedOptimization

} // namespace ElementaryStepOptimization
} // namespace Readuct
} // namespace Scine
