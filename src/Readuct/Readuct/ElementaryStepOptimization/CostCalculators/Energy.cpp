/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Energy.h"
#include "CommonTerms.h"
#include <Utils/Math/AutomaticDifferentiation/FirstND.h>
#include <Utils/Math/BSplines/BSpline.h>
#include <Utils/Typenames.h>

namespace Scine {
namespace Readuct {

namespace ElementaryStepOptimization {

namespace CostBasedOptimization {

Utils::AutomaticDifferentiation::FirstND Energy::calculateCostContribution(const Utils::BSplines::BSpline& spline,
                                                                           double u, double energy,
                                                                           const Utils::GradientCollection& gradients) const {
  Eigen::VectorXd energyGradientVector =
      Eigen::Map<const Eigen::VectorXd>(gradients.data(), gradients.cols() * gradients.rows());

  return CommonTerms::singlePointEnergy(spline, u, energy, energyGradientVector);
}

std::unique_ptr<ReactionPathCostCalculator> Energy::cloneImpl() const {
  return std::make_unique<Energy>(*this);
}

bool Energy::energiesRequiredImpl() const {
  return true;
}

} // namespace CostBasedOptimization

} // namespace ElementaryStepOptimization
} // namespace Readuct
} // namespace Scine
