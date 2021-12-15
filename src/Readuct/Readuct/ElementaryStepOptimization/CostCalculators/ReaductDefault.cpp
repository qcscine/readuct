/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ReaductDefault.h"

namespace Scine {
namespace Readuct {

namespace ElementaryStepOptimization {

namespace CostBasedOptimization {

ReaductDefault::ReaductDefault() {
  setTensionFactor(defaultTensionFactor);
}

void ReaductDefault::setTensionFactor(double factor) {
  double firstFactor = 1.0 - factor;
  combinedCostCalculator_.setFirstCalculatorContribution(firstFactor);
}

double ReaductDefault::getTensionFactor() const {
  double firstFactor = combinedCostCalculator_.getFirstCalculatorContribution();
  return 1.0 - firstFactor;
}

std::unique_ptr<ReactionPathCostCalculator> ReaductDefault::cloneImpl() const {
  return std::make_unique<ReaductDefault>(*this);
}

bool ReaductDefault::energiesRequiredImpl() const {
  return combinedCostCalculator_.energiesRequired();
}

void ReaductDefault::calculateCostImpl(const Utils::BSplines::BSpline& spline,
                                       const EnergiesAndGradientsAlongSpline& energyValues) {
  return combinedCostCalculator_.calculateCost(spline, energyValues);
}

double ReaductDefault::getCostImpl() const {
  return combinedCostCalculator_.getCost();
}

Eigen::MatrixXd ReaductDefault::getCostDerivativesImpl() const {
  return combinedCostCalculator_.getCostDerivatives();
}

} // namespace CostBasedOptimization

} // namespace ElementaryStepOptimization
} // namespace Readuct
} // namespace Scine
