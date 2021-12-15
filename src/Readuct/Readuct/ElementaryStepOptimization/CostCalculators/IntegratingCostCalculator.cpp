/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "IntegratingCostCalculator.h"
#include "Readuct/ElementaryStepOptimization/EnergiesAndGradientsAlongSpline.h"
#include <Utils/Math/BSplines/BSpline.h>

namespace Scine {
namespace Readuct {

using namespace Utils::AutomaticDifferentiation;

namespace ElementaryStepOptimization {

namespace CostBasedOptimization {

void IntegratingCostCalculator::calculateCostImpl(const Utils::BSplines::BSpline& spline,
                                                  const EnergiesAndGradientsAlongSpline& energyValues) {
  initializeCostCalculation(spline);

  const auto& points = energyValues.uValues;
  const auto& energies = energyValues.energies;
  const auto& gradients = energyValues.gradients;

  auto numberControlPoints = spline.controlPointCount();

  totalCost_ = FirstND{0, Eigen::MatrixXd::Zero(numberControlPoints, spline.getDim())};

  // Approximation: trapezoidal rule
  for (int i = 0; i < points.count(); ++i) {
    auto cost = calculateCostContribution(spline, points[i], energies[i], gradients[i]);

    if (i == 0 || i == points.count() - 1)
      cost /= 2;

    totalCost_ += cost;
  }

  totalCost_ *= points.interval();
}

double IntegratingCostCalculator::getCostImpl() const {
  return totalCost_.value();
}

Eigen::MatrixXd IntegratingCostCalculator::getCostDerivativesImpl() const {
  return totalCost_.derivatives();
}
} // namespace CostBasedOptimization

} // namespace ElementaryStepOptimization
} // namespace Readuct
} // namespace Scine
