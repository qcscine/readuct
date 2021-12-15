/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_ENERGY_H
#define READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_ENERGY_H

#include "IntegratingCostCalculator.h"

// TODO: Improve documentation

namespace Scine {
namespace Readuct {

namespace ElementaryStepOptimization {

namespace CostBasedOptimization {

/*!
 * Cost calculator integrating the energy along the b-spline.
 */
class Energy : public IntegratingCostCalculator {
  Utils::AutomaticDifferentiation::FirstND calculateCostContribution(const Utils::BSplines::BSpline& spline, double u,
                                                                     double energy,
                                                                     const Utils::GradientCollection& gradients) const override;
  std::unique_ptr<ReactionPathCostCalculator> cloneImpl() const override;
  bool energiesRequiredImpl() const override;
};

} // namespace CostBasedOptimization

} // namespace ElementaryStepOptimization

} // namespace Readuct
} // namespace Scine
#endif // ELEMENTARYSTEPOPTIMIZATION_COSTBASEDOPTIMIZATION_ENERGY_H