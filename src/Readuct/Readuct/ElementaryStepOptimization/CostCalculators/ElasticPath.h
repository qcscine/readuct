/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_ELASTICPATH_H
#define READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_ELASTICPATH_H

#include "IntegratingCostCalculator.h"

// TODO: Improve documentation

namespace Scine {
namespace Readuct {

namespace ElementaryStepOptimization {

namespace CostBasedOptimization {

/*!
 * Cost calculator for elasticity / tension.
 * Equivalent to the squared norm of the spline derivative at a given point.
 */
class SCINE_DLLEXPORT ElasticPath : public IntegratingCostCalculator {
  std::unique_ptr<ReactionPathCostCalculator> cloneImpl() const override;
  bool energiesRequiredImpl() const override;
  Utils::AutomaticDifferentiation::FirstND calculateCostContribution(const Utils::BSplines::BSpline& spline, double u,
                                                                     double energy,
                                                                     const Utils::GradientCollection& gradients) const override;
  Utils::AutomaticDifferentiation::FirstND elasticContribution(const Utils::BSplines::BSpline& spline, double u) const;
};

} // namespace CostBasedOptimization

} // namespace ElementaryStepOptimization

} // namespace Readuct
} // namespace Scine
#endif // ELEMENTARYSTEPOPTIMIZATION_COSTBASEDOPTIMIZATION_ELASTICPATH_H
