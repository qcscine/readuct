/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_INTEGRATINGCOSTCALCULATOR_H
#define READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_INTEGRATINGCOSTCALCULATOR_H

#include "ReactionPathCostCalculator.h"
#include <Utils/Math/AutomaticDifferentiation/FirstND.h>
#include <Utils/Typenames.h>

// TODO: Improve documentation

namespace Scine {
namespace Readuct {

namespace ElementaryStepOptimization {

namespace CostBasedOptimization {

/*!
 * ReactionPathCostCalculator relying on an integration along the b-spline.
 */
class SCINE_DLLEXPORT IntegratingCostCalculator : public ReactionPathCostCalculator {
 private:
  void calculateCostImpl(const Utils::BSplines::BSpline& spline, const EnergiesAndGradientsAlongSpline& energyValues) override;
  double getCostImpl() const override;
  Eigen::MatrixXd getCostDerivativesImpl() const override;

  /*! To be overwritten by derived classes if necessary; is called at the beginning of the cost calculation of a path.
   */
  virtual void initializeCostCalculation(const Utils::BSplines::BSpline& /*spline*/) {
  }
  virtual Utils::AutomaticDifferentiation::FirstND
  calculateCostContribution(const Utils::BSplines::BSpline& spline, double u, double energy,
                            const Utils::GradientCollection& gradients) const = 0;

  Utils::AutomaticDifferentiation::FirstND totalCost_;
};

} // namespace CostBasedOptimization

} // namespace ElementaryStepOptimization

} // namespace Readuct
} // namespace Scine
#endif // READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_INTEGRATINGCOSTCALCULATOR_H
