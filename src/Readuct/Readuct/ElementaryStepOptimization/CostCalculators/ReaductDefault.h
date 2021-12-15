/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_READUCTDEFAULT_H
#define READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_READUCTDEFAULT_H

#include "CostCombiner.h"
#include "ElasticPath.h"
#include "Energy.h"
#include "ReactionPathCostCalculator.h"

// TODO: Improve documentation

namespace Scine {
namespace Readuct {

namespace ElementaryStepOptimization {

namespace CostBasedOptimization {

/*!
 * Default ReaDuct cost calculator.
 */
class ReaductDefault : public ReactionPathCostCalculator {
 public:
  ReaductDefault();

  void setTensionFactor(double factor);
  double getTensionFactor() const;

  static constexpr double defaultTensionFactor = 1e-5;

 private:
  std::unique_ptr<ReactionPathCostCalculator> cloneImpl() const override;
  bool energiesRequiredImpl() const override;
  void calculateCostImpl(const Utils::BSplines::BSpline& spline, const EnergiesAndGradientsAlongSpline& energyValues) override;
  double getCostImpl() const override;
  Eigen::MatrixXd getCostDerivativesImpl() const override;

  CostCombiner<Energy, ElasticPath> combinedCostCalculator_;
};

} // namespace CostBasedOptimization

} // namespace ElementaryStepOptimization

} // namespace Readuct
} // namespace Scine
#endif // ELEMENTARYSTEPOPTIMIZATION_COSTBASEDOPTIMIZATION_READUCTDEFAULT_H
