/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_REACTIONPATHCOSTCALCULATOR_H
#define READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_REACTIONPATHCOSTCALCULATOR_H

#include "../MSVCCompatibility.h"
#include <Eigen/Core>
#include <memory>

namespace Scine {
namespace Utils {
namespace BSplines {
class BSpline;
}
} // namespace Utils
namespace Readuct {
namespace ElementaryStepOptimization {
struct EnergiesAndGradientsAlongSpline;

namespace CostBasedOptimization {

/**
 * @brief Interface for the cost calculation of reaction paths.
 *
 * Every B-spline curve is associated with a cost that is then minimized. For this, the B-spline is
 * discretized into a number of points, and a local cost is associated with each point. The total
 * cost is then obtained from integrating these local costs along the spline. In our case, the cost
 * is made up of two terms. The first term is the electronic energy, integrated from the start of the
 * spline to its end, while the second term is a tension term which ensures that the points are not
 * distributed unevenly along the spline.
 *
 * For more information, see the Ph.D. thesis of Alain Vaucher, DOI 10.3929/ethz-b-000265855.
 */
class SCINE_DLLEXPORT ReactionPathCostCalculator {
 public:
  std::unique_ptr<ReactionPathCostCalculator> clone() const;

  /**
   * @brief Default destructor.
   */
  virtual ~ReactionPathCostCalculator() = default;

  /**
   * @brief Whether energies (and their gradients) are required in the cost calculation.
   *
   * @return Whether energies (and their gradients) are required in the cost calculation.
   */
  bool energiesRequired() const;

  /**
   * @brief Evaluate the cost associated with a given spline.
   *
   * @param spline The B-spline for which the costs should be evaluated.
   * @param energyValues The values of the electronic energy and its gradient along the B-spline.
   */
  void calculateCost(const Utils::BSplines::BSpline& spline, const EnergiesAndGradientsAlongSpline& energyValues);

  /**
   * @brief Get the cost associated with a given spline (after it has been evaluated).
   *
   * @return The cost associated with a given spline.
   */
  double getCost() const;

  /**
   * @brief Get the cost derivatives associated with a given spline.
   * @return The cost derivatives associated with a given spline.
   */
  Eigen::MatrixXd getCostDerivatives() const;

 private:
  virtual std::unique_ptr<ReactionPathCostCalculator> cloneImpl() const = 0;
  virtual bool energiesRequiredImpl() const = 0;
  virtual void calculateCostImpl(const Utils::BSplines::BSpline& spline,
                                 const EnergiesAndGradientsAlongSpline& energyValues) = 0;
  virtual double getCostImpl() const = 0;
  virtual Eigen::MatrixXd getCostDerivativesImpl() const = 0;
};

inline std::unique_ptr<ReactionPathCostCalculator> ReactionPathCostCalculator::clone() const {
  return cloneImpl();
}

inline bool ReactionPathCostCalculator::energiesRequired() const {
  return energiesRequiredImpl();
}

inline void ReactionPathCostCalculator::calculateCost(const Utils::BSplines::BSpline& spline,
                                                      const EnergiesAndGradientsAlongSpline& energyValues) {
  calculateCostImpl(spline, energyValues);
}

inline double ReactionPathCostCalculator::getCost() const {
  return getCostImpl();
}

inline Eigen::MatrixXd ReactionPathCostCalculator::getCostDerivatives() const {
  return getCostDerivativesImpl();
}

} // namespace CostBasedOptimization

} // namespace ElementaryStepOptimization

} // namespace Readuct
} // namespace Scine
#endif // READUCT_ELEMENTARYSTEPOPTIMIZATION_COSTCALCULATORS_REACTIONPATHCOSTCALCULATOR_H
