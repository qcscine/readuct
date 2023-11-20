/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef READUCT_ELEMENTARYSTEPOPTIMIZATION_RECURRINGPROFILECALCULATOR_H
#define READUCT_ELEMENTARYSTEPOPTIMIZATION_RECURRINGPROFILECALCULATOR_H

#include "BSplineProfileCalculator.h"
#include "EnergiesAndGradientsAlongSpline.h"
#include <Core/Interfaces/Calculator.h>
#include <Utils/Technical/UniqueIdentifier.h>
#include <vector>

// TODO: Improve documentation

namespace Scine {
namespace Utils {
class State;
} // namespace Utils
namespace Readuct {

namespace ElementaryStepOptimization {
class ProfileEnergies;

/*!
 * This class extends BSplineProfileCalculator for the case that the same single point calculations are
 * performed along an evolving BSpline, which permits optimization by saving the electronic densities along the path.
 * TODO: Allow for non-equidistant u values -> remove deltaU etc.
 */
class SCINE_DLLEXPORT RecurringProfileCalculator {
 public:
  explicit RecurringProfileCalculator(Core::Calculator& calculator, int numberEquidistantPoints);

  int pointCount() const;
  const double& deltaU() const;
  const std::vector<double>& getCoordinates() const;
  const std::vector<double>& getEnergies() const;
  const PointSequence& pointSequence() const;
  const EnergiesAndGradientsAlongSpline& valuesAlongSpline() const;
  ProfileEnergies getProfileEnergies() const;

  void calculateEnergies(const Utils::BSplines::BSpline& spline);
  void calculateEnergies(const Utils::BSplines::BSpline& spline, std::vector<double>& energies);

  void calculateEnergiesAndGradients(const Utils::BSplines::BSpline& spline);
  void calculateEnergiesAndGradients(const Utils::BSplines::BSpline& spline, std::vector<double>& energies,
                                     std::vector<Utils::GradientCollection>& gradients);

  void calculateUpToSecondDerivative(const Utils::BSplines::BSpline& spline, std::vector<double>& energies,
                                     std::vector<Utils::GradientCollection>& gradients,
                                     std::vector<Eigen::MatrixXd>& hessians);

  void calculateEnergyAndGradients(const Utils::BSplines::BSpline& spline, double u, double& energy,
                                   Utils::GradientCollection& gradients);

 private:
  void injectDensity(int index);
  void saveDensity(int index);

  Core::Calculator& calculator_;
  BSplineProfileCalculator profileCalculator_;

  double deltaU_;
  EnergiesAndGradientsAlongSpline valuesAlongSpline_;
  std::vector<std::shared_ptr<Core::State>> densities_;
};

inline int RecurringProfileCalculator::pointCount() const {
  return pointSequence().count();
}

inline const PointSequence& RecurringProfileCalculator::pointSequence() const {
  return valuesAlongSpline_.uValues;
}

inline const EnergiesAndGradientsAlongSpline& RecurringProfileCalculator::valuesAlongSpline() const {
  return valuesAlongSpline_;
}

} // namespace ElementaryStepOptimization

} // namespace Readuct
} // namespace Scine
#endif // ELEMENTARYSTEPOPTIMIZATION_RECURRINGPROFILECALCULATOR_H
