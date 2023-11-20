/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "BSplineProfileCalculator.h"
#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/Math/BSplines/BSpline.h>

namespace Scine {
namespace Readuct {

namespace ElementaryStepOptimization {

BSplineProfileCalculator::BSplineProfileCalculator(Core::Calculator& calculator) : calculator_(calculator) {
}

void BSplineProfileCalculator::calculateAllEnergies(const Utils::BSplines::BSpline& spline,
                                                    const std::vector<double>& uValues, std::vector<double>& energies) {
  assert(energies.size() == uValues.size());

  for (int i = 0; i < static_cast<int>(uValues.size()); ++i) {
    calculateEnergy(spline, uValues[i], energies[i]);
  }
}

void BSplineProfileCalculator::calculateEnergy(const Utils::BSplines::BSpline& spline, double u, double& energy) {
  calculator_.modifyPositions(positionsForUValue(spline, u));
  calculator_.setRequiredProperties(Utils::Property::Energy);
  auto cr = calculator_.calculate("");
  energy = cr.get<Utils::Property::Energy>();
}

void BSplineProfileCalculator::calculateAllEnergiesAndGradients(const Utils::BSplines::BSpline& spline,
                                                                const std::vector<double>& uValues, std::vector<double>& energies,
                                                                std::vector<Utils::GradientCollection>& gradients) {
  assert(energies.size() == uValues.size());
  assert(gradients.size() == uValues.size());

  for (int i = 0; i < static_cast<int>(uValues.size()); ++i) {
    calculateEnergyAndGradients(spline, uValues[i], energies[i], gradients[i]);
  }
}

void BSplineProfileCalculator::calculateEnergyAndGradients(const Utils::BSplines::BSpline& spline, double u,
                                                           double& energy, Utils::GradientCollection& gradients) {
  calculator_.modifyPositions(positionsForUValue(spline, u));
  calculator_.setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
  auto cr = calculator_.calculate("");
  energy = cr.get<Utils::Property::Energy>();
  gradients = cr.get<Utils::Property::Gradients>();
}

void BSplineProfileCalculator::calculateAllUpToSecondDerivative(const Utils::BSplines::BSpline& spline,
                                                                const std::vector<double>& uValues, std::vector<double>& energies,
                                                                std::vector<Utils::GradientCollection>& gradients,
                                                                std::vector<Eigen::MatrixXd>& hessians) {
  assert(energies.size() == uValues.size());
  assert(gradients.size() == uValues.size());
  assert(hessians.size() == uValues.size());

  for (int i = 0; i < static_cast<int>(uValues.size()); ++i) {
    calculateUpToSecondDerivatives(spline, uValues[i], energies[i], gradients[i], hessians[i]);
  }
}

void BSplineProfileCalculator::calculateUpToSecondDerivatives(const Utils::BSplines::BSpline& /*spline*/, double /*u*/,
                                                              double& /*energy*/, Utils::GradientCollection& /*gradients*/,
                                                              Eigen::MatrixXd& /*hessians*/) {
  throw std::runtime_error("Hessian calculation currently not available");
}

Utils::PositionCollection BSplineProfileCalculator::positionsForUValue(const Utils::BSplines::BSpline& spline, double u) const {
  Eigen::VectorXd positions = spline.evaluate(u);
  return Eigen::Map<const Utils::PositionCollection>(positions.data(), positions.size() / 3, 3);
}

} // namespace ElementaryStepOptimization
} // namespace Readuct
} // namespace Scine
