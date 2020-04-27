/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "RecurringProfileCalculator.h"
#include "ProfileEnergies.h"
#include <Utils/CalculatorBasics/StatesHandler.h>
#include <Utils/Geometry/AtomCollection.h>
#include <cassert>

namespace Scine {
namespace Readuct {

namespace ElementaryStepOptimization {

RecurringProfileCalculator::RecurringProfileCalculator(Core::Calculator& calculator, int numberEquidistantPoints)
  : calculator_(calculator), profileCalculator_(calculator), densities_(numberEquidistantPoints, nullptr) {
  valuesAlongSpline_.uValues = PointSequence{0.0, 1.0, numberEquidistantPoints};
  valuesAlongSpline_.energies.resize(numberEquidistantPoints);
  valuesAlongSpline_.gradients.resize(numberEquidistantPoints);
  deltaU_ = pointSequence().interval();
}

const double& RecurringProfileCalculator::deltaU() const {
  return deltaU_;
}

const std::vector<double>& RecurringProfileCalculator::getCoordinates() const {
  return pointSequence().underlyingArray();
}

const std::vector<double>& RecurringProfileCalculator::getEnergies() const {
  return valuesAlongSpline().energies;
}

void RecurringProfileCalculator::calculateEnergies(const Utils::BSplines::BSpline& spline) {
  calculateEnergies(spline, valuesAlongSpline_.energies);
}

void RecurringProfileCalculator::calculateEnergies(const Utils::BSplines::BSpline& spline, std::vector<double>& energies) {
  assert(static_cast<int>(energies.size()) == pointCount());

  for (int i = 0; i < pointCount(); ++i) {
    initializeDensity(i, spline);
    injectDensity(i);
    profileCalculator_.calculateEnergy(spline, pointSequence()[i], energies[i]);
    saveDensity(i);
  }
}

void RecurringProfileCalculator::calculateEnergiesAndGradients(const Utils::BSplines::BSpline& spline) {
  calculateEnergiesAndGradients(spline, valuesAlongSpline_.energies, valuesAlongSpline_.gradients);
}

void RecurringProfileCalculator::calculateEnergiesAndGradients(const Utils::BSplines::BSpline& spline,
                                                               std::vector<double>& energies,
                                                               std::vector<Utils::GradientCollection>& gradients) {
  assert(static_cast<int>(energies.size()) == pointCount());
  assert(static_cast<int>(gradients.size()) == pointCount());

  for (int i = 0; i < pointCount(); ++i) {
    initializeDensity(i, spline);
    injectDensity(i);
    profileCalculator_.calculateEnergyAndGradients(spline, pointSequence()[i], energies[i], gradients[i]);
    saveDensity(i);
  }
}

void RecurringProfileCalculator::calculateUpToSecondDerivative(const Utils::BSplines::BSpline& spline,
                                                               std::vector<double>& energies,
                                                               std::vector<Utils::GradientCollection>& gradients,
                                                               std::vector<Eigen::MatrixXd>& hessians) {
  assert(static_cast<int>(energies.size()) == pointCount());
  assert(static_cast<int>(gradients.size()) == pointCount());
  assert(static_cast<int>(hessians.size()) == pointCount());

  for (int i = 0; i < pointCount(); ++i) {
    initializeDensity(i, spline);
    injectDensity(i);
    profileCalculator_.calculateUpToSecondDerivatives(spline, pointSequence()[i], energies[i], gradients[i], hessians[i]);
    saveDensity(i);
  }
}

void RecurringProfileCalculator::calculateEnergyAndGradients(const Utils::BSplines::BSpline& spline, double u,
                                                             double& energy, Utils::GradientCollection& gradients) {
  profileCalculator_.calculateEnergyAndGradients(spline, u, energy, gradients);
}

void RecurringProfileCalculator::injectDensity(int index) {
  assert(0 <= index && index < pointCount());
  calculator_.loadState(densities_[index]);
}

void RecurringProfileCalculator::initializeDensity(int index, const Utils::BSplines::BSpline& spline) {
  assert(0 <= index && index < pointCount());
  if (densities_[index] == nullptr) {
    auto atoms = calculator_.getStructure();
    auto positions = profileCalculator_.positionsForUValue(spline, valuesAlongSpline_.uValues[index]);
    atoms->setPositions(positions);
    calculator_.setStructure(*atoms);
    densities_[index] = calculator_.getState();
  }
}

void RecurringProfileCalculator::saveDensity(int index) {
  assert(0 <= index && index < pointCount());
  densities_[index] = calculator_.getState();
}

ProfileEnergies RecurringProfileCalculator::getProfileEnergies() const {
  const auto& x = valuesAlongSpline().uValues.underlyingArray();
  const auto& y = valuesAlongSpline().energies;
  return {x, y};
}

} // namespace ElementaryStepOptimization
} // namespace Readuct
} // namespace Scine
