/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "CommonTerms.h"
#include "Readuct/ElementaryStepOptimization/RecurringProfileCalculator.h"
#include <Utils/Math/AutomaticDifferentiation/FirstND.h>
#include <Utils/Math/BSplines/BSpline.h>
#include <Utils/Math/BSplines/ControlPointDerivatives.h>
#include <Utils/Typenames.h>
#include <Eigen/Dense>

namespace Scine {
namespace Readuct {

using namespace Utils::AutomaticDifferentiation;

namespace ElementaryStepOptimization {

namespace CostBasedOptimization {

FirstND CommonTerms::squaredNormOfSplineDerivative(const Utils::BSplines::BSpline& spline, double u) {
  Eigen::VectorXd splineDerivative = spline.evaluate(u, 1);
  double squaredNormOfDerivative = splineDerivative.squaredNorm();

  // d/dpAk (squaredNorm(C')) = 2 * C'k * dC'k/dpAk
  Eigen::MatrixXd cpTangentDerivatives = Utils::BSplines::ControlPointDerivatives::firstOrderCurveDerivatives(spline, u);
  Eigen::MatrixXd derivatives = 2 * cpTangentDerivatives * splineDerivative.asDiagonal();

  return FirstND(squaredNormOfDerivative, std::move(derivatives));
}

FirstND CommonTerms::singlePointEnergy(const Utils::BSplines::BSpline& spline, double u, double energy,
                                       const Eigen::VectorXd& gradients) {
  // dE/dpAk = dE/dCk * dCk/dpkA
  Eigen::MatrixXd cpDerivatives = Utils::BSplines::ControlPointDerivatives::curveDerivatives(spline, u);
  Eigen::MatrixXd derivatives = cpDerivatives * gradients.asDiagonal();

  return FirstND(energy, std::move(derivatives));
}

FirstND CommonTerms::squaredNormOfGradient(const Utils::BSplines::BSpline& spline, double u,
                                           const Eigen::VectorXd& gradients, const Eigen::MatrixXd& hessian) {
  double squaredNorm = gradients.squaredNorm();

  // d sqNorm(G) / dpAk = 2 * dCk/dpkA * G * H_k
  Eigen::VectorXd hgProduct = hessian * gradients;
  Eigen::MatrixXd cpDerivatives = Utils::BSplines::ControlPointDerivatives::curveDerivatives(spline, u);
  Eigen::MatrixXd derivatives = 2 * cpDerivatives * hgProduct.asDiagonal();

  return FirstND(squaredNorm, std::move(derivatives));
}

FirstND CommonTerms::dotProductOfGradientAndSplineDerivative(const Utils::BSplines::BSpline& spline, double u,
                                                             const Eigen::VectorXd& gradients, const Eigen::MatrixXd& hessian) {
  Eigen::VectorXd splineDerivative = spline.evaluate(u, 1);
  double dotProduct = gradients.dot(splineDerivative);

  // d (G.C) / dpAk = Gk * dC'k/dpAk + dCk/dpAk * C * H
  Eigen::MatrixXd hcProduct = hessian * splineDerivative;
  Eigen::MatrixXd cpTangentDerivatives = Utils::BSplines::ControlPointDerivatives::firstOrderCurveDerivatives(spline, u);
  Eigen::MatrixXd cpDerivatives = Utils::BSplines::ControlPointDerivatives::curveDerivatives(spline, u);
  Eigen::MatrixXd derivatives = cpTangentDerivatives * gradients.asDiagonal() + cpDerivatives * hcProduct.asDiagonal();

  return FirstND(dotProduct, std::move(derivatives));
}

FirstND CommonTerms::dotProductOfPositionAndSplineDerivative(const Utils::BSplines::BSpline& spline, double uPosition,
                                                             double uTangent) {
  Eigen::VectorXd position = spline.evaluate(uPosition);
  Eigen::VectorXd splineDerivative = spline.evaluate(uTangent, 1);
  double dotProduct = position.dot(splineDerivative);

  // d (C(uPos).C'(uTan)) / dpAk = dCk(uPos)/dPAk * C'k(uTan) + Ck(uPos) * dC'k(uTan)/dPAk
  Eigen::MatrixXd cpDerivatives = Utils::BSplines::ControlPointDerivatives::curveDerivatives(spline, uPosition);
  Eigen::MatrixXd cpTangentDerivatives =
      Utils::BSplines::ControlPointDerivatives::firstOrderCurveDerivatives(spline, uTangent);
  Eigen::MatrixXd derivatives = cpDerivatives * splineDerivative.asDiagonal() + cpTangentDerivatives * position.asDiagonal();

  return FirstND(dotProduct, std::move(derivatives));
}

FirstND CommonTerms::distanceSquared(const Utils::BSplines::BSpline& spline, double uA, double uB) {
  Eigen::VectorXd positionA = spline.evaluate(uA);
  Eigen::VectorXd positionB = spline.evaluate(uB);

  double distanceSquared = (positionA - positionB).squaredNorm();

  // d (Ck(uA) - Ck(uB))^2/dpk = 2 * (Ck(uA) - Ck(uB)) * (dCk(uA)/dpk - dCk(ub)/dpk)
  Eigen::MatrixXd cpDerivativesA = Utils::BSplines::ControlPointDerivatives::curveDerivatives(spline, uA);
  Eigen::MatrixXd cpDerivativesB = Utils::BSplines::ControlPointDerivatives::curveDerivatives(spline, uB);
  Eigen::MatrixXd derivatives = 2 * (cpDerivativesA - cpDerivativesB) * (positionA - positionB).asDiagonal();

  return FirstND(distanceSquared, std::move(derivatives));
}

FirstND CommonTerms::distanceSquaredToPosition(const Utils::BSplines::BSpline& spline, double u, const Eigen::VectorXd& pos) {
  Eigen::VectorXd splinePosition = spline.evaluate(u);

  double distanceSquared = (pos - splinePosition).squaredNorm();

  // d (Ck(u) - pos_k))^2/dpk = 2 * (Ck(uA) - pos_k) * dCk(uA)/dpk
  Eigen::MatrixXd cpDerivatives = Utils::BSplines::ControlPointDerivatives::curveDerivatives(spline, u);
  Eigen::MatrixXd derivatives = 2 * cpDerivatives * (splinePosition - pos).asDiagonal();

  return FirstND(distanceSquared, std::move(derivatives));
}

FirstND CommonTerms::dotProductOfPositions(const Utils::BSplines::BSpline& spline, double uA, double uB) {
  Eigen::VectorXd positionA = spline.evaluate(uA);
  Eigen::VectorXd positionB = spline.evaluate(uB);

  // d (C(uA).C(uB))/dpk = dCk(uA)/dpk*Ck(uB) + Ck(uA)*dCk(uB)/dpk
  Eigen::MatrixXd cpDerivativesA = Utils::BSplines::ControlPointDerivatives::curveDerivatives(spline, uA);
  Eigen::MatrixXd cpDerivativesB = Utils::BSplines::ControlPointDerivatives::curveDerivatives(spline, uB);

  Eigen::MatrixXd derivatives = cpDerivativesA * positionB.asDiagonal() + cpDerivativesB * positionA.asDiagonal();

  return FirstND(positionA.dot(positionB), std::move(derivatives));
}

FirstND CommonTerms::derivativeOfSquaredTangentNorm(const Utils::BSplines::BSpline& spline, double u) {
  Eigen::VectorXd firstDerivative = spline.evaluate(u, 1);
  Eigen::VectorXd secondDerivative = spline.evaluate(u, 2);

  /*! cost contribution:  */
  // d(||C'(u)||^2)/du) = 2 * dot(C',C'')
  // d/dpAk (d(||C'(u)||^2)/du)) = 2 * (dC'k/dpAk * C''k + C'k * dC''k/dpAk)
  double cost = 2 * firstDerivative.dot(secondDerivative);
  Eigen::MatrixXd cpFirstDerivatives = Utils::BSplines::ControlPointDerivatives::firstOrderCurveDerivatives(spline, u);
  Eigen::MatrixXd cpSecondDerivatives = Utils::BSplines::ControlPointDerivatives::secondOrderCurveDerivatives(spline, u);
  Eigen::MatrixXd derivatives =
      2.0 * (cpFirstDerivatives * secondDerivative.asDiagonal() + cpSecondDerivatives * firstDerivative.asDiagonal());

  return FirstND(cost, std::move(derivatives));
}

FirstND CommonTerms::energyOfHighestPointFromQuadraticInterpolation(RecurringProfileCalculator& profileCalculator,
                                                                    const std::vector<double>& energies,
                                                                    const Utils::BSplines::BSpline& spline) {
  const auto& coordinates = profileCalculator.getCoordinates();
  assert(coordinates.size() == energies.size());

  int maxIndex = getIndexForHighestEnergy(energies);
  auto numberPoints = static_cast<int>(energies.size());
  if (maxIndex == 0 || maxIndex == numberPoints - 1) {
    int numberControlPoints = spline.controlPointCount();
    auto cost = FirstND{energies[maxIndex], Eigen::MatrixXd::Zero(numberControlPoints, spline.getDim())};
    return cost; // I.e., all derivatives are zero
  }

  double uMax = getUValueWithMaxEnergy(coordinates, energies, maxIndex);

  return energyAlongSpline(profileCalculator, spline, uMax);
}

int CommonTerms::getIndexForHighestEnergy(const std::vector<double>& energies) {
  auto highestIndex = std::distance(energies.begin(), std::max_element(energies.begin(), energies.end()));
  return static_cast<int>(highestIndex);
}

FirstND CommonTerms::energyAlongSpline(RecurringProfileCalculator& profileCalculator,
                                       const Utils::BSplines::BSpline& spline, double u) {
  double energy;
  Utils::GradientCollection gc;
  profileCalculator.calculateEnergyAndGradients(spline, u, energy, gc);
  auto gradients = Eigen::Map<Eigen::VectorXd>(gc.data(), gc.cols() * gc.rows());

  // dE/ dpAk = dE/dCk * dCk/dpAk
  Eigen::MatrixXd cpDerivatives = Utils::BSplines::ControlPointDerivatives::curveDerivatives(spline, u);
  Eigen::MatrixXd derivatives = cpDerivatives * gradients.asDiagonal();

  return FirstND(energy, std::move(derivatives));
}

double CommonTerms::getUValueWithMaxEnergy(const std::vector<double>& coordinates, const std::vector<double>& energies,
                                           int maxIndex) {
  assert(coordinates.size() == energies.size());
  assert(maxIndex > 0 && maxIndex < static_cast<int>(energies.size()) - 1);

  // Solve a linear system of equations for the three points around the maximum
  // A x^2 + B x + C = y

  double u0 = coordinates[maxIndex - 1];
  double u1 = coordinates[maxIndex];
  double u2 = coordinates[maxIndex + 1];

  double e0 = energies[maxIndex - 1];
  double e1 = energies[maxIndex];
  double e2 = energies[maxIndex + 1];

  return CommonTerms::interpolateExtremumUValue(u0, u1, u2, e0, e1, e2);
}

double CommonTerms::interpolateExtremumUValue(double u0, double u1, double u2, double e0, double e1, double e2) {
  Eigen::Matrix3d A;
  A << u0 * u0, u0, 1, u1 * u1, u1, 1, u2 * u2, u2, 1;
  Eigen::Vector3d b;
  b << e0, e1, e2;

  // vector containing A, B, C
  Eigen::Vector3d coefficients = A.colPivHouseholderQr().solve(b);

  // vertex is (-B/2A, C - B^2/4A)
  double uMax = -coefficients(1) / (2 * coefficients(0));

  return uMax;
}

} // namespace CostBasedOptimization

} // namespace ElementaryStepOptimization
} // namespace Readuct
} // namespace Scine
