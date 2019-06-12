/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Readuct/BSplines/BSpline.h>
#include <Readuct/BSplines/InterpolationGenerator.h>
#include <Readuct/TransitionStateSearch/CostBasedOptimization/CostCalculators/CommonTerms.h>
#include <Utils/math/AutomaticDifferentiation/FirstND.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Readuct {
using namespace BSplines;
using namespace TransitionStateSearch::CostBasedOptimization;
namespace Tests {

class CommonCostTerms : public Test {
 public:
  BSpline randomBSpline;
  const int numberDimensions = 4;
  Eigen::MatrixXd randomTestFunctionCoefficients;

  double squaredNorm(double u) const;
  double squaredNorm(const BSplines::BSpline& s, double u) const;
  Eigen::MatrixXd numericalSquaredNormDerivatives(double u) const;

  double testEnergy(const Eigen::VectorXd& x) const;
  Eigen::VectorXd testGradients(const Eigen::VectorXd& x) const;
  Eigen::MatrixXd testHessian(const Eigen::VectorXd& x) const;
  Eigen::MatrixXd numericalEnergyDerivatives(double u) const;

  Eigen::MatrixXd numericalDerivativesOfSquaredGradient(double u) const;

  double dotProductOfGradientsAndSplineDerivative(const BSplines::BSpline& s, double u) const;
  Eigen::MatrixXd numericalDerivativesOfDotProductOfGradientAndSplineDerivative(double u) const;

  double dotProductOfPositionAndSplineDerivative(const BSplines::BSpline& s, double uPos, double uTan) const;
  Eigen::MatrixXd numericalDerivativesOfDotProductOfPositionAndSplineDerivative(double uPosition, double uTangent) const;
  double dotProductOfPositions(const BSplines::BSpline& s, double uA, double uB) const;
  Eigen::MatrixXd numericalDerivativesOfDotProductOfPositions(double uA, double uB) const;

  double squareDistance(const BSplines::BSpline& s, double uA, double uB) const;
  Eigen::MatrixXd numericalDerivativesOfSquaredDistance(double uA, double uB) const;

  double squareDistanceToPoint(const BSplines::BSpline& s, double u, const Eigen::VectorXd& p) const;
  Eigen::MatrixXd numericalDerivativesOfSquareDistanceToPoint(double u, const Eigen::VectorXd& p) const;

  double derivativeOfSquareTangentNorm(const BSplines::BSpline& s, double u) const;
  Eigen::MatrixXd numericalDerivativesOfderivativeOfSquareTangentNorm(double u) const;

 protected:
  void SetUp() override {
    const int numberPoints = 20;
    Eigen::MatrixXd randomPoints = Eigen::MatrixXd::Random(numberPoints, numberDimensions);
    InterpolationGenerator generator(randomPoints);
    randomBSpline = generator.generateBSpline();
    randomTestFunctionCoefficients = Eigen::MatrixXd::Random(numberDimensions, numberDimensions);
  }
};

double CommonCostTerms::squaredNorm(double u) const {
  return squaredNorm(randomBSpline, u);
}

double CommonCostTerms::squaredNorm(const BSplines::BSpline& s, double u) const {
  return s.evaluate(u, 1).squaredNorm();
}

Eigen::MatrixXd CommonCostTerms::numericalSquaredNormDerivatives(double u) const {
  int numberControlPoints = randomBSpline.controlPointCount();
  Eigen::MatrixXd derivatives(numberControlPoints, randomBSpline.getDim());

  auto originalNorm = squaredNorm(u);
  double epsilon = 1e-8;
  Eigen::VectorXd knotVector = randomBSpline.getKnotVector();

  for (int i = 0; i < numberControlPoints; ++i) {
    int controlPointIndex = i;
    for (int j = 0; j < randomBSpline.getDim(); ++j) {
      Eigen::MatrixXd controlPointsP = randomBSpline.getControlPointMatrix();
      controlPointsP(controlPointIndex, j) += epsilon;
      BSpline splineP{knotVector, controlPointsP, randomBSpline.getDegree()};
      double newNorm = squaredNorm(splineP, u);
      derivatives(i, j) = (newNorm - originalNorm) / epsilon;
    }
  }

  return derivatives;
}

double CommonCostTerms::testEnergy(const Eigen::VectorXd& x) const {
  // test function: f(x) = x' A x
  return x.transpose() * randomTestFunctionCoefficients * x;
}

Eigen::VectorXd CommonCostTerms::testGradients(const Eigen::VectorXd& x) const {
  // test function: f(x) = x' A x
  //                f'(x) = (A+A') x
  return (randomTestFunctionCoefficients + randomTestFunctionCoefficients.transpose()) * x;
}

Eigen::MatrixXd CommonCostTerms::testHessian(const Eigen::VectorXd& /*x*/) const {
  // test function: f(x) = x' A x
  //                f'(x) = (A+A') x
  //                f''(x) = (A+A')
  return (randomTestFunctionCoefficients + randomTestFunctionCoefficients.transpose());
}

Eigen::MatrixXd CommonCostTerms::numericalEnergyDerivatives(double u) const {
  int numberControlPoints = randomBSpline.controlPointCount();
  Eigen::MatrixXd derivatives(numberControlPoints, randomBSpline.getDim());

  auto originalEnergy = testEnergy(randomBSpline.evaluate(u));
  double epsilon = 1e-8;
  Eigen::VectorXd knotVector = randomBSpline.getKnotVector();

  for (int i = 0; i < numberControlPoints; ++i) {
    int controlPointIndex = i;
    for (int j = 0; j < randomBSpline.getDim(); ++j) {
      Eigen::MatrixXd controlPointsP = randomBSpline.getControlPointMatrix();
      controlPointsP(controlPointIndex, j) += epsilon;
      BSpline splineP{knotVector, controlPointsP, randomBSpline.getDegree()};
      double newEnergy = testEnergy(splineP.evaluate(u));
      derivatives(i, j) = (newEnergy - originalEnergy) / epsilon;
    }
  }

  return derivatives;
}

Eigen::MatrixXd CommonCostTerms::numericalDerivativesOfSquaredGradient(double u) const {
  int numberControlPoints = randomBSpline.controlPointCount();
  Eigen::MatrixXd derivatives(numberControlPoints, randomBSpline.getDim());

  auto originalGradientNorm2 = testGradients(randomBSpline.evaluate(u)).squaredNorm();
  double epsilon = 1e-8;
  Eigen::VectorXd knotVector = randomBSpline.getKnotVector();

  for (int i = 0; i < numberControlPoints; ++i) {
    int controlPointIndex = i;
    for (int j = 0; j < randomBSpline.getDim(); ++j) {
      Eigen::MatrixXd controlPointsP = randomBSpline.getControlPointMatrix();
      controlPointsP(controlPointIndex, j) += epsilon;
      BSpline splineP{knotVector, controlPointsP, randomBSpline.getDegree()};
      double newGradientNorm2 = testGradients(splineP.evaluate(u)).squaredNorm();
      derivatives(i, j) = (newGradientNorm2 - originalGradientNorm2) / epsilon;
    }
  }

  return derivatives;
}

double CommonCostTerms::dotProductOfGradientsAndSplineDerivative(const BSplines::BSpline& s, double u) const {
  auto gradients = testGradients(s.evaluate(u));
  auto splineDerivative = s.evaluate(u, 1);
  return gradients.dot(splineDerivative);
}

Eigen::MatrixXd CommonCostTerms::numericalDerivativesOfDotProductOfGradientAndSplineDerivative(double u) const {
  int numberControlPoints = randomBSpline.controlPointCount();
  Eigen::MatrixXd derivatives(numberControlPoints, randomBSpline.getDim());

  auto originalDotProduct = dotProductOfGradientsAndSplineDerivative(randomBSpline, u);

  double epsilon = 1e-8;
  Eigen::VectorXd knotVector = randomBSpline.getKnotVector();

  for (int i = 0; i < numberControlPoints; ++i) {
    int controlPointIndex = i;
    for (int j = 0; j < randomBSpline.getDim(); ++j) {
      Eigen::MatrixXd controlPointsP = randomBSpline.getControlPointMatrix();
      controlPointsP(controlPointIndex, j) += epsilon;
      BSpline splineP{knotVector, controlPointsP, randomBSpline.getDegree()};
      double newDotProduct = dotProductOfGradientsAndSplineDerivative(splineP, u);
      derivatives(i, j) = (newDotProduct - originalDotProduct) / epsilon;
    }
  }

  return derivatives;
}

double CommonCostTerms::dotProductOfPositionAndSplineDerivative(const BSplines::BSpline& s, double uPos, double uTan) const {
  auto position = s.evaluate(uPos);
  auto splineDerivative = s.evaluate(uTan, 1);
  return position.dot(splineDerivative);
}

Eigen::MatrixXd CommonCostTerms::numericalDerivativesOfDotProductOfPositionAndSplineDerivative(double uPosition,
                                                                                               double uTangent) const {
  int numberControlPoints = randomBSpline.controlPointCount();
  Eigen::MatrixXd derivatives(numberControlPoints, randomBSpline.getDim());

  auto originalDotProduct = dotProductOfPositionAndSplineDerivative(randomBSpline, uPosition, uTangent);

  double epsilon = 1e-8;
  Eigen::VectorXd knotVector = randomBSpline.getKnotVector();

  for (int i = 0; i < numberControlPoints; ++i) {
    int controlPointIndex = i;
    for (int j = 0; j < randomBSpline.getDim(); ++j) {
      Eigen::MatrixXd controlPointsP = randomBSpline.getControlPointMatrix();
      controlPointsP(controlPointIndex, j) += epsilon;
      BSpline splineP{knotVector, controlPointsP, randomBSpline.getDegree()};
      double newDotProduct = dotProductOfPositionAndSplineDerivative(splineP, uPosition, uTangent);
      derivatives(i, j) = (newDotProduct - originalDotProduct) / epsilon;
    }
  }

  return derivatives;
}

double CommonCostTerms::dotProductOfPositions(const BSplines::BSpline& s, double uA, double uB) const {
  auto pA = s.evaluate(uA);
  auto pB = s.evaluate(uB);
  return pA.dot(pB);
}

Eigen::MatrixXd CommonCostTerms::numericalDerivativesOfDotProductOfPositions(double uA, double uB) const {
  int numberControlPoints = randomBSpline.controlPointCount();
  Eigen::MatrixXd derivatives(numberControlPoints, randomBSpline.getDim());

  auto originalDotProduct = dotProductOfPositions(randomBSpline, uA, uB);

  double epsilon = 1e-8;
  Eigen::VectorXd knotVector = randomBSpline.getKnotVector();

  for (int i = 0; i < numberControlPoints; ++i) {
    int controlPointIndex = i;
    for (int j = 0; j < randomBSpline.getDim(); ++j) {
      Eigen::MatrixXd controlPointsP = randomBSpline.getControlPointMatrix();
      controlPointsP(controlPointIndex, j) += epsilon;
      BSpline splineP{knotVector, controlPointsP, randomBSpline.getDegree()};
      double newDotProduct = dotProductOfPositions(splineP, uA, uB);
      derivatives(i, j) = (newDotProduct - originalDotProduct) / epsilon;
    }
  }

  return derivatives;
}

double CommonCostTerms::squareDistance(const BSplines::BSpline& s, double uA, double uB) const {
  auto positionA = s.evaluate(uA);
  auto positionB = s.evaluate(uB);
  return (positionA - positionB).squaredNorm();
}

double CommonCostTerms::squareDistanceToPoint(const BSplines::BSpline& s, double u, const Eigen::VectorXd& p) const {
  auto splinePos = s.evaluate(u);
  return (splinePos - p).squaredNorm();
}

Eigen::MatrixXd CommonCostTerms::numericalDerivativesOfSquaredDistance(double uA, double uB) const {
  int numberControlPoints = randomBSpline.controlPointCount();
  Eigen::MatrixXd derivatives(numberControlPoints, randomBSpline.getDim());

  auto originalSquare = squareDistance(randomBSpline, uA, uB);

  double epsilon = 1e-8;
  Eigen::VectorXd knotVector = randomBSpline.getKnotVector();

  for (int i = 0; i < numberControlPoints; ++i) {
    int controlPointIndex = i;
    for (int j = 0; j < randomBSpline.getDim(); ++j) {
      Eigen::MatrixXd controlPointsP = randomBSpline.getControlPointMatrix();
      controlPointsP(controlPointIndex, j) += epsilon;
      BSpline splineP{knotVector, controlPointsP, randomBSpline.getDegree()};
      double newSquare = squareDistance(splineP, uA, uB);
      derivatives(i, j) = (newSquare - originalSquare) / epsilon;
    }
  }

  return derivatives;
}

Eigen::MatrixXd CommonCostTerms::numericalDerivativesOfSquareDistanceToPoint(double u, const Eigen::VectorXd& p) const {
  int numberControlPoints = randomBSpline.controlPointCount();
  Eigen::MatrixXd derivatives(numberControlPoints, randomBSpline.getDim());

  auto originalSquare = squareDistanceToPoint(randomBSpline, u, p);

  double epsilon = 1e-8;
  Eigen::VectorXd knotVector = randomBSpline.getKnotVector();

  for (int i = 0; i < numberControlPoints; ++i) {
    int controlPointIndex = i;
    for (int j = 0; j < randomBSpline.getDim(); ++j) {
      Eigen::MatrixXd controlPointsP = randomBSpline.getControlPointMatrix();
      controlPointsP(controlPointIndex, j) += epsilon;
      BSpline splineP{knotVector, controlPointsP, randomBSpline.getDegree()};
      double newSquare = squareDistanceToPoint(splineP, u, p);
      derivatives(i, j) = (newSquare - originalSquare) / epsilon;
    }
  }

  return derivatives;
}

double CommonCostTerms::derivativeOfSquareTangentNorm(const BSplines::BSpline& s, double u) const {
  const auto& splineWithDerivative = s;
  Eigen::VectorXd firstDerivative = splineWithDerivative.evaluate(u, 1);
  Eigen::VectorXd secondDerivative = splineWithDerivative.evaluate(u, 2);

  // d(||C'(u)||^2)/du = 2 * dot(C', C'')
  return 2 * firstDerivative.dot(secondDerivative);
}

Eigen::MatrixXd CommonCostTerms::numericalDerivativesOfderivativeOfSquareTangentNorm(double u) const {
  int numberControlPoints = randomBSpline.controlPointCount();
  Eigen::MatrixXd derivatives(numberControlPoints, randomBSpline.getDim());

  auto originalDerivative = derivativeOfSquareTangentNorm(randomBSpline, u);

  double epsilon = 1e-8;
  Eigen::VectorXd knotVector = randomBSpline.getKnotVector();

  for (int i = 0; i < numberControlPoints; ++i) {
    int controlPointIndex = i;
    for (int j = 0; j < randomBSpline.getDim(); ++j) {
      Eigen::MatrixXd controlPointsP = randomBSpline.getControlPointMatrix();
      controlPointsP(controlPointIndex, j) += epsilon;
      BSpline splineP{knotVector, controlPointsP, randomBSpline.getDegree()};
      double newDerivative = derivativeOfSquareTangentNorm(splineP, u);
      derivatives(i, j) = (newDerivative - originalDerivative) / epsilon;
    }
  }

  return derivatives;
}

TEST_F(CommonCostTerms, CalculateTheSquaredDerivativeNormCorrectly) {
  double testU = 0.5;

  auto squaredNormWithDerivatives = CommonTerms::squaredNormOfSplineDerivative(randomBSpline, testU);

  auto expectedSquaredNorm = squaredNorm(testU);
  auto expectedDerivatives = numericalSquaredNormDerivatives(testU);

  ASSERT_THAT(squaredNormWithDerivatives.value(), DoubleNear(expectedSquaredNorm, 1e-10));
  ASSERT_TRUE(squaredNormWithDerivatives.derivatives().isApprox(expectedDerivatives, 1e-5));
}

TEST_F(CommonCostTerms, CalculateTheEnergyCorrectly) {
  double testU = 0.5;

  auto energy = testEnergy(randomBSpline.evaluate(testU));
  auto gradients = testGradients(randomBSpline.evaluate(testU));

  auto energyWithDerivatives = CommonTerms::singlePointEnergy(randomBSpline, testU, energy, gradients);

  auto expectedDerivatives = numericalEnergyDerivatives(testU);

  ASSERT_THAT(energyWithDerivatives.value(), DoubleNear(energy, 1e-10));
  ASSERT_TRUE(energyWithDerivatives.derivatives().isApprox(expectedDerivatives, 1e-5));
}

TEST_F(CommonCostTerms, CalculateTheSquaredNormOfGradientCorrectly) {
  double testU = 0.5;

  auto gradients = testGradients(randomBSpline.evaluate(testU));
  auto hessian = testHessian(randomBSpline.evaluate(testU));
  auto gradientNorm2 = gradients.squaredNorm();

  auto gradientNorm2WithDerivatives = CommonTerms::squaredNormOfGradient(randomBSpline, testU, gradients, hessian);

  auto expectedDerivatives = numericalDerivativesOfSquaredGradient(testU);

  ASSERT_THAT(gradientNorm2WithDerivatives.value(), DoubleNear(gradientNorm2, 1e-10));
  ASSERT_TRUE(gradientNorm2WithDerivatives.derivatives().isApprox(expectedDerivatives, 1e-5));
}

TEST_F(CommonCostTerms, CalculateTheDotProductOfGradientsAndSplineDerivativeCorrectly) {
  double testU = 0.5;

  auto gradients = testGradients(randomBSpline.evaluate(testU));
  auto splineDerivative = randomBSpline.evaluate(testU, 1);
  auto hessian = testHessian(randomBSpline.evaluate(testU));
  auto dotProduct = gradients.dot(splineDerivative);

  auto dotProductWithDerivatives =
      CommonTerms::dotProductOfGradientAndSplineDerivative(randomBSpline, testU, gradients, hessian);

  auto expectedDerivatives = numericalDerivativesOfDotProductOfGradientAndSplineDerivative(testU);

  ASSERT_THAT(dotProductWithDerivatives.value(), DoubleNear(dotProduct, 1e-10));
  ASSERT_TRUE(dotProductWithDerivatives.derivatives().isApprox(expectedDerivatives, 1e-5));
}

TEST_F(CommonCostTerms, CalculateTheDotProductOfPositionAndSplineDerivativeCorrectly) {
  double uPosition = 0.5;
  double uTangent = 0.676;

  auto positions = randomBSpline.evaluate(uPosition);
  auto splineDerivative = randomBSpline.evaluate(uTangent, 1);
  auto dotProduct = positions.dot(splineDerivative);

  auto dotProductWithDerivatives = CommonTerms::dotProductOfPositionAndSplineDerivative(randomBSpline, uPosition, uTangent);

  auto expectedDerivatives = numericalDerivativesOfDotProductOfPositionAndSplineDerivative(uPosition, uTangent);

  ASSERT_THAT(dotProductWithDerivatives.value(), DoubleNear(dotProduct, 1e-10));
  ASSERT_TRUE(dotProductWithDerivatives.derivatives().isApprox(expectedDerivatives, 1e-5));
}

TEST_F(CommonCostTerms, CalculateTheDotProductOfPositionsCorrectly) {
  double uA = 0.5;
  double uB = 0.676;

  auto positionA = randomBSpline.evaluate(uA);
  auto positionB = randomBSpline.evaluate(uB);
  auto dotProduct = positionA.dot(positionB);

  auto dotProductWithDerivatives = CommonTerms::dotProductOfPositions(randomBSpline, uA, uB);

  auto expectedDerivatives = numericalDerivativesOfDotProductOfPositions(uA, uB);

  ASSERT_THAT(dotProductWithDerivatives.value(), DoubleNear(dotProduct, 1e-10));
  ASSERT_TRUE(dotProductWithDerivatives.derivatives().isApprox(expectedDerivatives, 1e-5));
}

TEST_F(CommonCostTerms, CalculateTheSquaredDistanceBetweenTwoPointsCorrectly) {
  double uA = 0.5;
  double uB = 0.676;

  auto positionsA = randomBSpline.evaluate(uA);
  auto positionsB = randomBSpline.evaluate(uB);
  auto distance = (positionsA - positionsB).squaredNorm();

  auto squaredDistanceWithDerivatives = CommonTerms::distanceSquared(randomBSpline, uA, uB);

  auto expectedDerivatives = numericalDerivativesOfSquaredDistance(uA, uB);

  ASSERT_THAT(squaredDistanceWithDerivatives.value(), DoubleNear(distance, 1e-10));
  ASSERT_TRUE(squaredDistanceWithDerivatives.derivatives().isApprox(expectedDerivatives, 1e-5));
}

TEST_F(CommonCostTerms, CalculateTheSquaredDistanceBetweenToAConstantPointCorrectly) {
  double u = 0.5;
  Eigen::VectorXd randomPoint = Eigen::VectorXd::Random(numberDimensions);

  auto splinePosition = randomBSpline.evaluate(u);
  auto distance = (splinePosition - randomPoint).squaredNorm();

  auto squaredDistanceWithDerivatives = CommonTerms::distanceSquaredToPosition(randomBSpline, u, randomPoint);

  auto expectedDerivatives = numericalDerivativesOfSquareDistanceToPoint(u, randomPoint);

  ASSERT_THAT(squaredDistanceWithDerivatives.value(), DoubleNear(distance, 1e-10));
  ASSERT_TRUE(squaredDistanceWithDerivatives.derivatives().isApprox(expectedDerivatives, 1e-5));
}

TEST_F(CommonCostTerms, CalculateTheDerivativeOfSquaredTangentNormCorrectly) {
  double testU = 0.5;

  auto derivOfSquaredTangentNorm = CommonTerms::derivativeOfSquaredTangentNorm(randomBSpline, testU);

  auto expectedValue = derivativeOfSquareTangentNorm(randomBSpline, testU);
  auto expectedDerivatives = numericalDerivativesOfderivativeOfSquareTangentNorm(testU);

  ASSERT_THAT(derivOfSquaredTangentNorm.value(), DoubleNear(expectedValue, 1e-10));
  ASSERT_TRUE(derivOfSquaredTangentNorm.derivatives().isApprox(expectedDerivatives, 1e-5));
}

} // namespace Tests
} // namespace Readuct
} // namespace Scine