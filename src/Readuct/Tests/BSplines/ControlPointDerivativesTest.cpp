/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Readuct/BSplines/BSpline.h>
#include <Readuct/BSplines/ControlPointDerivatives.h>
#include <Readuct/BSplines/InterpolationGenerator.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Readuct {
using namespace BSplines;
using namespace BSplines::ControlPointDerivatives;
namespace Tests {

class BSplineControlPointDerivatives : public Test {
 public:
  BSpline randomBSpline;
  const int numberDimensions = 4;

 protected:
  void SetUp() override {
    const int numberPoints = 20;
    Eigen::MatrixXd randomPoints = Eigen::MatrixXd::Random(numberPoints, numberDimensions);
    InterpolationGenerator generator(randomPoints);
    randomBSpline = generator.generateBSpline();
  }
};

TEST_F(BSplineControlPointDerivatives, AreCorrectInTheMiddle) {
  auto degree = randomBSpline.getDegree();
  Eigen::VectorXd knotVector = randomBSpline.getKnotVector();

  auto numberControlPoints = randomBSpline.controlPointCount();
  auto controlPointIndex = numberControlPoints / 2;
  auto u = 0.5;

  Eigen::VectorXd values = randomBSpline.evaluate(u);
  Eigen::VectorXd derivative = curveDerivatives(randomBSpline, u, controlPointIndex);

  auto delta = 1e-8;
  for (int i = 0; i < numberDimensions; ++i) {
    Eigen::MatrixXd controlPointsP = randomBSpline.getControlPointMatrix();
    controlPointsP(controlPointIndex, i) += delta;
    BSpline splineP{knotVector, controlPointsP, degree};

    Eigen::VectorXd valuesP = splineP.evaluate(u);

    auto expected = (valuesP[i] - values[i]) / delta;
    ASSERT_THAT(derivative[i], DoubleNear(expected, 1e-6));
  }
}

TEST_F(BSplineControlPointDerivatives, AreZeroIfTheControlPointIsNotLocal) {
  // The third control point plays no role at u = 0.9.
  auto controlPointIndex = 2;
  auto u = 0.9;

  Eigen::VectorXd values = randomBSpline.evaluate(u);
  Eigen::VectorXd derivative = curveDerivatives(randomBSpline, u, controlPointIndex);

  for (int i = 0; i < numberDimensions; ++i) {
    auto expected = 0.0;
    ASSERT_THAT(derivative[i], DoubleNear(expected, 1e-6));
  }
}

TEST_F(BSplineControlPointDerivatives, AreCorrectInTheMiddleForTangent) {
  auto degree = randomBSpline.getDegree();
  Eigen::VectorXd knotVector = randomBSpline.getKnotVector();

  auto numberControlPoints = randomBSpline.controlPointCount();
  auto controlPointIndex = numberControlPoints / 2;
  auto u = 0.5;

  Eigen::VectorXd tangent = randomBSpline.evaluate(u, 1);
  Eigen::VectorXd derivative = firstOrderCurveDerivatives(randomBSpline, u, controlPointIndex);

  auto delta = 1e-8;
  for (int i = 0; i < numberDimensions; ++i) {
    Eigen::MatrixXd controlPointsP = randomBSpline.getControlPointMatrix();
    controlPointsP(controlPointIndex, i) += delta;
    BSpline splineP{knotVector, controlPointsP, degree};

    Eigen::VectorXd tangentP = splineP.evaluate(u, 1);

    auto expected = (tangentP[i] - tangent[i]) / delta;
    ASSERT_THAT(derivative[i], DoubleNear(expected, 1e-5));
  }
}

TEST_F(BSplineControlPointDerivatives, AreZeroForTangentIfTheControlPointIsNotLocal) {
  // The third control point plays no role at u = 0.9.
  auto controlPointIndex = 2;
  auto u = 0.9;

  Eigen::VectorXd values = randomBSpline.evaluate(u);
  Eigen::VectorXd derivative = firstOrderCurveDerivatives(randomBSpline, u, controlPointIndex);

  for (int i = 0; i < numberDimensions; ++i) {
    auto expected = 0.0;
    ASSERT_THAT(derivative[i], DoubleNear(expected, 1e-6));
  }
}

TEST_F(BSplineControlPointDerivatives, AreCorrectInTheMiddleForSecondDerivative) {
  auto degree = randomBSpline.getDegree();
  Eigen::VectorXd knotVector = randomBSpline.getKnotVector();

  auto numberControlPoints = randomBSpline.controlPointCount();
  auto controlPointIndex = numberControlPoints / 2;
  auto u = 0.5;

  Eigen::VectorXd secondDerivative = randomBSpline.evaluate(u, 2);
  Eigen::VectorXd derivative = secondOrderCurveDerivatives(randomBSpline, u, controlPointIndex);

  auto delta = 1e-8;
  for (int i = 0; i < numberDimensions; ++i) {
    Eigen::MatrixXd controlPointsP = randomBSpline.getControlPointMatrix();
    controlPointsP(controlPointIndex, i) += delta;
    BSpline splineP{knotVector, controlPointsP, degree};

    Eigen::VectorXd secondDerivativeP = splineP.evaluate(u, 2);

    auto expected = (secondDerivativeP[i] - secondDerivative[i]) / delta;
    ASSERT_THAT(derivative[i], DoubleNear(expected, 1e-4));
  }
}

TEST_F(BSplineControlPointDerivatives, AreZeroForSecondDerivativeIfTheControlPointIsNotLocal) {
  // The third control point plays no role at u = 0.9.
  auto controlPointIndex = 2;
  auto u = 0.9;

  Eigen::VectorXd values = randomBSpline.evaluate(u);
  Eigen::VectorXd derivative = secondOrderCurveDerivatives(randomBSpline, u, controlPointIndex);

  for (int i = 0; i < numberDimensions; ++i) {
    auto expected = 0.0;
    ASSERT_THAT(derivative[i], DoubleNear(expected, 1e-6));
  }
}

} // namespace Tests
} // namespace Readuct
} // namespace Scine