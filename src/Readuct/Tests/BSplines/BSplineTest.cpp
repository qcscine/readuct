/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Readuct/BSplines/BSpline.h>
#include <Readuct/BSplines/Coefficients.h>
#include <Readuct/BSplines/Exceptions.h>
#include <Readuct/BSplines/InterpolationGenerator.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Readuct {
using namespace BSplines;
namespace Tests {

class ABspline : public Test {
 public:
  BSpline randomBSpline;

 protected:
  void SetUp() override {
    const int numberDimensions = 4;
    const int numberPoints = 20;
    Eigen::MatrixXd randomPoints = Eigen::MatrixXd::Random(numberPoints, numberDimensions);
    InterpolationGenerator generator(randomPoints);
    randomBSpline = generator.generateBSpline();
  }
};

namespace {
void compareAtRandomUValues(const BSpline& b1, const BSpline& b2) {
  for (auto u : {0.001111, 0.232323, 0.555, 0.8766, 0.999}) {
    Eigen::VectorXd v1 = b1.evaluate(u);
    Eigen::VectorXd v2 = b2.evaluate(u);
    ASSERT_TRUE(v1.isApprox(v2));
  }
}

void compareDividedSpline(const BSpline& b1, const BSpline& b2, const BSpline& b3) {
  for (auto u : {0.001111, 0.232323, 0.555, 0.8766, 0.999}) {
    Eigen::VectorXd v1 = b1.evaluate(u);
    Eigen::VectorXd v2 = b2.evaluate(u);
    Eigen::VectorXd v3 = b3.evaluate(u);

    Eigen::VectorXd v23{v2.size() + v3.size()};
    v23 << v2, v3;

    ASSERT_TRUE(v1.isApprox(v23));
  }
}
} // namespace

TEST_F(ABspline, CanBeGeneratedFromDegreeAndKnotsAndControlPoints) {
  auto degree = randomBSpline.getDegree();
  Eigen::VectorXd knotVector = randomBSpline.getKnotVector();
  Eigen::MatrixXd controlPoints = randomBSpline.getControlPointMatrix();

  BSpline copiedSpline{knotVector, controlPoints, degree};

  compareAtRandomUValues(randomBSpline, copiedSpline);
}

TEST_F(ABspline, CanBeDividedAlongDimensions) {
  auto degree = randomBSpline.getDegree();
  Eigen::VectorXd knotVector = randomBSpline.getKnotVector();
  Eigen::MatrixXd controlPoints = randomBSpline.getControlPointMatrix();

  BSpline bspline1d{knotVector, controlPoints.leftCols(1), degree};
  BSpline bspline3d{knotVector, controlPoints.rightCols(3), degree};

  compareDividedSpline(randomBSpline, bspline1d, bspline3d);
}

TEST_F(ABspline, CanBeReversed) {
  auto reversed = randomBSpline.reversed();
  for (auto u : {0.002, 0.3, 0.5436, 0.900}) {
    auto v = randomBSpline.evaluate(u);
    auto w = reversed.evaluate(1.0 - u);
    ASSERT_TRUE(v.isApprox(w));
  }
}

TEST_F(ABspline, ReversingTwoTimesDeliversOriginalSpline) {
  auto reversed = randomBSpline.reversed();
  auto reversedTwoTimes = reversed.reversed();
  compareAtRandomUValues(randomBSpline, reversedTwoTimes);
}

TEST_F(ABspline, ReturnsZeroDimensionalValueWhenUninitialized) {
  BSpline uninitialized;

  for (auto u : {0.01, 0.4, 0.8}) {
    auto v = uninitialized.evaluate(u);
    ASSERT_THAT(v.size(), Eq(0));
  }
}

TEST_F(ABspline, CalculatesNoDerivativesByDefault) {
  ASSERT_THAT(randomBSpline.getHighestCalculatedDerivative(), Eq(0));
}

TEST_F(ABspline, CalculatesMissingDerivatives) {
  double u = 0.333;
  ASSERT_THAT(randomBSpline.getHighestCalculatedDerivative(), Eq(0));
  randomBSpline.evaluate(u, 2);
  ASSERT_THAT(randomBSpline.getHighestCalculatedDerivative(), Eq(2));
}

TEST_F(ABspline, EvaluationIsSameAsLinearCombinationOfCoefficientsAndControlPoints) {
  double u = 0.5678;

  const auto& controlPoints = randomBSpline.getControlPointMatrix();
  Eigen::VectorXd bsplineCoefficients = randomBSpline.calculateBSplineCoefficientVector(u);

  Eigen::VectorXd evaluated = randomBSpline.evaluate(u);
  Eigen::VectorXd linearCombination = controlPoints.transpose() * bsplineCoefficients;

  ASSERT_TRUE(evaluated.isApprox(linearCombination));
}

TEST_F(ABspline, CoefficientsClassCorrespondsToCoefficientVector) {
  double u = 0.5678;

  Eigen::VectorXd fullBsplineCoefficients = randomBSpline.calculateBSplineCoefficientVector(u);
  auto coefficientClass = randomBSpline.calculateBSplineCoefficients(u);

  for (int i = 0; i < fullBsplineCoefficients.size(); ++i) {
    ASSERT_THAT(coefficientClass.get(i), DoubleEq(fullBsplineCoefficients(i)));
  }
}

} // namespace Tests
} // namespace Readuct
} // namespace Scine
