/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Readuct/BSplines/Exceptions.h>
#include <Readuct/BSplines/FixedEndsFixedKnotSplineGenerator.h>
#include <Readuct/BSplines/LooseEndsFixedKnotSplineGenerator.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Readuct {
using namespace BSplines;
namespace Tests {

class AFixedKnotSplineGenerator : public Test {
 public:
  int splineDegree = 3;
  Eigen::VectorXd knotVector;
  Eigen::VectorXd dataPoints;
  Eigen::VectorXd dataCoords;

  int numberDataPoints = 12;

 protected:
  void SetUp() override {
    knotVector.resize(10); // 10 = nControlPoints + degree + 1 --> 6 control points
    knotVector << 0, 0, 0, 0, 0.33, 0.66, 1, 1, 1, 1;

    dataCoords.resize(numberDataPoints);
    dataPoints.resize(numberDataPoints);
    dataCoords << 0.0, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0;
    dataPoints << 9.9, 9.9, 9.80, 9.6, 8.2, 8.3, 9.0, 9.1, 8.0, 5.0, 4.0, 1.0;
  }
};

TEST_F(AFixedKnotSplineGenerator, LooseEndsDoesNotModifyTheKnotVector) {
  LooseEndsFixedKnotSplineGenerator generator(dataPoints, dataCoords, knotVector, splineDegree);
  auto spline = generator.generateBSpline();
  ASSERT_THAT(spline.getKnotVector(), Eq(knotVector));
}

TEST_F(AFixedKnotSplineGenerator, LooseEndsDeliversApproximatelyCorrectCurve) {
  LooseEndsFixedKnotSplineGenerator generator(dataPoints, dataCoords, knotVector, splineDegree);
  auto spline = generator.generateBSpline();

  for (int i = 0; i < numberDataPoints; ++i) {
    double u = dataCoords[i];
    double expected = dataPoints[i];
    double obtained = spline.evaluate(u)(0);
    // value must be within +- 1.0
    ASSERT_THAT(obtained, DoubleNear(expected, 1.0));
  }
}

TEST_F(AFixedKnotSplineGenerator, LooseEndsIsSymmetric) {
  // Here we take the opposite of every point and add ten.
  // Then we look if this is also the case for the spline
  LooseEndsFixedKnotSplineGenerator generator(dataPoints, dataCoords, knotVector, splineDegree);
  auto spline = generator.generateBSpline();

  for (int i = 0; i < numberDataPoints; ++i) {
    dataPoints[i] = 10 - dataPoints[i];
  }

  LooseEndsFixedKnotSplineGenerator generator2(dataPoints, dataCoords, knotVector, splineDegree);
  auto spline2 = generator2.generateBSpline();

  for (int i = 0; i < numberDataPoints; ++i) {
    double u = dataCoords[i];
    double sum = spline.evaluate(u)(0) + spline2.evaluate(u)(0);
    ASSERT_THAT(sum, DoubleNear(10.0, 1e-8));
  }
}

TEST_F(AFixedKnotSplineGenerator, FixedEndsDoesNotModifyTheKnotVector) {
  FixedEndsFixedKnotSplineGenerator generator(dataPoints, dataCoords, knotVector, splineDegree);
  auto spline = generator.generateBSpline();
  ASSERT_THAT(spline.getKnotVector(), Eq(knotVector));
}

TEST_F(AFixedKnotSplineGenerator, FixedEndsDoesNotModifyTheFirstAndLastPoints) {
  FixedEndsFixedKnotSplineGenerator generator(dataPoints, dataCoords, knotVector, splineDegree);
  auto spline = generator.generateBSpline();
  auto lastIndex = dataPoints.size() - 1;
  ASSERT_THAT(spline.evaluate(0)(0), DoubleEq(dataPoints(0)));
  ASSERT_THAT(spline.evaluate(1.0)(0), DoubleEq(dataPoints(lastIndex)));
}

TEST_F(AFixedKnotSplineGenerator, FixedEndsDeliversApproximatelyCorrectCurve) {
  FixedEndsFixedKnotSplineGenerator generator(dataPoints, dataCoords, knotVector, splineDegree);
  auto spline = generator.generateBSpline();

  for (int i = 0; i < numberDataPoints; ++i) {
    double u = dataCoords[i];
    double expected = dataPoints[i];
    double obtained = spline.evaluate(u)(0);
    // value must be within +- 1.0
    ASSERT_THAT(obtained, DoubleNear(expected, 1.0));
  }
}

TEST_F(AFixedKnotSplineGenerator, FixedEndsIsSymmetric) {
  // Here we take the opposite of every point and add ten.
  // Then we look if this is also the case for the spline
  FixedEndsFixedKnotSplineGenerator generator(dataPoints, dataCoords, knotVector, splineDegree);
  auto spline = generator.generateBSpline();

  for (int i = 0; i < numberDataPoints; ++i) {
    dataPoints[i] = 10 - dataPoints[i];
  }

  FixedEndsFixedKnotSplineGenerator generator2(dataPoints, dataCoords, knotVector, splineDegree);
  auto spline2 = generator2.generateBSpline();

  for (int i = 0; i < numberDataPoints; ++i) {
    double u = dataCoords[i];
    double sum = spline.evaluate(u)(0) + spline2.evaluate(u)(0);
    ASSERT_THAT(sum, DoubleNear(10.0, 1e-8));
  }
}

TEST_F(AFixedKnotSplineGenerator, FixedEndsThrowsIfFirstCoordinateIsNotZero) {
  // The first data point must have coordinate 0.0 and the last one must have coordinate 1.0
  dataCoords(0) = 0.001;
  ASSERT_THROW(FixedEndsFixedKnotSplineGenerator generator(dataPoints, dataCoords, knotVector, splineDegree), Exception);
}

TEST_F(AFixedKnotSplineGenerator, FixedEndsThrowsIfLastCoordinateIsNotOne) {
  // The first data point must have coordinate 0.0 and the last one must have coordinate 1.0
  dataCoords(numberDataPoints - 1) = 0.9999;
  ASSERT_THROW(FixedEndsFixedKnotSplineGenerator generator(dataPoints, dataCoords, knotVector, splineDegree), Exception);
}

} // namespace Tests
} // namespace Readuct
} // namespace Scine