/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Readuct/BSplines/BSpline.h>
#include <Readuct/BSplines/DimensionManipulator.h>
#include <Readuct/BSplines/Exceptions.h>
#include <Readuct/BSplines/InterpolationGenerator.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Readuct {
using namespace BSplines;
namespace Tests {

class ADimensionManipulator : public Test {
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

TEST_F(ADimensionManipulator, SplitSplinesHaveCorrectDimensions) {
  int numberDimensions = randomBSpline.getDim();
  int splitIndex = 1;

  auto splitSplines = DimensionManipulator::split(randomBSpline, splitIndex);

  ASSERT_THAT(splitSplines.first.getDim(), Eq(splitIndex));
  ASSERT_THAT(splitSplines.second.getDim(), Eq(numberDimensions - splitIndex));
}

TEST_F(ADimensionManipulator, SplitSplinesHaveSameKnotVector) {
  int splitIndex = 1;

  auto splitSplines = DimensionManipulator::split(randomBSpline, splitIndex);

  ASSERT_THAT(splitSplines.first.getKnotVector(), Eq(randomBSpline.getKnotVector()));
  ASSERT_THAT(splitSplines.second.getKnotVector(), Eq(randomBSpline.getKnotVector()));
}

TEST_F(ADimensionManipulator, CanBeDividedAlongDimensions) {
  auto splitSplines = DimensionManipulator::split(randomBSpline, 1);

  compareDividedSpline(randomBSpline, splitSplines.first, splitSplines.second);
}

TEST_F(ADimensionManipulator, ThrowsForInvalidDimensionIndex) {
  int dims = randomBSpline.getDim();
  ASSERT_THROW(DimensionManipulator::split(randomBSpline, 0), InvalidSplitDimension);
  ASSERT_THROW(DimensionManipulator::split(randomBSpline, dims), InvalidSplitDimension);
}

TEST_F(ADimensionManipulator, MergedSplineHasSameKnotVector) {
  auto mergedSpline = DimensionManipulator::merge(randomBSpline, randomBSpline);

  ASSERT_THAT(mergedSpline.getKnotVector(), Eq(randomBSpline.getKnotVector()));
}

TEST_F(ADimensionManipulator, MergedSplineHasCorrectDimension) {
  auto mergedSpline = DimensionManipulator::merge(randomBSpline, randomBSpline);

  ASSERT_THAT(mergedSpline.getDim(), Eq(2 * randomBSpline.getDim()));
}

TEST_F(ADimensionManipulator, MergeIsWorking) {
  auto mergedSpline = DimensionManipulator::merge(randomBSpline, randomBSpline);

  compareDividedSpline(mergedSpline, randomBSpline, randomBSpline);
}

TEST_F(ADimensionManipulator, ThrowsIfSplinesToMergeDoNotHaveTheSameKnotVector) {
  auto degree = randomBSpline.getDegree();
  Eigen::VectorXd knotVector = randomBSpline.getKnotVector();
  const Eigen::MatrixXd& controlPoints = randomBSpline.getControlPointMatrix();

  // slightly modify one of the knot vector elements
  knotVector(5) *= 0.99;
  BSpline otherSpline{knotVector, controlPoints, degree};

  ASSERT_THROW(DimensionManipulator::merge(randomBSpline, otherSpline), IncompatibleKnotVectorsForDimensionalMerge);
}

TEST_F(ADimensionManipulator, ExtractionTest) {
  // extract the third dimension
  auto extractedSpline = DimensionManipulator::extractDimensions(randomBSpline, 2, 3);

  for (auto u : {0.001111, 0.232323, 0.555, 0.8766, 0.999}) {
    auto extractedSecondDim = extractedSpline.evaluate(u)(0);
    auto expected = randomBSpline.evaluate(u)(2);

    ASSERT_THAT(extractedSecondDim, DoubleEq(expected));
  }
}

} // namespace Tests
} // namespace Readuct
} // namespace Scine