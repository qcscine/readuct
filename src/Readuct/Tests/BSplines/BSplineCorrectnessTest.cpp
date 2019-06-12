/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Readuct/BSplines/BSpline.h>
#include <Readuct/BSplines/ContainerConverter.h>
#include <Readuct/BSplines/ControlPolygonGenerator.h>
#include <gmock/gmock.h>
#include <Eigen/Core>

using namespace testing;
namespace Scine {
namespace Readuct {
using namespace BSplines;
namespace Tests {

class ABSplineCorrectnessTest : public Test {
 public:
  unsigned p;
  unsigned dim;
  BSpline bs;

  void SetUp() override {
    p = 3;
    Eigen::MatrixXd data(5, 2);
    data.row(0) = Eigen::Vector2d(1, 1);
    data.row(1) = Eigen::Vector2d(1.5, -0.5);
    data.row(2) = Eigen::Vector2d(3, 1);
    data.row(3) = Eigen::Vector2d(4.5, 0.5);
    data.row(4) = Eigen::Vector2d(5, -1);

    ControlPolygonGenerator bsGenerator(data, p, true);
    bs = bsGenerator.generateBSpline();
  }
};

TEST_F(ABSplineCorrectnessTest, CorrectDerivatives0) {
  Eigen::Vector2d ref(2);

  ref = Eigen::Vector2d(1., 1.);
  ASSERT_TRUE(bs.evaluate(0.0, 0).isApprox(ref));
  ref = Eigen::Vector2d(2.552, 0.312);
  ASSERT_TRUE(bs.evaluate(0.4, 0).isApprox(ref));
  ref = Eigen::Vector2d(3.448, 0.584);
  ASSERT_TRUE(bs.evaluate(0.6, 0).isApprox(ref));
  ref = Eigen::Vector2d(5., -1.);
  ASSERT_TRUE(bs.evaluate(1.0, 0).isApprox(ref));
}

TEST_F(ABSplineCorrectnessTest, CorrectDerivatives1) {
  Eigen::Vector2d ref(2);

  ref = Eigen::Vector2d(3., -9.);
  ASSERT_TRUE(bs.evaluate(0.0, 1).isApprox(ref));
  ref = Eigen::Vector2d(4.44, 2.04);
  ASSERT_TRUE(bs.evaluate(0.4, 1).isApprox(ref));
  ref = Eigen::Vector2d(4.44, 0.12);
  ASSERT_TRUE(bs.evaluate(0.6, 1).isApprox(ref));
  ref = Eigen::Vector2d(3., -9.);
  ASSERT_TRUE(bs.evaluate(1.0, 1).isApprox(ref));
}

TEST_F(ABSplineCorrectnessTest, CorrectDerivatives2) {
  Eigen::Vector2d ref(2);

  ref = Eigen::Vector2d(6., 54.);
  ASSERT_TRUE(bs.evaluate(0.0, 2).isApprox(ref));
  ref = Eigen::Vector2d(1.2, 1.2);
  ASSERT_TRUE(bs.evaluate(0.4, 2).isApprox(ref));
  ref = Eigen::Vector2d(-1.2, -15.6);
  ASSERT_TRUE(bs.evaluate(0.6, 2).isApprox(ref));
  ref = Eigen::Vector2d(-6., -30.);
  ASSERT_TRUE(bs.evaluate(1.0, 2).isApprox(ref));
}

TEST_F(ABSplineCorrectnessTest, CorrectDerivatives3) {
  Eigen::Vector2d ref(2);

  ref = Eigen::Vector2d(-12., -132.);
  ASSERT_TRUE(bs.evaluate(0.0, 3).isApprox(ref));
  ref = Eigen::Vector2d(-12., -132.);
  ASSERT_TRUE(bs.evaluate(0.4, 3).isApprox(ref));
  ref = Eigen::Vector2d(-12., -36.);
  ASSERT_TRUE(bs.evaluate(0.6, 3).isApprox(ref));
  ref = Eigen::Vector2d(-12., -36.);
  ASSERT_TRUE(bs.evaluate(1.0, 3).isApprox(ref));
}

TEST_F(ABSplineCorrectnessTest, CorrectDerivatives4) {
  Eigen::Vector2d ref(2);

  ref = Eigen::Vector2d(0., 0.);
  ASSERT_TRUE(bs.evaluate(0.0, 4).isApprox(ref));
  ref = Eigen::Vector2d(0., 0.);
  ASSERT_TRUE(bs.evaluate(0.4, 4).isApprox(ref));
  ref = Eigen::Vector2d(0., 0.);
  ASSERT_TRUE(bs.evaluate(0.6, 4).isApprox(ref));
  ref = Eigen::Vector2d(0., 0.);
  ASSERT_TRUE(bs.evaluate(1.0, 4).isApprox(ref));
}

TEST_F(ABSplineCorrectnessTest, CorrectDerivatives5) {
  Eigen::Vector2d ref(2);

  ref = Eigen::Vector2d(0., 0.);
  ASSERT_TRUE(bs.evaluate(0.0, 5).isApprox(ref));
  ref = Eigen::Vector2d(0., 0.);
  ASSERT_TRUE(bs.evaluate(0.4, 5).isApprox(ref));
  ref = Eigen::Vector2d(0., 0.);
  ASSERT_TRUE(bs.evaluate(0.6, 5).isApprox(ref));
  ref = Eigen::Vector2d(0., 0.);
  ASSERT_TRUE(bs.evaluate(1.0, 5).isApprox(ref));
}

TEST_F(ABSplineCorrectnessTest, CompareEvaluationMethods) {
  for (unsigned k = 0; k < p; ++k) {
    ASSERT_TRUE(bs.evaluateNaive(0.0, k).isApprox(bs.evaluate(0.0, k)));
  }
  for (unsigned k = 0; k < p; ++k) {
    ASSERT_TRUE(bs.evaluateNaive(1 / 3.0, k).isApprox(bs.evaluate(1 / 3.0, k)));
  }
  for (unsigned k = 0; k < p; ++k) {
    ASSERT_TRUE(bs.evaluateNaive(0.5, k).isApprox(bs.evaluate(0.5, k)));
  }
  for (unsigned k = 0; k < p; ++k) {
    ASSERT_TRUE(bs.evaluateNaive(0.6, k).isApprox(bs.evaluate(0.6, k)));
  }
  for (unsigned k = 0; k < p; ++k) {
    ASSERT_TRUE(bs.evaluateNaive(1.0, k).isApprox(bs.evaluate(1.0, k)));
  }
}

} // namespace Tests
} // namespace Readuct
} // namespace Scine