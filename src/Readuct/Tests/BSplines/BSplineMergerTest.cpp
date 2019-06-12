/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Readuct/BSplines/BSpline.h>
#include <Readuct/BSplines/ContainerConverter.h>
#include <Readuct/BSplines/ControlPolygonGenerator.h>
#include <Readuct/BSplines/FixedEndsPenalizedLeastSquaresGenerator.h>
#include <Readuct/BSplines/Merger.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Readuct {
using namespace BSplines;
namespace Tests {

class ABSplineMergerTest : public Test {
 public:
  unsigned p;
  BSpline bsMerged;

  void SetUp() override {
    p = 3;
    Eigen::MatrixXd data1(6, 2);
    data1.row(0) = Eigen::Vector2d(-4.0, -4.00);
    data1.row(1) = Eigen::Vector2d(-3.0, -2.00);
    data1.row(2) = Eigen::Vector2d(-2.0, -0.75);
    data1.row(3) = Eigen::Vector2d(-1.0, -0.25);
    data1.row(4) = Eigen::Vector2d(-0.5, -0.05);
    data1.row(5) = Eigen::Vector2d(+0.0, -0.00);

    Eigen::MatrixXd data2(6, 2);
    data2.row(0) = Eigen::Vector2d(+0.0, +0.10);
    data2.row(1) = Eigen::Vector2d(+0.5, +0.05);
    data2.row(2) = Eigen::Vector2d(+1.0, +0.25);
    data2.row(3) = Eigen::Vector2d(+2.0, +0.75);
    data2.row(4) = Eigen::Vector2d(+3.0, +2.00);
    data2.row(5) = Eigen::Vector2d(+4.0, +4.00);

    BSpline bs1 = ControlPolygonGenerator(data1, p, true).generateBSpline();
    BSpline bs2 = ControlPolygonGenerator(data2, p, true).generateBSpline();

    Merger bsMerger_(bs1, bs2);
    bsMerged = bsMerger_.connectToBSpline1();
  }
};

TEST_F(ABSplineMergerTest, CheckConnectToFirst) {
  ASSERT_THAT(bsMerged.evaluate(0.5)(0), DoubleNear(0, 10E-15));
  ASSERT_THAT(bsMerged.evaluate(0.5)(1), DoubleNear(0, 10E-15));
}

} // namespace Tests
} // namespace Readuct
} // namespace Scine