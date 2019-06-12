/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Readuct/BSplines/BSpline.h>
#include <Readuct/BSplines/ContainerConverter.h>
#include <Readuct/BSplines/Splitter.h>
#include <gmock/gmock.h>
#include <Eigen/Core>

using namespace testing;
namespace Scine {
namespace Readuct {
using namespace BSplines;
namespace Tests {

class ABSplineSplitterTest : public Test {
 public:
  unsigned p;
  BSpline bs, bsLeft, bsRight;
  Splitter bsSplitter;

  void SetUp() override {
  }
};

TEST_F(ABSplineSplitterTest, DisconnectednessDegree3) {
  unsigned p = 3;
  Eigen::MatrixXd controlPoints(5, 1);
  controlPoints << 1, 2, 3, 4, 5;

  Eigen::MatrixXd knotVector(9, 1);
  knotVector << 0, 0, 0, 0, 0.5, 1, 1, 1, 1;

  BSpline bs = BSpline(knotVector, controlPoints, p);

  double uSplit = 0.5;
  auto splitResult = bsSplitter.split(uSplit, bs, {false, false});

  bsLeft = splitResult.first;
  bsRight = splitResult.second;

  Eigen::VectorXd knotVectorLeft = bsLeft.getKnotVector();
  Eigen::MatrixXd controlPointsLeft = bsLeft.getControlPointMatrix();

  /*
  std::cout << "{UpLeft,PpLeft}={";
  ContainerConverter::writeWithOuterBraces(knotVectorLeft);
  std::cout <<"," << std::endl;
  ContainerConverter::writeWithOuterBraces(controlPointsLeft);
  std::cout <<"};" << std::endl << std::endl;
  */

  Eigen::VectorXd knotVectorRight = bsRight.getKnotVector();
  Eigen::MatrixXd controlPointsRight = bsRight.getControlPointMatrix();

  /*
  std::cout << "{UpRight,PpRight}={";
  ContainerConverter::writeWithOuterBraces(knotVectorRight);
  std::cout <<"," << std::endl;
  ContainerConverter::writeWithOuterBraces(controlPointsRight);
  std::cout <<"};" << std::endl << std::endl;
  */

  Eigen::VectorXd refKnotVectorLeft;
  refKnotVectorLeft.resize(8);
  refKnotVectorLeft << 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5;
  ASSERT_TRUE(knotVectorLeft.isApprox(refKnotVectorLeft));

  Eigen::VectorXd refControlPointsLeft;
  refControlPointsLeft.resize(4);
  refControlPointsLeft << 1, 2, 2.5, 3;
  ASSERT_TRUE(controlPointsLeft.isApprox(refControlPointsLeft));

  Eigen::VectorXd refKnotVectorRight;
  refKnotVectorRight.resize(8);
  refKnotVectorRight << 0.5, 0.5, 0.5, 0.5, 1, 1, 1, 1;
  ASSERT_TRUE(knotVectorRight.isApprox(refKnotVectorRight));

  Eigen::VectorXd refControlPointsRight;
  refControlPointsRight.resize(4);
  refControlPointsRight << 3, 3.5, 4, 5;
  ASSERT_TRUE(controlPointsRight.isApprox(refControlPointsRight));
}

TEST_F(ABSplineSplitterTest, Multiplicity0Degree3) {
  unsigned p = 3;
  Eigen::MatrixXd controlPoints(5, 1);
  controlPoints << 1, 2, 3, 4, 5;

  Eigen::MatrixXd knotVector(9, 1);
  knotVector << 0, 0, 0, 0, 0.5, 1, 1, 1, 1;

  BSpline bs = BSpline(knotVector, controlPoints, p);

  double uSplit = 0.75;
  auto splitResult = bsSplitter.split(uSplit, bs, {false, false});

  bsLeft = splitResult.first;
  bsRight = splitResult.second;

  Eigen::VectorXd knotVectorLeft = bsLeft.getKnotVector();
  Eigen::MatrixXd controlPointsLeft = bsLeft.getControlPointMatrix();

  /*
  std::cout << "{UpLeft,PpLeft}={";
  ContainerConverter::writeWithOuterBraces(knotVectorLeft);
  std::cout <<"," << std::endl;
  ContainerConverter::writeWithOuterBraces(controlPointsLeft);
  std::cout <<"};" << std::endl << std::endl;
  */

  Eigen::VectorXd knotVectorRight = bsRight.getKnotVector();
  Eigen::MatrixXd controlPointsRight = bsRight.getControlPointMatrix();

  /*
  std::cout << "{UpRight,PpRight}={";
  ContainerConverter::writeWithOuterBraces(knotVectorRight);
  std::cout <<"," << std::endl;
  ContainerConverter::writeWithOuterBraces(controlPointsRight);
  std::cout <<"};" << std::endl << std::endl;
  */

  Eigen::VectorXd refKnotVectorLeft;
  refKnotVectorLeft.resize(9);
  refKnotVectorLeft << 0, 0, 0, 0, 0.5, 0.75, 0.75, 0.75, 0.75;
  ASSERT_TRUE(knotVectorLeft.isApprox(refKnotVectorLeft));

  Eigen::VectorXd refControlPointsLeft;
  refControlPointsLeft.resize(5);
  refControlPointsLeft << 1, 2, 2.75, 3.5, 3.8125;
  ASSERT_TRUE(controlPointsLeft.isApprox(refControlPointsLeft));

  Eigen::VectorXd refKnotVectorRight;
  refKnotVectorRight.resize(8);
  refKnotVectorRight << 0.75, 0.75, 0.75, 0.75, 1, 1, 1, 1;
  ASSERT_TRUE(knotVectorRight.isApprox(refKnotVectorRight));

  Eigen::VectorXd refControlPointsRight;
  refControlPointsRight.resize(4);
  refControlPointsRight << 3.8125, 4.125, 4.5, 5;
  ASSERT_TRUE(controlPointsRight.isApprox(refControlPointsRight));
}

TEST_F(ABSplineSplitterTest, DisconnectednessDegree1) {
  unsigned p = 1;
  Eigen::MatrixXd controlPoints(5, 1);
  controlPoints << 1, 2, 3, 4, 5;

  Eigen::MatrixXd knotVector(7, 1);
  knotVector << 0, 0, 0.25, 0.5, 0.75, 1, 1;

  BSpline bs = BSpline(knotVector, controlPoints, p);

  double uSplit = 0.5;
  auto splitResult = bsSplitter.split(uSplit, bs, {false, false});

  bsLeft = splitResult.first;
  bsRight = splitResult.second;

  Eigen::VectorXd knotVectorLeft = bsLeft.getKnotVector();
  Eigen::MatrixXd controlPointsLeft = bsLeft.getControlPointMatrix();

  /*
  std::cout << "{UpLeft,PpLeft}={";
  ContainerConverter::writeWithOuterBraces(knotVectorLeft);
  std::cout <<"," << std::endl;
  ContainerConverter::writeWithOuterBraces(controlPointsLeft);
  std::cout <<"};" << std::endl << std::endl;
  */

  Eigen::VectorXd knotVectorRight = bsRight.getKnotVector();
  Eigen::MatrixXd controlPointsRight = bsRight.getControlPointMatrix();

  /*
  std::cout << "{UpRight,PpRight}={";
  ContainerConverter::writeWithOuterBraces(knotVectorRight);
  std::cout <<"," << std::endl;
  ContainerConverter::writeWithOuterBraces(controlPointsRight);
  std::cout <<"};" << std::endl << std::endl;
  */

  Eigen::VectorXd refKnotVectorLeft;
  refKnotVectorLeft.resize(5);
  refKnotVectorLeft << 0, 0, 0.25, 0.5, 0.5;
  ASSERT_TRUE(knotVectorLeft.isApprox(refKnotVectorLeft));

  Eigen::VectorXd refControlPointsLeft;
  refControlPointsLeft.resize(3);
  refControlPointsLeft << 1, 2, 3;
  ASSERT_TRUE(controlPointsLeft.isApprox(refControlPointsLeft));

  Eigen::VectorXd refKnotVectorRight;
  refKnotVectorRight.resize(5);
  refKnotVectorRight << 0.5, 0.5, 0.75, 1, 1;
  ASSERT_TRUE(knotVectorRight.isApprox(refKnotVectorRight));

  Eigen::VectorXd refControlPointsRight;
  refControlPointsRight.resize(3);
  refControlPointsRight << 3, 4, 5;
  ASSERT_TRUE(controlPointsRight.isApprox(refControlPointsRight));
}

} // namespace Tests
} // namespace Readuct
} // namespace Scine