/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Readuct/BSplines/BSpline.h>
#include <Readuct/BSplines/ContainerConverter.h>
#include <Readuct/BSplines/ControlPolygonGenerator.h>
#include <Readuct/BSplines/RootFinder.h>
#include <Readuct/BSplines/Splitter.h>
#include <gmock/gmock.h>
#include <Eigen/Core>

using namespace testing;
namespace Scine {
namespace Readuct {
using namespace BSplines;
namespace Tests {

/*! Test for the root-finding algorithm.
 * Remaining problem: The root-finding algorithm cannot identify repeated roots that do not originate from p repeated
 * control points being zero but from a special arrangement of control points (see difficult repeated root case tests).
 * In this case, the algorithm finds single roots instead. Accordingly, some repeated roots can be misinterpreted as
 * single roots but such cases can be expected to occur rarely.
 * For the use in BSplineStationaryPointFinder this is no problem because we check for saddle points also for single
 * roots. Therefore the motifs are identified correctly anyways.
 * */
class ABSplineRootFinderTest : public Test {
 public:
  void SetUp() override {
  }
};

TEST_F(ABSplineRootFinderTest, AllCasesDegree1) {
  Eigen::Matrix<double, 88, 1> mat;
  // clang-format off
  mat <<
     1,-1,
    -1, 0, 1, 0,
     1, 1, 0, 0,-1,-1, 0, 0,
    -1,-1,-1, 0, 0, 0, 1, 1, 1, 0, 0, 0,
     1, 1, 1, 1, 0, 0, 0, 0,-1,-1,-1,-1, 0, 0, 0, 0,
    -1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
     1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0,
     1,-1;
  // clang-format on
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 1;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(2);
  refSingleRoots << 5 / 870.0, 865 / 870.0;
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(2);
  refRepeatedRoot << 3 / 87.0, 5 / 87.0;
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(10, 2);
  refRootIntervals << 8 / 87.0, 9 / 87.0, 12 / 87.0, 13 / 87.0, 17 / 87.0, 19 / 87.0, 23 / 87.0, 25 / 87.0, 30 / 87.0,
      33 / 87.0, 38 / 87.0, 41 / 87.0, 47 / 87.0, 51 / 87.0, 57 / 87.0, 61 / 87.0, 68 / 87.0, 73 / 87.0, 80 / 87.0,
      85 / 87.0;
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, AllCasesDegree2) {
  Eigen::Matrix<double, 88, 1> mat;
  // clang-format off
  mat <<
     1,-1,
    -1, 0, 1, 0,
     1, 1, 0, 0,-1,-1, 0, 0,
    -1,-1,-1, 0, 0, 0, 1, 1, 1, 0, 0, 0,
     1, 1, 1, 1, 0, 0, 0, 0,-1,-1,-1,-1, 0, 0, 0, 0,
    -1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
     1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0,
     1,-1;
  // clang-format on
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 2;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  std::cout.precision(17);
  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(3);
  refSingleRoots << 0.0034057351024820047, 5 / 172.0, 0.99639416219251897;
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(2);
  refRepeatedRoot << 4 / 43.0, 6 / 43.0;
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(8, 2);
  refRootIntervals << 17 / 86.0, 9 / 43.0, 23 / 86.0, 12 / 43.0, 15 / 43.0, 16 / 43.0, 19 / 43.0, 20 / 43.0, 47 / 86.0,
      25 / 43.0, 57 / 86.0, 30 / 43.0, 34 / 43.0, 36 / 43.0, 40 / 43.0, 42 / 43.0;
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, AllCasesDegree3) {
  Eigen::Matrix<double, 88, 1> mat;
  // clang-format off
  mat <<
     1,-1,
    -1, 0, 1, 0,
     1, 1, 0, 0,-1,-1, 0, 0,
    -1,-1,-1, 0, 0, 0, 1, 1, 1, 0, 0, 0,
     1, 1, 1, 1, 0, 0, 0, 0,-1,-1,-1,-1, 0, 0, 0, 0,
    -1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
     1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0,
     1,-1;
  // clang-format on
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 3;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }
  std::cout.precision(17);
  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(4);
  refSingleRoots << 0.0024316354319178341, 2 / 85.0, 3 / 34.0, 0.99735926339414371;
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(2);
  refRepeatedRoot << 1 / 5.0, 23 / 85.0;
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(6, 2);
  refRootIntervals << 6 / 17.0, 31 / 85.0, 38 / 85.0, 39 / 85.0, 47 / 85.0, 49 / 85.0, 57 / 85.0, 59 / 85.0, 4 / 5.0,
      71 / 85.0, 16 / 17.0, 83 / 85.0;
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, SingleRootsDegree1) {
  Eigen::Matrix<double, 5, 1> mat;
  mat << 1, -1, 1, -1, 1;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 1;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(4);
  refSingleRoots << 1 / 8.0, 3 / 8.0, 5 / 8.0, 7 / 8.0;
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(0);
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(0, 2);
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, SingleRootsDegree2) {
  Eigen::Matrix<double, 5, 1> mat;
  mat << 1, -1, 1, -1, 1;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 2;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(4);
  refSingleRoots << 1 / 9.0, 1 / 3.0, 2 / 3.0, 8 / 9.0;
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(0);
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(0, 2);
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, SingleRootsDegree3) {
  Eigen::Matrix<double, 5, 1> mat;
  mat << 1, -1, 1, -1, 1;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 3;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(2);
  refSingleRoots << 1 / 8.0, 7 / 8.0;
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(0);
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(0, 2);
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, SingleRootsDegree4) {
  Eigen::Matrix<double, 5, 1> mat;
  mat << 1, -1, -1, -1, 1;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 4;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(2);
  refSingleRoots << 0.1593749806833933, 0.8406250193166066;
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(0);
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(0, 2);
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, repeatedRootDegree1) {
  Eigen::Matrix<double, 3, 1> mat;
  mat << 1, 0, 1;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 1;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen); std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(0);
  EXPECT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(1);
  refRepeatedRoot << 0.5;
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(0, 2);
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, repeatedRootDegree2) {
  Eigen::Matrix<double, 4, 1> mat;
  mat << 1, 0, 0, 1;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 2;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(0);
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(1);
  refRepeatedRoot << 0.5;
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(0, 2);
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, repeatedRootDegree3) {
  Eigen::Matrix<double, 5, 1> mat;
  mat << 1, 0, 0, 0, 1;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 3;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(0);
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(1);
  refRepeatedRoot << 0.5;
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(0, 2);
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, repeatedRootDegree4) {
  Eigen::Matrix<double, 6, 1> mat;
  mat << 1, 0, 0, 0, 0, 1;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 4;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(0);
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(1);
  refRepeatedRoot << 0.5;
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(0, 2);
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, RootIntervalDegree1) {
  Eigen::Matrix<double, 4, 1> mat;
  mat << 1, 0, 0, 1;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 1;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(0);
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(0);
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(1, 2);
  refRootIntervals << 1 / 3.0, 2 / 3.0;
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, RootIntervalDegree2) {
  Eigen::Matrix<double, 5, 1> mat;
  mat << 1, 0, 0, 0, 1;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 2;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(0);
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(0);
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(1, 2);
  refRootIntervals << 1 / 3.0, 2 / 3.0;
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, RootIntervalDegree3) {
  Eigen::Matrix<double, 6, 1> mat;
  mat << 1, 0, 0, 0, 0, 1;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 3;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(0);
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(0);
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(1, 2);
  refRootIntervals << 1 / 3.0, 2 / 3.0;
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, RootIntervalDegree4) {
  Eigen::Matrix<double, 7, 1> mat;
  mat << 1, 0, 0, 0, 0, 0, 1;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 4;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(0);
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(0);
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(1, 2);
  refRootIntervals << 1 / 3.0, 2 / 3.0;
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, SingleRootsAtEndsDegree1) {
  Eigen::Matrix<double, 4, 1> mat;
  mat << 0, 1, -1, 0;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 1;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(3);
  refSingleRoots << 0.0, 0.5, 1.0;
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(0);
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(0, 2);
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, SingleRootsAtEndsDegree2) {
  Eigen::Matrix<double, 4, 1> mat;
  mat << 0, 1, -1, 0;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 2;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(3);
  refSingleRoots << 0.0, 0.5, 1.0;
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(0);
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(0, 2);
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, SingleRootsAtEndsDegree3) {
  Eigen::Matrix<double, 4, 1> mat;
  mat << 0, 1, -1, 0;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 3;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(3);
  refSingleRoots << 0.0, 0.5, 1.0;
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(0);
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(0, 2);
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, SingleRootsAtEndsDegree4) {
  Eigen::Matrix<double, 5, 1> mat;
  mat << 0, 1, -1, 1, 0;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 4;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(2);
  refSingleRoots << 0.0, 1.0;
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(0);
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(0, 2);
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, repeatedRootsAtEndsDegree2) {
  Eigen::Matrix<double, 6, 1> mat;
  mat << 0, 0, 1, -1, 0, 0;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 2;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(3);
  refSingleRoots << 0.0, 0.5, 1.0;
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(0);
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(0, 2);
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, repeatedRootsAtEndsDegree3) {
  Eigen::Matrix<double, 8, 1> mat;
  mat << 0, 0, 0, 1, -1, 0, 0, 0;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 3;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(3);
  refSingleRoots << 0.0, 0.5, 1.0;
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(0);
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(0, 2);
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, repeatedRootsAtEndsDegree4) {
  Eigen::Matrix<double, 10, 1> mat;
  mat << 0, 0, 0, 0, 1, -1, 0, 0, 0, 0;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 4;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(3);
  refSingleRoots << 0.0, 0.5, 1.0;
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(0);
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(0, 2);
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, RootIntervalAtEndsDegree1) {
  Eigen::Matrix<double, 6, 1> mat;
  mat << 0, 0, 1, -1, 0, 0;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 1;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(1);
  refSingleRoots << 0.5;
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(0);
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(2, 2);
  refRootIntervals << 0.0, 1 / 5.0, 4 / 5.0, 1.0;
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, RootIntervalAtEndsDegree2) {
  Eigen::Matrix<double, 8, 1> mat;
  mat << 0, 0, 0, 1, -1, 0, 0, 0;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 2;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(1);
  refSingleRoots << 0.5;
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(0);
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(2, 2);
  refRootIntervals << 0.0, 1 / 6.0, 5 / 6.0, 1.0;
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, RootIntervalAtEndsDegree3) {
  Eigen::Matrix<double, 10, 1> mat;
  mat << 0, 0, 0, 0, 1, -1, 0, 0, 0, 0;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 3;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(1);
  refSingleRoots << 0.5;
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(0);
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(2, 2);
  refRootIntervals << 0.0, 0.14285714285714285, 0.8571428571428571, 1.0;
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, RootIntervalAtEndsDegree4) {
  Eigen::Matrix<double, 12, 1> mat;
  mat << 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 4;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  std::vector<double> repeatedRoots = bsRootFinder.getRepeatedRoots();
  Eigen::VectorXd repeatedRootsEigen;
  repeatedRootsEigen.resize(repeatedRoots.size());
  for (unsigned i = 0; i < repeatedRoots.size(); ++i) {
    repeatedRootsEigen(i) = repeatedRoots[i];
  }

  std::vector<std::pair<double, double>> rootIntervals = bsRootFinder.getRootIntervals();
  Eigen::MatrixXd rootIntervalsEigen;
  rootIntervalsEigen.resize(rootIntervals.size(), 2);
  for (unsigned i = 0; i < rootIntervals.size(); ++i) {
    rootIntervalsEigen(i, 0) = rootIntervals[i].first;
    rootIntervalsEigen(i, 1) = rootIntervals[i].second;
  }

  // ContainerConverter::writeWithOuterBraces(singleRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(repeatedRootsEigen);   std::cout << std::endl;
  // ContainerConverter::writeWithOuterBraces(rootIntervalsEigen); std::cout << std::endl;

  Eigen::VectorXd refSingleRoots;
  refSingleRoots.resize(1);
  refSingleRoots << 0.5;
  ASSERT_TRUE(singleRootsEigen.isApprox(refSingleRoots));

  Eigen::VectorXd refRepeatedRoot;
  refRepeatedRoot.resize(0);
  ASSERT_TRUE(repeatedRootsEigen.isApprox(refRepeatedRoot));

  Eigen::MatrixXd refRootIntervals;
  refRootIntervals.resize(2, 2);
  refRootIntervals << 0.0, 1 / 8.0, 7 / 8.0, 1.0;
  ASSERT_TRUE(rootIntervalsEigen.isApprox(refRootIntervals));
}

TEST_F(ABSplineRootFinderTest, PureRootIntervalDegree1) {
  Eigen::Matrix<double, 6, 1> mat;
  mat << 0, 0, 0, 0, 0, 0;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 1;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  ASSERT_THAT(singleRootsEigen.size(), Eq(0));
}

TEST_F(ABSplineRootFinderTest, PureRootIntervalDegree2) {
  Eigen::Matrix<double, 6, 1> mat;
  mat << 0, 0, 0, 0, 0, 0;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 2;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  ASSERT_THAT(singleRootsEigen.size(), Eq(0));
}

TEST_F(ABSplineRootFinderTest, PureRootIntervalDegree3) {
  Eigen::Matrix<double, 6, 1> mat;
  mat << 0, 0, 0, 0, 0, 0;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 3;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  ASSERT_THAT(singleRootsEigen.size(), Eq(0));
}

TEST_F(ABSplineRootFinderTest, PureRootIntervalDegree4) {
  Eigen::Matrix<double, 6, 1> mat;
  mat << 0, 0, 0, 0, 0, 0;
  Eigen::MatrixXd data;
  data = mat;

  unsigned p = 4;
  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  BSpline bs = myBSplineGenerator_.generateBSpline();

  RootFinder bsRootFinder(bs);
  bsRootFinder.findRoots(0, 1e-13);

  std::vector<double> singleRoots = bsRootFinder.getSingleRoots();
  Eigen::VectorXd singleRootsEigen;
  singleRootsEigen.resize(singleRoots.size());
  for (unsigned i = 0; i < singleRoots.size(); ++i) {
    singleRootsEigen(i) = singleRoots[i];
  }

  ASSERT_THAT(singleRootsEigen.size(), Eq(0));
}

} // namespace Tests
} // namespace Readuct
} // namespace Scine