/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Readuct/BSplines/BSpline.h>
#include <Readuct/BSplines/BSplineTools.h>
#include <Readuct/BSplines/ControlPolygonGenerator.h>
#include <gmock/gmock.h>
#include <Eigen/Core>

using namespace testing;
namespace Scine {
namespace Readuct {
using namespace BSplines;
namespace Tests {

class ABSplineToolsTest : public Test {
 public:
  unsigned p;
  BSpline bs;

  unsigned l;
  unsigned idx;

  void SetUp() override {
  }
};

TEST_F(ABSplineToolsTest, LeftEqualCorrectIdDegree1) {
  p = 1;

  Eigen::MatrixXd data(3, 2);
  data.row(0) = Eigen::Vector2d(1, 1);
  data.row(1) = Eigen::Vector2d(2, 0);
  data.row(2) = Eigen::Vector2d(3, 1);

  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  bs = myBSplineGenerator_.generateBSpline();

  // U =  { 0 , 0 ,0.5, 1 , 1 };
  //      | 0 | 1 | 2 | 3 | 4 |

  l = BSplineTools::findIdxOfLeftOrEqualDomainKnot(0.0, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(1u));

  l = BSplineTools::findIdxOfLeftOrEqualDomainKnot(0.4999, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(1u));

  l = BSplineTools::findIdxOfLeftOrEqualDomainKnot(0.5, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(2u));

  l = BSplineTools::findIdxOfLeftOrEqualDomainKnot(0.5001, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(2u));

  l = BSplineTools::findIdxOfLeftOrEqualDomainKnot(1.0, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(3u));
}

TEST_F(ABSplineToolsTest, LeftCorrectIdDegree1) {
  p = 1;

  Eigen::MatrixXd data(3, 2);
  data.row(0) = Eigen::Vector2d(1, 1);
  data.row(1) = Eigen::Vector2d(2, 0);
  data.row(2) = Eigen::Vector2d(3, 1);

  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  bs = myBSplineGenerator_.generateBSpline();

  // U =  { 0 , 0 ,0.5, 1 , 1 };
  //      | 0 | 1 | 2 | 3 | 4 |

  l = BSplineTools::findIdxOfLeftDomainKnot(0.0, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(1u));

  l = BSplineTools::findIdxOfLeftDomainKnot(0.4999, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(1u));

  l = BSplineTools::findIdxOfLeftDomainKnot(0.5, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(1u));

  l = BSplineTools::findIdxOfLeftDomainKnot(0.5001, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(2u));

  l = BSplineTools::findIdxOfLeftDomainKnot(1.0, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(2u));
}

TEST_F(ABSplineToolsTest, RightCorrectIdDegree1) {
  p = 1;

  Eigen::MatrixXd data(3, 2);
  data.row(0) = Eigen::Vector2d(1, 1);
  data.row(1) = Eigen::Vector2d(2, 0);
  data.row(2) = Eigen::Vector2d(3, 1);

  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  bs = myBSplineGenerator_.generateBSpline();

  // U =  { 0 , 0 ,0.5, 1 , 1 };
  //      | 0 | 1 | 2 | 3 | 4 |

  l = BSplineTools::findIdxOfRightDomainKnot(0.0, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(2u));

  l = BSplineTools::findIdxOfRightDomainKnot(0.4999, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(2u));

  l = BSplineTools::findIdxOfRightDomainKnot(0.5, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(3u));

  l = BSplineTools::findIdxOfRightDomainKnot(0.5001, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(3u));

  l = BSplineTools::findIdxOfRightDomainKnot(1.0, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(3u));
}

TEST_F(ABSplineToolsTest, RightOrEqualCorrectIdDegree1) {
  p = 1;

  Eigen::MatrixXd data(3, 2);
  data.row(0) = Eigen::Vector2d(1, 1);
  data.row(1) = Eigen::Vector2d(2, 0);
  data.row(2) = Eigen::Vector2d(3, 1);

  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  bs = myBSplineGenerator_.generateBSpline();

  // U =  { 0 , 0 ,0.5, 1 , 1 };
  //      | 0 | 1 | 2 | 3 | 4 |

  l = BSplineTools::findIdxOfRightOrEqualDomainKnot(0.0, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(1u));

  l = BSplineTools::findIdxOfRightOrEqualDomainKnot(0.4999, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(2u));

  l = BSplineTools::findIdxOfRightOrEqualDomainKnot(0.5, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(2u));

  l = BSplineTools::findIdxOfRightOrEqualDomainKnot(0.5001, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(3u));

  l = BSplineTools::findIdxOfRightOrEqualDomainKnot(1.0, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(3u));
}

TEST_F(ABSplineToolsTest, LeftEqualCorrectIdDegree3) {
  p = 3;

  Eigen::MatrixXd data(5, 2);
  data.row(0) = Eigen::Vector2d(1, 1);
  data.row(1) = Eigen::Vector2d(2, 0);
  data.row(2) = Eigen::Vector2d(3, 1);
  data.row(3) = Eigen::Vector2d(4, 0);
  data.row(4) = Eigen::Vector2d(5, 1);

  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  bs = myBSplineGenerator_.generateBSpline();

  // U =  { 0 , 0 , 0 , 0 ,0.5, 1 , 1 , 1 , 1 };
  //      | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |

  l = BSplineTools::findIdxOfLeftOrEqualDomainKnot(0.0, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(3u));

  l = BSplineTools::findIdxOfLeftOrEqualDomainKnot(0.4999, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(3u));

  l = BSplineTools::findIdxOfLeftOrEqualDomainKnot(0.5, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(4u));

  l = BSplineTools::findIdxOfLeftOrEqualDomainKnot(0.5001, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(4u));

  l = BSplineTools::findIdxOfLeftOrEqualDomainKnot(1.0, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(5u));
}

TEST_F(ABSplineToolsTest, LeftCorrectIdDegree3) {
  p = 3;

  Eigen::MatrixXd data(5, 2);
  data.row(0) = Eigen::Vector2d(1, 1);
  data.row(1) = Eigen::Vector2d(2, 0);
  data.row(2) = Eigen::Vector2d(3, 1);
  data.row(3) = Eigen::Vector2d(4, 0);
  data.row(4) = Eigen::Vector2d(5, 1);

  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  bs = myBSplineGenerator_.generateBSpline();

  // U =  { 0 , 0 , 0 , 0 ,0.5, 1 , 1 , 1 , 1 };
  //      | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |

  l = BSplineTools::findIdxOfLeftDomainKnot(0.0, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(3u));

  l = BSplineTools::findIdxOfLeftDomainKnot(0.4999, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(3u));

  l = BSplineTools::findIdxOfLeftDomainKnot(0.5, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(3u));

  l = BSplineTools::findIdxOfLeftDomainKnot(0.5001, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(4u));

  l = BSplineTools::findIdxOfLeftDomainKnot(1.0, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(4u));
}

TEST_F(ABSplineToolsTest, RightCorrectIdDegree3) {
  p = 3;

  Eigen::MatrixXd data(5, 2);
  data.row(0) = Eigen::Vector2d(1, 1);
  data.row(1) = Eigen::Vector2d(2, 0);
  data.row(2) = Eigen::Vector2d(3, 1);
  data.row(3) = Eigen::Vector2d(4, 0);
  data.row(4) = Eigen::Vector2d(5, 1);

  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  bs = myBSplineGenerator_.generateBSpline();

  // U =  { 0 , 0 , 0 , 0 ,0.5, 1 , 1 , 1 , 1 };
  //      | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |

  l = BSplineTools::findIdxOfRightDomainKnot(0.0, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(4u));

  l = BSplineTools::findIdxOfRightDomainKnot(0.4999, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(4u));

  l = BSplineTools::findIdxOfRightDomainKnot(0.5, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(5u));

  l = BSplineTools::findIdxOfRightDomainKnot(0.5001, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(5u));

  l = BSplineTools::findIdxOfRightDomainKnot(1.0, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(5u));
}

TEST_F(ABSplineToolsTest, RightOrEqualCorrectIdDegree3) {
  p = 3;

  Eigen::MatrixXd data(5, 2);
  data.row(0) = Eigen::Vector2d(1, 1);
  data.row(1) = Eigen::Vector2d(2, 0);
  data.row(2) = Eigen::Vector2d(3, 1);
  data.row(3) = Eigen::Vector2d(4, 0);
  data.row(4) = Eigen::Vector2d(5, 1);

  ControlPolygonGenerator myBSplineGenerator_(data, p, true);
  bs = myBSplineGenerator_.generateBSpline();

  // U =  { 0 , 0 , 0 , 0 ,0.5, 1 , 1 , 1 , 1 };
  //      | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |

  l = BSplineTools::findIdxOfRightOrEqualDomainKnot(0.0, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(3u));

  l = BSplineTools::findIdxOfRightOrEqualDomainKnot(0.4999, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(4u));

  l = BSplineTools::findIdxOfRightOrEqualDomainKnot(0.5, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(4u));

  l = BSplineTools::findIdxOfRightOrEqualDomainKnot(0.5001, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(5u));

  l = BSplineTools::findIdxOfRightOrEqualDomainKnot(1.0, p, bs.getKnotVector());
  ASSERT_THAT(l, Eq(5u));
}

} // namespace Tests
} // namespace Readuct
} // namespace Scine