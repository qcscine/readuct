/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Readuct/BSplines/BSpline.h>
#include <Readuct/BSplines/GeometricManipulations.h>
#include <Readuct/BSplines/InterpolationGenerator.h>
#include <Utils/Geometry.h>
#include <Utils/QuaternionFit.h>
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Readuct {
using namespace BSplines;
namespace Tests {

class ABSplineGeometryManipulation : public Test {
 public:
  BSpline randomBSpline;
  const int numberDimensions = 19;
  const int numberPoints = 20;

 protected:
  void SetUp() override {
    Eigen::MatrixXd randomPoints = Eigen::MatrixXd::Random(numberPoints, numberDimensions);
    InterpolationGenerator generator(randomPoints);
    randomBSpline = generator.generateBSpline();
  }
};

TEST_F(ABSplineGeometryManipulation, CanTranslateAWholeBSpline) {
  Eigen::VectorXd translation = Eigen::VectorXd::Random(numberDimensions);

  auto translatedBSpline = GeometricManipulations::translate(randomBSpline, translation);

  for (auto u : {0.001, 0.009, 0.44, 0.77, 0.9110}) {
    Eigen::VectorXd expected = randomBSpline.evaluate(u) + translation;
    Eigen::VectorXd result = translatedBSpline.evaluate(u);
    ASSERT_TRUE(result.isApprox(expected));
  }
}

TEST_F(ABSplineGeometryManipulation, CanTranslateASingleDimension) {
  double randomNumber = -0.6757;
  Eigen::VectorXd translation(1);
  translation << randomNumber;

  auto translatedBSpline = GeometricManipulations::translate(randomBSpline, 2, 3, translation);
  for (auto u : {0.001, 0.009, 0.44, 0.77, 0.9110}) {
    Eigen::VectorXd expected = randomBSpline.evaluate(u);
    expected(2) += randomNumber;

    Eigen::VectorXd result = translatedBSpline.evaluate(u);
    ASSERT_TRUE(result.isApprox(expected));
  }
}

TEST_F(ABSplineGeometryManipulation, CanTranslateMultipleTimesNDimensions) {
  Eigen::Vector2d translation = Eigen::Vector2d::Random();

  // Translate two times two dimensions
  auto translatedBSpline = GeometricManipulations::translate(randomBSpline, 5, 9, translation);

  for (auto u : {0.001, 0.009, 0.44, 0.77, 0.9110}) {
    Eigen::VectorXd expected = randomBSpline.evaluate(u);
    expected.segment(5, 2) += translation;
    expected.segment(7, 2) += translation;
    Eigen::VectorXd result = translatedBSpline.evaluate(u);
    ASSERT_TRUE(result.isApprox(expected));
  }
}

TEST_F(ABSplineGeometryManipulation, CanRecenter3NDSpline) {
  auto recenteredSpline = GeometricManipulations::recenter3d(randomBSpline, 1, 19);

  // All points are linear combinations and must also be centered at the origin
  for (auto u : {0.001, 0.009, 0.44, 0.77, 0.9110}) {
    Eigen::VectorXd oldV = randomBSpline.evaluate(u);
    Eigen::VectorXd newV = recenteredSpline.evaluate(u);
    ASSERT_THAT(newV(0), DoubleEq(oldV(0)));

    Eigen::Vector3d coordinateSum = Eigen::Vector3d::Zero();
    for (int i = 1; i < 19; i += 3) {
      Eigen::Vector3d vec = newV.segment(i, 3);
      coordinateSum += vec;
    }

    ASSERT_TRUE(coordinateSum.isZero());
  }
}

TEST_F(ABSplineGeometryManipulation, AllowsThreeDimensionalTransformations) {
  double randomScaling = 2.333;
  double randomAngle = 1.56789;
  Eigen::Vector3d randomTranslation = Eigen::Vector3d::Random();
  Eigen::Vector3d randomAxis = Eigen::Vector3d::Random().normalized();

  using MyTransform = Eigen::Transform<double, 3, Eigen::Affine>;
  MyTransform t = MyTransform::Identity();
  t *= Eigen::AngleAxisd(randomAngle, randomAxis);
  t *= Eigen::Scaling(randomScaling);
  t *= Eigen::Translation<double, 3>(randomTranslation);

  auto transformedSpline = GeometricManipulations::transform3d(randomBSpline, 1, 19, t);
  for (auto u : {0.001, 0.009, 0.44, 0.77, 0.9110}) {
    Eigen::VectorXd reference = randomBSpline.evaluate(u);
    Eigen::VectorXd result = transformedSpline.evaluate(u);

    ASSERT_THAT(result(0), DoubleEq(reference(0)));
    for (int i = 1; i < 19; i += 3) {
      Eigen::Vector3d vec = reference.segment(i, 3);
      Eigen::Vector3d expected = t * vec;
      Eigen::Vector3d obtained = result.segment(i, 3);
      ASSERT_TRUE(obtained.isApprox(expected));
    }
  }
}

TEST_F(ABSplineGeometryManipulation, AllowsThreeDimensionalRotation) {
  double randomAngle = 1.56789;
  Eigen::Vector3d randomAxis = Eigen::Vector3d::Random().normalized();

  Eigen::AngleAxisd rotation = Eigen::AngleAxisd(randomAngle, randomAxis);

  auto transformedSpline = GeometricManipulations::rotate3d(randomBSpline, 1, 19, rotation);
  for (auto u : {0.001, 0.009, 0.44, 0.77, 0.9110}) {
    Eigen::VectorXd reference = randomBSpline.evaluate(u);
    Eigen::VectorXd result = transformedSpline.evaluate(u);

    ASSERT_THAT(result(0), DoubleEq(reference(0)));
    for (int i = 1; i < 19; i += 3) {
      Eigen::Vector3d vec = reference.segment(i, 3);
      Eigen::Vector3d expected = rotation.toRotationMatrix() * vec;
      Eigen::Vector3d obtained = result.segment(i, 3);
      ASSERT_TRUE(obtained.isApprox(expected));
    }
  }
}

TEST_F(ABSplineGeometryManipulation, CanRemoveRotoTranslation) {
  BSpline spline;
  const int numberDimensions = 18;
  const int numberPoints = 20;
  Eigen::MatrixXd randomPoints = Eigen::MatrixXd::Random(numberPoints, numberDimensions);
  InterpolationGenerator generator(randomPoints);
  spline = generator.generateBSpline();

  auto splineWithoutRotoTranslation = GeometricManipulations::removeRotoTranslation(spline);
  const auto& controlPoints = splineWithoutRotoTranslation.getControlPointMatrix();
  for (int i = 1; i < spline.controlPointCount(); ++i) {
    const Eigen::VectorXd& refCP = controlPoints.row(i - 1);
    const Eigen::VectorXd& fitCP = controlPoints.row(i);

    Eigen::MatrixXd refMat = Utils::Geometry::positionVectorToMatrix(refCP);
    Eigen::MatrixXd fitMat = Utils::Geometry::positionVectorToMatrix(fitCP);

    Utils::QuaternionFit fit(refMat, fitMat);

    Eigen::Vector3d translationVector = fit.getTransVector();
    Eigen::Matrix3d rotation = fit.getRotationMatrix();

    Eigen::AngleAxisd angleRotation{rotation};

    ASSERT_THAT(angleRotation.angle(), DoubleNear(0.0, 1e-6));
    ASSERT_TRUE(translationVector.isZero(1e-6));
  }
}

} // namespace Tests
} // namespace Readuct
} // namespace Scine