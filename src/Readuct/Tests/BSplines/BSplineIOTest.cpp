/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Readuct/BSplines/IO.h>
#include <Readuct/BSplines/InterpolationGenerator.h>
#include <Readuct/BSplines/MolecularSpline.h>
#include <gmock/gmock.h>
#include <sstream>

using namespace testing;
namespace Scine {
namespace Readuct {
using namespace BSplines;
namespace Tests {

class ABSplineIO : public Test {
 public:
  BSpline randomBSpline;
  Utils::ElementTypeCollection ec;

 protected:
  void SetUp() override {
    const int numberDimensions = 4;
    const int numberPoints = 20;
    Eigen::MatrixXd randomPoints = Eigen::MatrixXd::Random(numberPoints, numberDimensions);
    InterpolationGenerator generator(randomPoints);
    randomBSpline = generator.generateBSpline();

    ec.push_back(Utils::ElementType::Ra);
    ec.push_back(Utils::ElementType::Ne);
    ec.push_back(Utils::ElementType::Br);
    ec.push_back(Utils::ElementType::H);
  }
};

TEST_F(ABSplineIO, WritingAndReadingDeliversSameMolecularSpline) {
  MolecularSpline molecularSpline{ec, randomBSpline};

  std::stringstream stream;
  IO::writeToStream(stream, molecularSpline);

  auto newMolecularSpline = IO::readFromStream(stream);
  const auto& newElements = newMolecularSpline.getElements();
  const auto& newSpline = newMolecularSpline.getBSpline();

  ASSERT_TRUE(ec == newElements);

  ASSERT_THAT(newSpline.getDegree(), Eq(randomBSpline.getDegree()));
  ASSERT_THAT(newSpline.getKnotVector(), Eq(randomBSpline.getKnotVector()));
  ASSERT_THAT(newSpline.getControlPointMatrix(), Eq(randomBSpline.getControlPointMatrix()));
}

} // namespace Tests
} // namespace Readuct
} // namespace Scine