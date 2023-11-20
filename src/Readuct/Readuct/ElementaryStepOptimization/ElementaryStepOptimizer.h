/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef READUCT_ELEMENTARYSTEPOPTIMIZATION_ELEMENTARYSTEPOPTIMIZER_H
#define READUCT_ELEMENTARYSTEPOPTIMIZATION_ELEMENTARYSTEPOPTIMIZER_H

#include "CostCalculators/ReaductDefault.h"
#include "ElementaryStepOptimizerSettings.h"
#include "EnergiesAndGradientsAlongSpline.h"
#include "ReactionProfile.h"
#include "RecurringProfileCalculator.h"
#include "TypeConverter.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Optimizer/GradientBased/SteepestDescent.h"
#include "Utils/Settings.h"
#include "Utils/UniversalSettings/SettingsNames.h"
#include <Core/Interfaces/Calculator.h>
#include <Eigen/Core>
#include <utility>

namespace Scine {
namespace Readuct {
namespace ElementaryStepOptimization {

/**
 * @brief The base class for all reaction path optimizers. The only purpose of this base class is to hide the template
 *        parameters.
 */
class ElementaryStepOptimizerBase {
 public:
  /// @brief Default constructor.
  ElementaryStepOptimizerBase() = default;

  /// @brief Virtual default destructor.
  virtual ~ElementaryStepOptimizerBase() = default;

  /**
   * @brief Optimize an interpolated elementary step.
   *
   * @return int  The final number of optimization cycles carried out.
   */
  virtual int optimize(Core::Log& log) = 0;

  /**
   * @brief Get the current reaction profile (at the end of an optimization, this is the final reaction profile).
   *
   * @return ReactionProfile The reaction profile.
   */
  virtual ReactionProfile& getReactionProfile() = 0;

  /**
   * @brief Apply the given settings.
   * @param settings The new settings to be applied.
   */
  virtual void setSettings(const Utils::Settings& settings) = 0;

  /**
   * @brief Get the public settings as a Utils::Settings object.
   * @return Settings A settings object with the current settings.
   */
  virtual Utils::Settings getSettings() const = 0;
  /**
   * @brief Get the settings of the calculator used for the energy calculations during the optimization.
   * @return std::shared_ptr<Settings> The settings of the calculator.
   */
  virtual const std::shared_ptr<Utils::Settings> getCalculatorSettings() const = 0;
  /**
   * @brief The underlying convergence check
   *
   * @return GradientBasedCheck the class holding all convergence thresholds.
   */
  virtual Utils::GradientBasedCheck getConvergenceCheck() const = 0;
};

/**
 * @brief A class to optimize reaction paths based on the B-Splines approach described in Alain Vaucher's Ph. D. thesis.
 *
 * @tparam OptimizerType  Expects any of the Optimizer classes.
 */
template<class OptimizerType>
class ElementaryStepOptimizer : public ElementaryStepOptimizerBase {
 public:
  /**
   * @brief Construct a new TestOptimizer object.
   *
   * @param calculator The calculator to be used for the single point/gradient calculations.
   * @param initialElementaryStep The elementary step to be optimized.
   */
  explicit ElementaryStepOptimizer(Core::Calculator& calculator, ReactionProfile initialReactionProfile)
    : _calculator(calculator), _profile(std::move(initialReactionProfile)){};

  /**
   * @brief Apply the given settings.
   * @param settings The new settings to be applied.
   */
  virtual void setSettings(const Utils::Settings& settings) override {
    check.applySettings(settings);
    optimizer.applySettings(settings);
    numberEquidistantPoints = settings.getInt("num_integration_points");
  };

  /**
   * @brief Get the public settings as a Utils::Settings object.
   * @return Settings A settings object with the current settings.
   */
  virtual Utils::Settings getSettings() const override {
    auto settings = ElementaryStepOptimizerSettings<OptimizerType, Utils::GradientBasedCheck>(optimizer, check);
    settings.modifyInt("num_integration_points", numberEquidistantPoints);
    return settings;
  };

  /**
   * @brief Get the settings of the calculator used for the energy calculations during the optimization.
   * @return std::shared_ptr<Settings> The settings of the calculator.
   */
  const std::shared_ptr<Utils::Settings> getCalculatorSettings() const override {
    return std::make_shared<Utils::Settings>(_calculator.settings());
  };

  /**
   * @brief
   *
   * @return int  The final number of optimization cycles carried out.
   */
  int optimize(Core::Log& log) override {
    // Preparations
    EnergiesAndGradientsAlongSpline valuesAlongSpline;
    const auto& elements = _profile.getMolecularSpline().getElements();
    const auto structure = _calculator.getStructure();
    if (!structure || structure->getElements() != elements) {
      Utils::PositionCollection positions = Utils::PositionCollection::Zero(elements.size(), 3); // dummy positions
      Utils::AtomCollection atoms(elements, positions);
      _calculator.setStructure(atoms);
    }

    // Get initial variable values and map to format required by update function
    auto& spline = _profile.getMolecularSpline().getBSpline();
    auto variables = TypeConverter::getInnerControlPointMatrix(spline);
    Eigen::VectorXd variablesVector;
    int nRows = variables.rows();
    int nColumns = variables.cols();
    variablesVector = Eigen::Map<const Eigen::VectorXd>(variables.data(), nRows * nColumns);

    // Create a profile calculator
    RecurringProfileCalculator profileCalculator(_calculator, numberEquidistantPoints);

    // Define update function
    auto const update = [&](const Eigen::VectorXd& parameters, double& value, Eigen::VectorXd& gradients) {
      Eigen::MatrixXd variables;
      variables = Eigen::Map<const Eigen::MatrixXd>(parameters.data(), nRows, nColumns);

      // Insert variables into Spline
      TypeConverter::setInnerControlPoints(spline, variables);

      // Evaluate cost value
      profileCalculator.calculateEnergiesAndGradients(spline);
      auto valuesAlongSpline = profileCalculator.valuesAlongSpline();
      CostBasedOptimization::ReaductDefault costCalculator;
      costCalculator.calculateCost(spline, valuesAlongSpline);
      value = costCalculator.getCost();

      // Evaluate gradients
      Eigen::MatrixXd fullDerivatives = costCalculator.getCostDerivatives();
      int numberInnerControlPoints = spline.controlPointCount() - 2;
      Eigen::MatrixXd derivatives = fullDerivatives.middleRows(1, numberInnerControlPoints);
      gradients = Eigen::Map<const Eigen::VectorXd>(derivatives.data(), nRows * nColumns);

      // Update reaction profile
      _profile.getProfileEnergies() = ProfileEnergies{profileCalculator.getCoordinates(), profileCalculator.getEnergies()};
    };

    // Optimize
    auto cycles = optimizer.optimize(variablesVector, update, check, log);

    return cycles;
  }

  /**
   * @brief Get the current reaction profile (at the end of an optimization, this is the final reaction profile).
   *
   * @return ReactionProfile The reaction profile.
   */
  ReactionProfile& getReactionProfile() override {
    return _profile;
  }

  /**
   * @brief The underlying convergence check
   *
   * @note getter to be accessible via base class
   * @return GradientBasedCheck the class holding all convergence thresholds.
   */
  Utils::GradientBasedCheck getConvergenceCheck() const override {
    return check;
  };

  /// @brief The underlying optimizer (public in order to change its settings).
  OptimizerType optimizer;

  /// @brief The underlying convergence check (public in order to change its settings).
  Utils::GradientBasedCheck check;

  /// @brief The number of interpolation points on the spline.
  int numberEquidistantPoints = numIntegrationPointsDefaultValue;

 private:
  Core::Calculator& _calculator;
  ReactionProfile _profile;
};
} // namespace ElementaryStepOptimization
} // namespace Readuct
} // namespace Scine

#endif // READUCT_ELEMENTARYSTEPOPTIMIZER_H
