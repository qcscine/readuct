/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef READUCT_TSOPTIMIZATIONTASK_H_
#define READUCT_TSOPTIMIZATIONTASK_H_

/* Readuct */
#include "Tasks/Task.h"
/* Scine */
#include "Utils/Geometry/GeometryUtilities.h"
#include "Utils/Geometry/InternalCoordinates.h"
#include <Utils/GeometricDerivatives/HessianUtilities.h>
#include <Utils/GeometricDerivatives/NormalModeAnalysis.h>
#include <Utils/GeometryOptimization/GeometryOptimization.h>
#include <Utils/GeometryOptimization/GeometryOptimizer.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Optimizer/GradientBased/Dimer.h>
#include <Utils/Optimizer/HessianBased/Bofill.h>
#include <Utils/Optimizer/HessianBased/EigenvectorFollowing.h>
/* External */
#include <boost/exception/diagnostic_information.hpp>
#include <boost/filesystem.hpp>
/* std c++ */
#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>

namespace Scine {
namespace Readuct {

class TsOptimizationTask : public Task {
 public:
  /**
   * @brief Construct a new TsOptimizationTask.
   * @param input  The input system names for the task.
   * @param output The output system names for the task.
   * @param logger The logger to/through which all text output will be handled.
   */
  TsOptimizationTask(std::vector<std::string> input, std::vector<std::string> output, std::shared_ptr<Core::Log> logger = nullptr)
    : Task(std::move(input), std::move(output), std::move(logger)) {
  }

  std::string name() const override {
    return "TS Optimization";
  }

  bool run(SystemsMap& systems, Utils::UniversalSettings::ValueCollection taskSettings, bool testRunOnly = false,
           std::vector<std::function<void(const int&, const Utils::AtomCollection&, const Utils::Results&, const std::string&)>>
               observers = {}) const final {
    warningIfMultipleInputsGiven();
    warningIfMultipleOutputsGiven();

    bool silentCalculator = taskSettings.extract("silent_stdout_calculator", true);
    std::shared_ptr<Core::Calculator> calc;
    if (!testRunOnly) { // leave out in case of task chaining --> attention calc is NULL
      // Note: _input is guaranteed not to be empty by Task constructor
      calc = copyCalculator(systems, _input.front(), name());
      Utils::CalculationRoutines::setLog(*calc, true, true, !silentCalculator);

      // Check system size
      if (calc->getStructure()->size() == 1) {
        throw std::runtime_error("Cannot calculate transition state for monoatomic systems.");
      }
    }

    // Generate optimizer
    std::string optimizertype = taskSettings.extract("optimizer", std::string{"bofill"});
    std::transform(optimizertype.begin(), optimizertype.end(), optimizertype.begin(), ::tolower);

    // Check for automatic mode selection
    // sanity check
    if (taskSettings.valueExists("automatic_mode_selection") &&
        (taskSettings.valueExists(optimizertype + "_follow_mode") || taskSettings.valueExists("ev_follow_mode"))) {
      throw std::logic_error("You specified an automatic selection of the mode and gave a specific mode yourself. "
                             "This is not possible. Only give one of those options.");
    }
    std::vector<int> relevantAtoms = taskSettings.extract("automatic_mode_selection", std::vector<int>{});
    // another sanity check
    if (!relevantAtoms.empty()) {
      for (const auto& index : relevantAtoms) {
        if (index < 0) {
          throw std::logic_error("You gave an atom index smaller than 0 in automatic_mode_selection. "
                                 "This does not make sense.");
        }
      }
    }

    using XyzHandler = Utils::XyzStreamHandler;
    auto cout = _logger->output;

    std::shared_ptr<Utils::GeometryOptimizerBase> optimizer;
    if (optimizertype == "bofill") {
      auto tmp = std::make_shared<Utils::GeometryOptimizer<Utils::Bofill>>(*calc);
      optimizer = std::move(tmp);
    }
    else if (optimizertype == "ef" || optimizertype == "ev" || optimizertype == "evf" ||
             optimizertype == "eigenvectorfollowing" || optimizertype == "eigenvector_following") {
      auto tmp = std::make_shared<Utils::GeometryOptimizer<Utils::EigenvectorFollowing>>(*calc);
      optimizer = std::move(tmp);
    }
    else if (optimizertype == "dimer") {
      auto tmp = std::make_shared<Utils::GeometryOptimizer<Utils::Dimer>>(*calc);

      /* get if coordinates are transformed to know whether a provided guess vector has to be transformed */
      auto coordinateSystem = tmp->coordinateSystem;
      if (taskSettings.valueExists("geoopt_coordinate_system")) {
        coordinateSystem = Utils::CoordinateSystemInterpreter::getCoordinateSystemFromString(
            taskSettings.getString("geoopt_coordinate_system"));
      }

      /* Check for possible guess vectors */
      bool useEigenvectorForFirstStep = taskSettings.extract("dimer_calculate_hessian_once", false);
      if (!useEigenvectorForFirstStep && !relevantAtoms.empty()) {
        _logger->warning << "Specified Dimer optimizer with automatic mode selection, this means that a Hessian has to "
                            "be calculated.\n";
      }
      /* set to true if mode is specified or automated selection of mode is set in Settings */
      if (taskSettings.valueExists("dimer_follow_mode") || !relevantAtoms.empty()) {
        useEigenvectorForFirstStep = true;
      }
      int followedMode = taskSettings.extract("dimer_follow_mode", 0);
      std::string guessVectorFileName = taskSettings.extract("dimer_guess_vector_file", std::string(""));
      std::vector<std::string> discreteGuessesFileNames =
          taskSettings.extract("dimer_discrete_guesses", std::vector<std::string>{"", ""});

      /* sanity checks for Dimer init options */
      if (discreteGuessesFileNames.size() != 2) {
        throw std::runtime_error("Problem with discrete guesses. You may only give list with two filepaths.");
      }
      // calc cannot be checked in test run
      if (!testRunOnly && !useEigenvectorForFirstStep && calc->results().has<Utils::Property::Hessian>()) {
        _logger->warning << "Did not specify to use eigenvector for initialization of Dimer, although your system "
                            "already includes a Hessian. The Hessian will therefore be ignored. "
                            "If you want to change that, activate the setting 'dimer_calculate_hessian_once': true"
                         << Core::Log::endl;
      }
      if (useEigenvectorForFirstStep) {
        if (!guessVectorFileName.empty()) {
          throw std::logic_error("You specified to use the eigenvector of the Hessian as Dimer initialization, "
                                 "but also specified a file to read a guess vector for the initialization. "
                                 "Only select one.");
        }
        if (!discreteGuessesFileNames[0].empty()) {
          throw std::logic_error("You specified to use the eigenvector of the Hessian as Dimer initialization, "
                                 "but also specified files to read discrete guesses for the initialization. "
                                 "Only select one.");
        }
      }
      else if (!guessVectorFileName.empty() && !discreteGuessesFileNames[0].empty()) {
        throw std::logic_error("You specified a file to read a guess vector for the Dimer initialization. "
                               "but also specified files to read discrete guesses for the initialization. "
                               "Only select one.");
      }

      if (testRunOnly) {
        ; // skip all dimer if statements here, because they require a Hessian calculation or file look-up
      }
      /* Generate and pass initial eigenvector and tell optimizer to not perform the first rotation */
      else if (useEigenvectorForFirstStep) {
        tmp->optimizer.skipFirstRotation = true;
        auto modes = getNormalModes(*calc);
        auto system = calc->getStructure();
        if (!relevantAtoms.empty()) { // clash with dimer_follow_mode option is already checked prior
          followedMode = selectMode(*calc, relevantAtoms);
        }
        auto mode = modes.getMode(followedMode); // range check for followedMode is included here
        /* Currently Hessian not possible in internal coordinates */
        if (coordinateSystem == Utils::CoordinateSystem::Internal) {
          /* Transform vector of normal mode into internal coordinates */
          Utils::InternalCoordinates internalCoordinates(*system);
          tmp->optimizer.guessVector = std::make_shared<Eigen::VectorXd>(internalCoordinates.coordinatesToInternal(mode));
        }
        else if (coordinateSystem == Utils::CoordinateSystem::CartesianWithoutRotTrans) {
          /* Can give proper invH for all but internal */
          /* Calculate inverse Hessian with one mode flipped for BFGS within Dimer optimizer */
          Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
          auto hessian = calc->results().get<Utils::Property::Hessian>(); // guaranteed to have it from getNormalModes
          es.compute(hessian);
          int nEigenValues = es.eigenvalues().size();
          Eigen::MatrixXd eigenValueMatrix = Eigen::MatrixXd::Identity(nEigenValues, nEigenValues);
          eigenValueMatrix.diagonal() = es.eigenvalues();
          eigenValueMatrix(followedMode, followedMode) = -es.eigenvalues()[followedMode];
          Eigen::MatrixXd switchedHessian = es.eigenvectors() * eigenValueMatrix * es.eigenvectors().transpose();
          /* Transform vector of normal mode into internal coordinates */
          Utils::InternalCoordinates internalCoordinates(*system, true);
          tmp->optimizer.guessVector = std::make_shared<Eigen::VectorXd>(internalCoordinates.coordinatesToInternal(mode));
          /* Get rid of translation and rotation in Hessian */
          Eigen::MatrixXd switchedTransformedHessian = internalCoordinates.hessianToInternal(switchedHessian);
          tmp->optimizer.invH = switchedTransformedHessian.inverse();
        }
        else if (coordinateSystem == Utils::CoordinateSystem::Cartesian) {
          /* no not give invH, because this leads to a lot of translations */
          tmp->optimizer.guessVector =
              std::make_shared<Eigen::VectorXd>(Eigen::Map<const Eigen::VectorXd>(mode.data(), mode.cols() * mode.rows()));
        }
        else {
          throw std::logic_error("Unknown coordinate system " +
                                 Utils::CoordinateSystemInterpreter::getStringFromCoordinateSystem(coordinateSystem));
        }
      }
      /* Read vector from file to give as guess vector */
      else if (!guessVectorFileName.empty()) {
        /* push data into std::vector */
        std::vector<double> container;
        double num = 0.0;
        std::ifstream file(guessVectorFileName, std::ios::in);
        if (!file.is_open()) {
          throw std::runtime_error("Problem when opening file " + guessVectorFileName);
        }
        while (file >> num) {
          container.push_back(num);
        }
        file.close();
        Eigen::Map<Eigen::VectorXd> guessVector(container.data(), container.size());
        auto system = calc->getStructure();
        /* Provided vector is Cartesian and internal shall be used -> transform to internal */
        if (coordinateSystem != Utils::CoordinateSystem::Cartesian &&
            static_cast<unsigned long>(guessVector.size()) == 3 * system->getElements().size()) {
          std::shared_ptr<Utils::InternalCoordinates> transformation = nullptr;
          if (coordinateSystem == Utils::CoordinateSystem::Internal) {
            transformation = std::make_shared<Utils::InternalCoordinates>(*system);
          }
          else if (coordinateSystem == Utils::CoordinateSystem::CartesianWithoutRotTrans) {
            transformation = std::make_shared<Utils::InternalCoordinates>(*system, true);
          }
          // vector to PositionCollection
          Utils::PositionCollection guessPosition =
              Eigen::Map<Utils::PositionCollection>(guessVector.data(), guessVector.size() / 3, 3);
          try {
            tmp->optimizer.guessVector =
                std::make_shared<Eigen::VectorXd>(transformation->coordinatesToInternal(guessPosition));
          }
          catch (const Utils::InternalCoordinatesException& e) {
            _logger->warning << "Could not transform given Cartesian coordinates to Internal coordinates.\n"
                             << "Omitting guess and initializing dimer with random vector. " << Core::Log::endl;
          }
        }
        /* Provided vector is not Cartesian and Cartesian shall be used -> transform to Cartesian */
        else if (coordinateSystem == Utils::CoordinateSystem::Cartesian &&
                 static_cast<unsigned long>(guessVector.size()) != 3 * system->getElements().size()) {
          try {
            std::shared_ptr<Utils::InternalCoordinates> transformation = nullptr;
            if (static_cast<unsigned long>(guessVector.size()) == 3 * system->getElements().size() - 6 ||
                static_cast<unsigned long>(guessVector.size()) == 3 * system->getElements().size() - 5) {
              // cartesianWithoutRotTrans because 5 or 6 missing
              transformation = std::make_shared<Utils::InternalCoordinates>(*system, true);
            }
            else {
              transformation = std::make_shared<Utils::InternalCoordinates>(*system);
            }
            /* single guess vector cannot get backtransformed directly to Cartesian coords */
            /* first get current position and transform to internal */
            auto cartPos0 = system->getPositions();
            auto ircVec0 = transformation->coordinatesToInternal(cartPos0);
            /* add guessvector */
            Eigen::VectorXd ircVec1 = ircVec0 + guessVector;
            /* transform this new vector back to Cartesian and diff of 2 Cartesian PositionCollection is guessvector */
            auto cartPos1 = transformation->coordinatesToCartesian(ircVec1);
            Utils::PositionCollection posDiff = cartPos1 - cartPos0;
            tmp->optimizer.guessVector =
                std::make_shared<Eigen::VectorXd>(Eigen::Map<Eigen::VectorXd>(posDiff.data(), posDiff.size()));
          }
          catch (const Utils::InternalCoordinatesException& e) {
            _logger->warning << "Could not transform given non-Cartesian coordinates to Cartesian.\n"
                             << "Omitting guess and initializing dimer with random vector. " << Core::Log::endl;
          }
        }
        /* Provided vector should fit with wanted coordinates, optimizer later checks the correct length again */
        else {
          tmp->optimizer.guessVector = std::make_shared<Eigen::VectorXd>(guessVector);
        }
      }
      /* Read 2 structures from files to form guess vector from discrete difference of the two */
      else if (!discreteGuessesFileNames.at(0).empty() && !discreteGuessesFileNames.at(1).empty()) {
        Utils::AtomCollection system;
        /* read first */
        std::string filepath1 = discreteGuessesFileNames.at(0);
        std::ifstream file1(filepath1, std::ios::in);
        if (!file1.is_open()) {
          throw std::runtime_error("Problem when opening file " + filepath1);
        }
        system = XyzHandler::read(file1);
        auto position1 = system.getPositions();
        file1.close();

        /* read second */
        std::string filepath2 = discreteGuessesFileNames.at(1);
        std::ifstream file2(filepath2, std::ios::in);
        if (!file2.is_open()) {
          throw std::runtime_error("Problem when opening file " + filepath2);
        }
        system = XyzHandler::read(file2);
        auto position2 = system.getPositions();
        file2.close();

        Utils::PositionCollection positionDiff = position2 - position1;
        Eigen::Map<Eigen::VectorXd> mode(positionDiff.data(), positionDiff.size());
        std::shared_ptr<Utils::InternalCoordinates> transformation = nullptr;
        if (coordinateSystem == Utils::CoordinateSystem::Internal) {
          transformation = std::make_shared<Utils::InternalCoordinates>(system);
        }
        else if (coordinateSystem == Utils::CoordinateSystem::CartesianWithoutRotTrans) {
          transformation = std::make_shared<Utils::InternalCoordinates>(system, true);
        }
        if (transformation) {
          tmp->optimizer.guessVector = std::make_shared<Eigen::VectorXd>(transformation->coordinatesToInternal(positionDiff));
        }
        else {
          tmp->optimizer.guessVector = std::make_shared<Eigen::VectorXd>(mode);
        }
      }
      optimizer = std::move(tmp);
    } // end of Dimer handling
    else {
      throw std::runtime_error(
          "Unknown Optimizer requested for TS optimization, available are: Bofill, EVF and Dimer!");
    }

    // Read and delete special settings
    bool stopOnError = stopOnErrorExtraction(taskSettings);
    // Apply settings
    auto settings = optimizer->getSettings();
    settings.merge(taskSettings);
    // set automatically selected mode, if requested
    // for Dimer this happened already, because it is GradientBased and does not know anything about modes / Hessians
    // do not do in testRun because calc would be empty
    if (!relevantAtoms.empty() && !testRunOnly) {
      // necessary because multiple strings are accepted for eigenvectorfollowing optimizer
      if (optimizertype != "dimer" && optimizertype != "bofill") {
        settings.modifyValue("ev_follow_mode", selectMode(*calc, relevantAtoms));
      }
      else if (optimizertype == "bofill") {
        settings.modifyValue("bofill_follow_mode", selectMode(*calc, relevantAtoms));
      }
    }
    if (!settings.valid()) {
      settings.throwIncorrectSettings();
    }
    optimizer->setSettings(settings);
    if (!testRunOnly && !Utils::GeometryOptimization::settingsMakeSense(*optimizer)) {
      throw std::logic_error("The given calculator settings are too inaccurate for the given convergence criteria of "
                             "this optimization Task");
    }

    // If no errors encountered until here, the basic settings should be alright
    if (testRunOnly) {
      return true;
    }

    // Add observer
    // Trajectory stream
    Utils::XyzStreamHandler writer;
    const std::string& partialOutput = (!_output.empty() ? _output[0] : _input[0]);
    boost::filesystem::path dir(partialOutput);
    boost::filesystem::create_directory(dir);
    boost::filesystem::path trjfile(partialOutput + ".tsopt.trj.xyz");
    std::ofstream trajectory((dir / trjfile).string(), std::ofstream::out);
    cout.printf("%7s %16s %16s %16s %16s\n", "Cycle", "Energy", "Energy Diff.", "Step RMS", "Max. Step");
    double oldEnergy = 0.0;
    Eigen::VectorXd oldParams;
    auto func = [&](const int& cycle, const double& energy, const Eigen::VectorXd& params) {
      if (oldParams.size() != params.size()) {
        oldParams = params;
      }
      auto diff = (params - oldParams).eval();
      cout.printf("%7d %+16.9f %+16.9f %+16.9f %+16.9f\n", cycle, energy, energy - oldEnergy,
                  sqrt(diff.squaredNorm() / diff.size()), diff.cwiseAbs().maxCoeff());
      oldEnergy = energy;
      oldParams = params;
      auto structure = calc->getStructure();
      XyzHandler::write(trajectory, *structure, std::to_string(energy));
    };
    optimizer->addObserver(func);

    // Add custom observers
    auto customObservers = [&calc, &observers](const int& cycle, const double& /*energy*/, const Eigen::VectorXd& /*params*/) {
      for (auto& observer : observers) {
        auto atoms = calc->getStructure();
        Utils::Results& results = calc->results();
        observer(cycle, *atoms, results, "ts_optimization");
      }
    };
    optimizer->addObserver(customObservers);

    // Run optimization
    auto structure = calc->getStructure();
    int cycles = 0;
    try {
      cycles = optimizer->optimize(*structure, *_logger);
    }
    catch (...) {
      XyzHandler::write(trajectory, *calc->getStructure());
      trajectory.close();
      _logger->error << "TS Optimization failed with error!" << Core::Log::endl;
      if (stopOnError) {
        throw;
      }
      _logger->error << boost::current_exception_diagnostic_information() << Core::Log::endl;
      return false;
    }
    trajectory.close();

    int maxiter = settings.getInt("convergence_max_iterations");
    if (cycles < maxiter) {
      cout << Core::Log::endl
           << "    Converged after " << cycles << " iterations." << Core::Log::endl
           << Core::Log::endl;
    }
    else {
      cout << Core::Log::endl
           << "    Stopped after " << maxiter << " iterations." << Core::Log::endl
           << Core::Log::endl;
      if (stopOnError) {
        throw std::runtime_error("Problem: TS optimization did not converge.");
      }
    }

    // Print/Store results
    systems[partialOutput] = calc;
    boost::filesystem::path xyzfile(partialOutput + ".xyz");
    std::ofstream xyz((dir / xyzfile).string(), std::ofstream::out);
    XyzHandler::write(xyz, *(calc->getStructure()));
    xyz.close();

    return cycles < maxiter;
  }

 private:
  int selectMode(Core::Calculator& calc, const std::vector<int>& relevantAtoms) const {
    _logger->output << "  Automatically selecting normal mode " << Core::Log::endl;
    int nAtoms = calc.getStructure()->size();
    for (const auto& index : relevantAtoms) {
      if (index >= nAtoms) {
        throw std::logic_error("You gave an atom index larger than the number of atoms in automatic_mode_selection. "
                               "This does not make sense.");
      }
    }
    Utils::NormalModesContainer modes = getNormalModes(calc);
    auto waveNumbers = modes.getWaveNumbers();
    std::vector<double> contributions;
    // cycle all imag. frequencies
    for (int i = 0; i < static_cast<int>(waveNumbers.size()); ++i) {
      if (waveNumbers[i] >= 0.0) {
        if (i == 0) {
          throw std::runtime_error("Structure has no imaginary normal mode.");
        }
        break;
      }
      auto mode = modes.getMode(i);
      double contribution = 0.0;
      // cycle all atoms
      for (int j = 0; j < mode.size(); ++j) {
        if (std::find(relevantAtoms.begin(), relevantAtoms.end(), j) != relevantAtoms.end()) {
          contribution += mode.row(j).norm(); // add norm of vector of relevant atom contributing to mode
        }
      }
      contributions.push_back(contribution);
    }
    int selection = std::distance(contributions.begin(), std::max_element(contributions.begin(), contributions.end()));
    double tolerance = contributions[selection] * 0.1;
    // Contributions are sorted by wave number
    //   -> pick the first one that fits into the tolerance.
    for (unsigned int i = 0; i < contributions.size(); i++) {
      if ((contributions[selection] - contributions[i]) < tolerance) {
        selection = i;
        break;
      }
    }
    _logger->output << "    Automatically selected mode " << std::to_string(selection) << " with wave number "
                    << std::to_string(waveNumbers[selection]) << " cm^-1" << Core::Log::endl;
    return selection;
  }

  static Utils::NormalModesContainer getNormalModes(Core::Calculator& calc) {
    if (!calc.results().has<Utils::Property::Hessian>()) {
      calc.setRequiredProperties(Utils::Property::Hessian);
      calc.calculate("Hessian Calculation");
    }
    auto hessian = calc.results().get<Utils::Property::Hessian>();
    auto system = calc.getStructure();
    return Utils::NormalModeAnalysis::calculateNormalModes(hessian, system->getElements(), system->getPositions());
  }
};

} // namespace Readuct
} // namespace Scine

#endif // READUCT_TSOPTIMIZATIONTASK_H_
