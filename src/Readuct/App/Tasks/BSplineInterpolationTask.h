/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef READUCT_BSPLINEINTERPOLATIONTASK_H_
#define READUCT_BSPLINEINTERPOLATIONTASK_H_

#include "../../Readuct/ElementaryStepOptimization/ElementaryStepOptimizer.h"
#include "../../Readuct/ElementaryStepOptimization/ReactionProfile.h"
#include "Tasks/Task.h"
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/MolecularTrajectoryIO.h>
#include <Utils/Math/BSplines/FixedEndsPenalizedLeastSquaresGenerator.h>
#include <Utils/Math/BSplines/InterpolationGenerator.h>
#include <Utils/Math/BSplines/LinearInterpolator.h>
#include <Utils/Math/BSplines/MolecularSpline.h>
#include <Utils/Math/QuaternionFit.h>
#include <Utils/MolecularTrajectory.h>
#include <Utils/Optimizer/GradientBased/Lbfgs.h>
#include <Utils/Optimizer/GradientBased/SteepestDescent.h>
#include <boost/filesystem.hpp>
#include <fstream>

namespace Scine {
namespace Readuct {

class BSplineInterpolationTask : public Task {
 public:
  /**
   * @brief Construct a new BSplineInterpolationTask.
   * @param input  The input system names for the task.
   * @param output The output system names for the task.
   * @param logger The logger to/through which all text output will be handled.
   */
  BSplineInterpolationTask(std::vector<std::string> input, std::vector<std::string> output,
                           std::shared_ptr<Core::Log> logger = nullptr)
    : Task(std::move(input), std::move(output), std::move(logger)) {
  }

  std::string name() const override {
    return "BSpline Interpolation";
  }

  bool run(std::map<std::string, std::shared_ptr<Core::Calculator>>& systems,
           Utils::UniversalSettings::ValueCollection taskSettings, bool testRunOnly = false,
           std::vector<std::function<void(const int&, const Utils::AtomCollection&, const Utils::Results&, const std::string&)>>
               observers = {}) const final {
    if (_input.size() != 2) {
      throw std::logic_error("The B-Spline task needs two input systems.");
    }
    warningIfMultipleOutputsGiven();
    if (observers.size() > 0) {
      throw std::logic_error(
          "BSplineInterpolationTask does not feature algorithm accepting observers, yet one was given");
    }

    // Read and set user-specified settings
    alignStructuresBeforeInterpolation_ = taskSettings.extract("align_structures", alignStructuresBeforeInterpolation_);
    numberStructuresForMolecularTrajectory_ = taskSettings.extract("num_structures", numberStructuresForMolecularTrajectory_);
    optimize_ = taskSettings.extract("optimize", optimize_);
    numberControlPointsForInterpolation_ = taskSettings.extract("num_control_points", numberControlPointsForInterpolation_);
    optimizer_ = taskSettings.extract("optimizer", optimizer_);
    std::transform(optimizer_.begin(), optimizer_.end(), optimizer_.begin(), ::tolower);
    extractTsGuess_ = taskSettings.extract("extract_ts_guess", extractTsGuess_);
    extractTsGuessNeighbours_ = taskSettings.extract("extract_ts_guess_neighbours", extractTsGuessNeighbours_);
    tangentFileName_ = taskSettings.extract("tangent_file", tangentFileName_);
    coordinateThresholdForMaximumExtraction_ =
        taskSettings.extract("extract_threshold", coordinateThresholdForMaximumExtraction_);
    bool silentCalculator = taskSettings.extract("silent_stdout_calculator", true);

    // If no errors encountered until here, the basic settings should be alright
    if (testRunOnly) {
      return true; // leave out rest in case of task chaining
    }

    // Note: _input is guaranteed not to be empty by Task constructor
    auto calc = copyCalculator(systems, _input.front(), name());
    auto secondCalculator = systems.at(_input.back());
    Utils::CalculationRoutines::setLog(*calc, true, true, !silentCalculator);
    Utils::CalculationRoutines::setLog(*secondCalculator, true, true, !silentCalculator);
    if (calc->settings() != secondCalculator->settings()) {
      _logger->warning
          << "  Warning: The given systems have different settings. Only taking first and ignoring second.\n";
    }

    // Get inputs
    auto start = calc->getStructure();
    auto end = secondCalculator->getStructure();
    if (start->getElements() != end->getElements()) {
      throw std::logic_error("The provided structures have different elements. Impossible to find path between them.");
    }

    // Check system size
    if (start->size() == 1) {
      throw std::runtime_error("Cannot perform " + name() + " for monoatomic systems.");
    }

    Utils::MolecularTrajectory trajectoryGuess = Utils::MolecularTrajectory();
    trajectoryGuess.setElementTypes(start->getElements());
    if (taskSettings.valueExists("trajectory_guess")) {
      std::vector<std::string> paths = taskSettings.getStringList("trajectory_guess");
      for (unsigned long i = 0; i < paths.size(); ++i) {
        Utils::MolecularTrajectory traj = readTrajectory(paths[i]);
        if (traj.getElementTypes() != trajectoryGuess.getElementTypes()) {
          throw std::logic_error("The provided trajectory guess no. " + std::to_string(i) +
                                 " has different elements than the provided input structures.");
        }
        for (auto& step : traj) {
          trajectoryGuess.push_back(step);
        }
      }
      taskSettings.dropValue("trajectory_guess");
    }

    // Interpolate between start and end structure
    auto cout = _logger->output;
    cout << "  Interpolating Reaction Path\n";
    ElementaryStepOptimization::ReactionProfile interpolatedProfile = interpolateElementaryStep(*start, *end, trajectoryGuess);
    Utils::MolecularTrajectory interpolatedTrajectory = convertProfileToTrajectory(interpolatedProfile);
    // Print/Store results
    using Writer = Utils::ChemicalFileHandler;
    const std::string& outputSystem = ((!_output.empty()) ? _output[0] : _input[0]);
    boost::filesystem::path dir(outputSystem);
    boost::filesystem::create_directory(dir);
    systems[outputSystem] = calc;
    boost::filesystem::path xyzfile(outputSystem + "_interpolated.xyz");
    writeTrajectory(interpolatedTrajectory, (dir / xyzfile).string());
    cout << "  Interpolating Complete\n\n";

    // Optimize the interpolated path
    Utils::MolecularTrajectory optimizedTrajectory;
    ElementaryStepOptimization::ReactionProfile optimizedProfile;
    bool converged = true;
    if (optimize_) {
      cout << "  Optimizing Reaction Path\n";
      auto result = optimizeElementaryStep(std::move(interpolatedProfile), calc, taskSettings);
      optimizedProfile = result.first;
      converged = result.second;
      optimizedTrajectory = convertProfileToTrajectory(optimizedProfile);
      // Output
      boost::filesystem::path xyzfile2(outputSystem + "_optimized.xyz");
      writeTrajectory(optimizedTrajectory, (dir / xyzfile2).string());
      std::string info = (converged) ? "Complete" : "Stopped";
      cout << "  Optimization " + info + "\n\n";
    }

    std::vector<Utils::AtomCollection> tsGuess;
    if (optimize_ && extractTsGuess_) {
      cout << "  Extracting Transition State Guess\n";
      tsGuess = extractTransitionStateGuessStructure(optimizedProfile, calc, dir);
      // Output
      boost::filesystem::path xyzfile3 = dir / (outputSystem + "_tsguess.xyz");
      Writer::write(xyzfile3.string(), tsGuess.at(0));

      if (extractTsGuessNeighbours_) {
        boost::filesystem::path xyzfile3_1 = dir / (outputSystem + "_tsguess-1.xyz");
        Writer::write(xyzfile3_1.string(), tsGuess.at(1));

        boost::filesystem::path xyzfile3_2 = dir / (outputSystem + "_tsguess+1.xyz");
        Writer::write(xyzfile3_2.string(), tsGuess.at(2));
      }

      cout << "  Extraction Complete\n\n";
    }
    cout << Core::Log::nl << Core::Log::endl;
    return converged;
  }

 private:
  ElementaryStepOptimization::ReactionProfile
  interpolateElementaryStep(const Utils::AtomCollection& start, Utils::AtomCollection end,
                            Utils::MolecularTrajectory trajectoryGuess = Utils::MolecularTrajectory()) const {
    auto cout = _logger->output;
    if (alignStructuresBeforeInterpolation_) {
      auto positionsToAlign = end.getPositions();
      Utils::Geometry::Manipulations::alignPositions(start.getPositions(), positionsToAlign);
      end.setPositions(std::move(positionsToAlign));
    }

    Eigen::VectorXd startVector = Eigen::Map<const Eigen::VectorXd>(
        start.getPositions().data(), start.getPositions().cols() * start.getPositions().rows());
    Eigen::VectorXd endVector =
        Eigen::Map<const Eigen::VectorXd>(end.getPositions().data(), end.getPositions().cols() * end.getPositions().rows());
    Utils::BSplines::BSpline spline;
    if (trajectoryGuess.empty()) {
      spline = Utils::BSplines::LinearInterpolator::generate(startVector, endVector, numberControlPointsForInterpolation_);
    }
    else {
      if (start.size() != trajectoryGuess.molecularSize()) {
        throw std::logic_error("The given trajectory guess has a different size than the given start system.");
      }
      /* Remove frames from trajectory at end and beginning which are too similar to given endpoints */
      bool viableTrajectory = true;
      double rmsdThreshold = 1.0;
      bool stop = false;
      while (!stop) {
        Utils::QuaternionFit fit(start.getPositions(), trajectoryGuess.front());
        if (fit.getRMSD() <= rmsdThreshold) {
          trajectoryGuess.erase(trajectoryGuess.begin());
          if (trajectoryGuess.empty()) {
            _logger->warning
                << "Warning: no viable trajectory guess provided, all frames are too similar to the endpoints, now "
                   "performing linear interpolation."
                << Core::Log::endl;
            viableTrajectory = false;
            stop = true;
          }
        }
        else {
          stop = true;
        }
      }
      if (viableTrajectory) {
        stop = false;
      }
      while (!stop) {
        Utils::QuaternionFit fit(end.getPositions(), trajectoryGuess.back());
        if (fit.getRMSD() <= rmsdThreshold) {
          trajectoryGuess.erase(trajectoryGuess.end());
          if (trajectoryGuess.empty()) {
            _logger->warning
                << "WARNING, no viable trajectory guess provided, all frames are too similar to the endpoints, now "
                   "performing linear interpolation."
                << Core::Log::endl;
            viableTrajectory = false;
            stop = true;
          }
        }
        else {
          stop = true;
        }
      }
      /* Combine input for spline generator or perform linear interpolation if no viable trajectory given */
      if (viableTrajectory) {
        Eigen::MatrixXd interpolationPoints(trajectoryGuess.size() + 2, startVector.size());
        interpolationPoints.row(0) = startVector;
        interpolationPoints.row(trajectoryGuess.size() + 1) = endVector;
        for (int i = 1; i < trajectoryGuess.size() + 1; ++i) {
          interpolationPoints.row(i) = Eigen::Map<const Eigen::VectorXd>(
              trajectoryGuess.at(i - 1).data(), trajectoryGuess.at(i - 1).cols() * trajectoryGuess.at(i - 1).rows());
        }
        Utils::BSplines::FixedEndsPenalizedLeastSquaresGenerator generator(interpolationPoints,
                                                                           numberControlPointsForInterpolation_);
        spline = generator.generateBSpline();
      }
      else {
        spline = Utils::BSplines::LinearInterpolator::generate(startVector, endVector, numberControlPointsForInterpolation_);
      }
    }
    auto molecularSpline = Utils::BSplines::MolecularSpline{start.getElements(), std::move(spline)};
    auto profile = ElementaryStepOptimization::ReactionProfile{std::move(molecularSpline)};

    return profile;
  }

  std::pair<ElementaryStepOptimization::ReactionProfile, bool>
  optimizeElementaryStep(ElementaryStepOptimization::ReactionProfile interpolatedProfile,
                         std::shared_ptr<Core::Calculator> calculator,
                         Utils::UniversalSettings::ValueCollection& taskSettings) const {
    auto cout = _logger->output;
    // Read and delete special settings
    bool stopOnError = stopOnErrorExtraction(taskSettings);
    std::shared_ptr<ElementaryStepOptimization::ElementaryStepOptimizerBase> optimizer;
    if (optimizer_ == "steepestdescent" || optimizer_ == "sd") {
      auto tmp = std::make_shared<ElementaryStepOptimization::ElementaryStepOptimizer<Utils::SteepestDescent>>(
          *calculator, std::move(interpolatedProfile));
      // The original ReaDuct implementation's settings converted into the new form.
      tmp->check.requirement = 4;
      tmp->check.gradRMS = 1.0e-3;
      tmp->check.stepMaxCoeff = 1.0e+10;
      tmp->check.stepRMS = 1.0e+10;
      tmp->check.gradMaxCoeff = 1.0e+10;
      tmp->check.deltaValue = 1.0e+10;
      optimizer = std::move(tmp);
    }
    else if (optimizer_ == "lbfgs") {
      auto tmp = std::make_shared<ElementaryStepOptimization::ElementaryStepOptimizer<Utils::Lbfgs>>(
          *calculator, std::move(interpolatedProfile));
      tmp->optimizer.linesearch = false;
      // The original ReaDuct implementation's settings converted into the new form.
      tmp->check.requirement = 4;
      tmp->check.gradRMS = 1.0e-3;
      tmp->check.stepMaxCoeff = 1.0e+10;
      tmp->check.stepRMS = 1.0e+10;
      tmp->check.gradMaxCoeff = 1.0e+10;
      tmp->check.deltaValue = 1.0e+10;
      optimizer = std::move(tmp);
    }
    else {
      throw std::runtime_error("Optimizer not supported.");
    }

    // Apply user settings
    auto settings = optimizer->getSettings();
    settings.merge(taskSettings);
    if (!settings.valid()) {
      settings.throwIncorrectSettings();
    }
    optimizer->setSettings(settings);
    if (!Utils::GeometryOptimization::settingsMakeSense(*optimizer)) {
      throw std::logic_error("The given calculator settings are too inaccurate for the given convergence criteria of "
                             "this optimization Task");
    }

    const int cycles = optimizer->optimize(*_logger);
    const int maxiter = settings.getInt("convergence_max_iterations");

    bool converged = cycles < maxiter;
    if (converged) {
      cout << Core::Log::endl
           << "    Converged after " << cycles << " iterations." << Core::Log::endl
           << Core::Log::endl;
    }
    else {
      cout << Core::Log::endl
           << "    Stopped after " << maxiter << " iterations." << Core::Log::endl
           << Core::Log::endl;
      if (stopOnError) {
        throw std::runtime_error("Problem: Path optimization did not converge.");
      }
    }
    const auto& optimizedProfile = optimizer->getReactionProfile();

    return std::make_pair(optimizedProfile, converged);
  }

  std::vector<Utils::AtomCollection> extractTransitionStateGuessStructure(ElementaryStepOptimization::ReactionProfile profile,
                                                                          std::shared_ptr<Core::Calculator> calculator,
                                                                          const boost::filesystem::path& dir) const {
    auto cout = _logger->output;
    const auto& molecularSpline = profile.getMolecularSpline();

    // Generate an empty molecule with the correct element types
    const Utils::ElementTypeCollection& ec = molecularSpline.getElements();
    Utils::PositionCollection pc(ec.size(), 3);

    calculator->setStructure(Utils::AtomCollection(ec, pc));
    calculator->setRequiredProperties(Utils::Property::Energy);

    std::vector<double> coordinates;
    std::vector<double> energies;
    int maxEnergyIndex;
    getInitialValues(maxEnergyIndex, coordinates, energies, calculator, profile);

    while (pointDistance(coordinates, maxEnergyIndex) > coordinateThresholdForMaximumExtraction_) {
      coordinates = getNewCoordinates(coordinates, maxEnergyIndex);
      energies = getNewEnergies(energies, coordinates, maxEnergyIndex, calculator, molecularSpline);
      maxEnergyIndex = getIndexForMaxEnergyAndCheckValidity(energies);
    }

    /* Write out tangent */
    if (!tangentFileName_.empty()) {
      Utils::BSplines::BSpline spline = molecularSpline.getBSpline();
      Utils::BSplines::BSpline derivative = spline.getDerivativeBSpline(1);
      Utils::BSplines::MolecularSpline molecularDerivative(ec, derivative);
      Utils::PositionCollection derivativePosition = molecularDerivative.getPositions(coordinates.at(maxEnergyIndex));
      Eigen::Map<Eigen::VectorXd> tangent(derivativePosition.data(), derivativePosition.size());
      writeTangentToFile(tangent, dir);
    }

    std::vector<Utils::AtomCollection> result;
    result.push_back(molecularSpline.at(coordinates.at(maxEnergyIndex)));
    if (extractTsGuessNeighbours_) {
      result.push_back(molecularSpline.at(coordinates.at(maxEnergyIndex - 1)));
      result.push_back(molecularSpline.at(coordinates.at(maxEnergyIndex + 1)));
    }

    return result;
  }

  void writeTangentToFile(const Eigen::VectorXd& tangent, const boost::filesystem::path& generalOutputDir) const {
    // if absolute path given by user, path is directly used
    // else the relative path is used relative to the general output of the calculation
    boost::filesystem::path tangentPath(tangentFileName_);
    if (tangentFileName_.at(0) != '/') {
      tangentPath = boost::filesystem::absolute(generalOutputDir / tangentFileName_);
    }
    // create necessary directories for file if necessary
    boost::filesystem::path parent = tangentPath.parent_path();
    if (!parent.empty()) {
      // does not give error if directories already exist
      boost::filesystem::detail::create_directories(parent);
    }
    std::ofstream fout(tangentPath.string());
    if (!fout.is_open()) {
      throw std::runtime_error("Problem when opening/creating file: " + tangentPath.string());
    }
    fout.imbue(std::locale("C"));
    for (int i = 0; i < tangent.size(); ++i) {
      fout << tangent.row(i) << '\n';
    }
    fout.close();
  };

  Utils::MolecularTrajectory convertProfileToTrajectory(ElementaryStepOptimization::ReactionProfile profile) const {
    auto elements = profile.getMolecularSpline().getElements();
    auto spline = profile.getMolecularSpline().getBSpline();
    Utils::MolecularTrajectory trajectory = discretizeSpline(elements, spline, numberStructuresForMolecularTrajectory_);
    return trajectory;
  }

  static Utils::MolecularTrajectory discretizeSpline(const Utils::ElementTypeCollection& elements,
                                                     const Utils::BSplines::BSpline& spline, int numberPoints) {
    Utils::MolecularTrajectory t;
    t.setElementTypes(elements);
    double deltaU = 1.0 / (numberPoints - 1);
    for (int i = 0; i < numberPoints; ++i) {
      double u = i * deltaU;
      Utils::PositionCollection pc =
          Eigen::Map<const Utils::PositionCollection>(spline.evaluate(u).data(), elements.size(), 3);
      t.push_back(std::move(pc));
    }
    return t;
  }

  static void writeTrajectory(const Utils::MolecularTrajectory& trajectory, const std::string& filepath) {
    std::ofstream ostream(filepath, std::ofstream::out);
    Utils::MolecularTrajectoryIO::write(Utils::MolecularTrajectoryIO::format::xyz, ostream, trajectory);
    ostream.close();
  }

  static Utils::MolecularTrajectory readTrajectory(const std::string& filepath) {
    Utils::MolecularTrajectory trajectory;
    std::filebuf fb;
    if (fb.open(filepath, std::istream::in)) {
      std::istream is(&fb);
      trajectory = Utils::MolecularTrajectoryIO::read(Utils::MolecularTrajectoryIO::format::xyz, is);
      fb.close();
    }
    return trajectory;
  }

  static void getInitialValues(int& maxEnergyIndex, std::vector<double>& coordinates, std::vector<double>& energies,
                               std::shared_ptr<Core::Calculator> calculator,
                               const ElementaryStepOptimization::ReactionProfile& profile) {
    const auto& profileEnergies = profile.getProfileEnergies();
    if (!profile.getProfileEnergies().empty()) {
      coordinates = profileEnergies.getCoordinates();
      energies = profileEnergies.getEnergies();
      maxEnergyIndex = getIndexForMaxEnergyAndCheckValidity(energies);
    }
    else {
      const auto& molecularSpline = profile.getMolecularSpline();
      coordinates = {0.0, 0.5, 1.0};
      energies.resize(3);

      calculator->modifyPositions(molecularSpline.getPositions(coordinates[0]));
      energies[0] = calculator->calculate("").get<Utils::Property::Energy>();

      calculator->modifyPositions(molecularSpline.getPositions(coordinates[1]));
      energies[1] = calculator->calculate("").get<Utils::Property::Energy>();

      calculator->modifyPositions(molecularSpline.getPositions(coordinates[2]));
      energies[2] = calculator->calculate("").get<Utils::Property::Energy>();

      maxEnergyIndex = 1;
    }
  }

  static double pointDistance(const std::vector<double>& coordinates, int maxEnergyIndex) {
    double diff = coordinates[maxEnergyIndex + 1] - coordinates[maxEnergyIndex - 1];
    return diff / 2;
  }

  static std::vector<double> getNewCoordinates(const std::vector<double>& oldCoordinates, int maxEnergyIndex) {
    std::vector<double> u(5);
    u[0] = oldCoordinates[maxEnergyIndex - 1];
    u[2] = oldCoordinates[maxEnergyIndex];
    u[4] = oldCoordinates[maxEnergyIndex + 1];
    u[1] = 0.5 * (u[0] + u[2]);
    u[3] = 0.5 * (u[2] + u[4]);
    return u;
  }

  static std::vector<double> getNewEnergies(const std::vector<double>& oldEnergies, const std::vector<double>& coordinates,
                                            int maxEnergyIndex, std::shared_ptr<Core::Calculator> calculator,
                                            const Utils::BSplines::MolecularSpline& spline) {
    std::vector<double> e(5);
    e[0] = oldEnergies[maxEnergyIndex - 1];
    e[2] = oldEnergies[maxEnergyIndex];
    e[4] = oldEnergies[maxEnergyIndex + 1];

    calculator->modifyPositions(spline.getPositions(coordinates[1]));
    e[1] = calculator->calculate("").get<Utils::Property::Energy>();

    calculator->modifyPositions(spline.getPositions(coordinates[3]));
    e[3] = calculator->calculate("").get<Utils::Property::Energy>();

    return e;
  }

  static int getIndexForMaxEnergyAndCheckValidity(const std::vector<double>& energies) {
    auto maxIt = std::max_element(energies.begin(), energies.end());
    auto dist = std::distance(energies.begin(), maxIt);
    if (dist == 0 || dist == static_cast<int>(energies.size()) - 1) {
      throw std::runtime_error("Problem: End point has maximal energy.");
    }
    return static_cast<int>(dist);
  }

  // Default settings
  mutable bool alignStructuresBeforeInterpolation_ = true;
  mutable int numberControlPointsForInterpolation_ = 5;
  mutable int numberStructuresForMolecularTrajectory_ = 10;
  mutable bool optimize_ = true;
  mutable std::string optimizer_ = "lbfgs";
  mutable bool extractTsGuess_ = false;
  mutable bool extractTsGuessNeighbours_ = false;
  mutable std::string tangentFileName_;
  mutable double coordinateThresholdForMaximumExtraction_ = 1e-3;
};

} // namespace Readuct
} // namespace Scine

#endif // READUCT_BSPLINEINTERPOLATIONTASK_H_
