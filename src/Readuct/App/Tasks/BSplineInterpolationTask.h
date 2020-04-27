/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef READUCT_BSPLINEINTERPOLATIONTASK_H_
#define READUCT_BSPLINEINTERPOLATIONTASK_H_

#include "../../Readuct/ElementaryStepOptimization/ElementaryStepOptimizer.h"
#include "../../Readuct/ElementaryStepOptimization/ReactionProfile.h"
#include "Tasks/Task.h"
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/MolecularTrajectoryIO.h>
#include <Utils/IO/Yaml.h>
#include <Utils/Math/BSplines/FixedEndsPenalizedLeastSquaresGenerator.h>
#include <Utils/Math/BSplines/InterpolationGenerator.h>
#include <Utils/Math/BSplines/LinearInterpolator.h>
#include <Utils/Math/BSplines/MolecularSpline.h>
#include <Utils/Math/QuaternionFit.h>
#include <Utils/MolecularTrajectory.h>
#include <Utils/Optimizer/GradientBased/Lbfgs.h>
#include <Utils/Optimizer/GradientBased/SteepestDescent.h>
#include <fstream>

namespace Scine {
namespace Readuct {

class BSplineInterpolationTask : public Task {
 public:
  BSplineInterpolationTask(std::vector<std::string> input, std::vector<std::string> output) : Task(input, output) {
  }

  std::string name() const override {
    return "BSpline Interpolation";
  }

  virtual bool run(std::map<std::string, std::shared_ptr<Core::Calculator>>& systems, const YAML::Node& taskSettings) const final {
    // Get/Copy Calculator
    std::shared_ptr<Core::Calculator> calc;
    if (systems.find(_input[0]) != systems.end()) {
      calc = std::shared_ptr<Core::Calculator>(systems.at(_input[0])->clone().release());
    }
    else {
      throw std::runtime_error("Missing system '" + _input[0] + "' in BSpline interpolation.");
    }

    // Get inputs
    auto start = systems.at(_input[0])->getStructure();
    auto end = systems.at(_input[1])->getStructure();

    if (start->getElements() != end->getElements()) {
      throw std::runtime_error("The given structures do not have the same element types.");
    }

    Utils::MolecularTrajectory trajectoryGuess = Utils::MolecularTrajectory();
    // Read and set user-specified settings
    auto settingsCopy = taskSettings;
    if (auto type = settingsCopy["align_structures"]) {
      alignStructuresBeforeInterpolation_ = type.as<bool>();
      settingsCopy.remove("align_structures");
    }
    if (auto type = settingsCopy["num_structures"]) {
      numberStructuresForMolecularTrajectory_ = type.as<int>();
      settingsCopy.remove("num_structures");
    }
    if (auto type = settingsCopy["optimize"]) {
      optimize_ = type.as<bool>();
      settingsCopy.remove("optimize");
    }
    if (auto type = settingsCopy["num_control_points"]) {
      numberControlPointsForInterpolation_ = type.as<int>();
      settingsCopy.remove("num_control_points");
    }
    if (auto type = settingsCopy["optimizer"]) {
      optimizer_ = type.as<std::string>();
      settingsCopy.remove("optimizer");
    }
    if (auto type = settingsCopy["extract_ts_guess"]) {
      extractTsGuess_ = type.as<bool>();
      settingsCopy.remove("extract_ts_guess");
    }
    if (auto type = settingsCopy["extract_ts_guess_neighbours"]) {
      extractTsGuessNeighbours_ = type.as<bool>();
      settingsCopy.remove("extract_ts_guess_neighbours");
    }
    if (auto type = settingsCopy["tangent_file"]) {
      tangentFileName_ = type.as<std::string>();
      extractTangent_ = true;
      settingsCopy.remove("tangent_file");
    }
    if (auto type = settingsCopy["extract_threshold"]) {
      coordinateThresholdForMaximumExtraction_ = type.as<double>();
      settingsCopy.remove("extract_threshold");
    }
    if (auto type = settingsCopy["trajectory_guess"]) {
      std::vector<std::string> paths = type.as<std::vector<std::string>>();
      for (int i = 0; i < paths.size(); ++i) {
        Utils::MolecularTrajectory traj = readTrajectory(paths.at(i));
        for (int j = 0; j < traj.size(); ++j) {
          trajectoryGuess.push_back(traj.at(j));
        }
      }
      settingsCopy.remove("trajectory_guess");
    }

    // Interpolate between start and end structure
    printf("  Interpolating Reaction Path\n");
    ElementaryStepOptimization::ReactionProfile interpolatedProfile = interpolateElementaryStep(*start, *end, trajectoryGuess);
    Utils::MolecularTrajectory interpolatedTrajectory = convertProfileToTrajectory(interpolatedProfile);
    // Print/Store results
    Utils::XyzStreamHandler writer;
    boost::filesystem::path dir(((_output.size() > 0) ? _output[0] : _input[0]));
    boost::filesystem::create_directory(dir);
    if (_output.size() > 0) {
      systems[_output[0]] = calc;
      boost::filesystem::path xyzfile(_output[0] + "_interpolated.xyz");
      writeTrajectory(interpolatedTrajectory, (dir / xyzfile).string());
    }
    else {
      systems[_input[0]] = calc;
      boost::filesystem::path xyzfile(_input[0] + "_interpolated.xyz");
      writeTrajectory(interpolatedTrajectory, (dir / xyzfile).string());
    }
    printf("  Interpolating Complete\n\n");

    // Optimize the interpolated path
    Utils::MolecularTrajectory optimizedTrajectory;
    ElementaryStepOptimization::ReactionProfile optimizedProfile;
    if (optimize_) {
      printf("  Optimizing Reaction Path\n");
      optimizedProfile = optimizeElementaryStep(interpolatedProfile, calc, settingsCopy);
      optimizedTrajectory = convertProfileToTrajectory(optimizedProfile);
      // Output
      if (_output.size() > 0) {
        boost::filesystem::path xyzfile2(_output[0] + "_optimized.xyz");
        writeTrajectory(optimizedTrajectory, (dir / xyzfile2).string());
      }
      else {
        boost::filesystem::path xyzfile2(_input[0] + "_optimized.xyz");
        writeTrajectory(optimizedTrajectory, (dir / xyzfile2).string());
      }
      printf("  Optimization Complete\n\n");
    }

    std::vector<Utils::AtomCollection> tsGuess;
    if (optimize_ && extractTsGuess_) {
      printf("  Extracting Transition State Guess\n");
      tsGuess = extractTransitionStateGuessStructure(optimizedProfile, calc, dir);
      // Output
      if (_output.size() > 0) {
        boost::filesystem::path xyzfile3(_output[0] + "_tsguess.xyz");
        std::ofstream xyz((dir / xyzfile3).string(), std::ofstream::out);
        writer.write(xyz, tsGuess.at(0));
        xyz.close();
        if (extractTsGuessNeighbours_) {
          boost::filesystem::path xyzfile3_1(_output[0] + "_tsguess-1.xyz");
          std::ofstream xyz_1((dir / xyzfile3_1).string(), std::ofstream::out);
          writer.write(xyz_1, tsGuess.at(1));
          xyz_1.close();
          boost::filesystem::path xyzfile3_2(_output[0] + "_tsguess+1.xyz");
          std::ofstream xyz_2((dir / xyzfile3_2).string(), std::ofstream::out);
          writer.write(xyz_2, tsGuess.at(2));
          xyz_2.close();
        }
      }
      else {
        boost::filesystem::path xyzfile3(_input[0] + "_tsguess.xyz");
        std::ofstream xyz((dir / xyzfile3).string(), std::ofstream::out);
        writer.write(xyz, tsGuess.at(0));
        xyz.close();
        if (extractTsGuessNeighbours_) {
          boost::filesystem::path xyzfile3_1(_input[0] + "_tsguess-1.xyz");
          std::ofstream xyz_1((dir / xyzfile3_1).string(), std::ofstream::out);
          writer.write(xyz_1, tsGuess.at(1));
          xyz_1.close();
          boost::filesystem::path xyzfile3_2(_input[0] + "_tsguess+1.xyz");
          std::ofstream xyz_2((dir / xyzfile3_2).string(), std::ofstream::out);
          writer.write(xyz_2, tsGuess.at(2));
          xyz_2.close();
        }
      }
      printf("  Extraction Complete\n\n");
    }
    printf("\n");
    return true;
  }

 private:
  ElementaryStepOptimization::ReactionProfile
  interpolateElementaryStep(const Utils::AtomCollection& start, Utils::AtomCollection end,
                            Utils::MolecularTrajectory trajectoryGuess = Utils::MolecularTrajectory()) const {
    if (alignStructuresBeforeInterpolation_) {
      auto positionsToAlign = end.getPositions();
      Utils::Geometry::alignPositions(start.getPositions(), positionsToAlign);
      end.setPositions(std::move(positionsToAlign));
    }

    Eigen::VectorXd startVector = Eigen::Map<const Eigen::VectorXd>(
        start.getPositions().data(), start.getPositions().cols() * start.getPositions().rows());
    Eigen::VectorXd endVector =
        Eigen::Map<const Eigen::VectorXd>(end.getPositions().data(), end.getPositions().cols() * end.getPositions().rows());
    Utils::BSplines::BSpline spline;
    if (trajectoryGuess.empty())
      spline = Utils::BSplines::LinearInterpolator::generate(startVector, endVector, numberControlPointsForInterpolation_);
    else {
      assert(start.size() == trajectoryGuess.molecularSize());
      /* Remove frames from trajectory at end and beginning which are too similar to given endpoints */
      bool viableTrajectory = true;
      double rmsdThreshold = 1.0;
      bool stop = false;
      while (!stop) {
        Utils::QuaternionFit fit(start.getPositions(), trajectoryGuess.front());
        if (fit.getRMSD() <= rmsdThreshold) {
          trajectoryGuess.erase(trajectoryGuess.begin());
          if (trajectoryGuess.empty()) {
            std::cout << "ERROR, no viable trajectory guess provided, all frames are too similar to the endpoints, now "
                         "performing linear interpolation"
                      << std::endl;
            viableTrajectory = false;
            stop = true;
          }
        }
        else {
          stop = true;
        }
      }
      if (viableTrajectory)
        stop = false;
      while (!stop) {
        Utils::QuaternionFit fit(end.getPositions(), trajectoryGuess.back());
        if (fit.getRMSD() <= rmsdThreshold) {
          trajectoryGuess.erase(trajectoryGuess.end());
          if (trajectoryGuess.empty()) {
            std::cout << "ERROR, no viable trajectory guess provided, all frames are too similar to the endpoints, now "
                         "performing linear interpolation"
                      << std::endl;
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

  ElementaryStepOptimization::ReactionProfile
  optimizeElementaryStep(ElementaryStepOptimization::ReactionProfile interpolatedProfile,
                         std::shared_ptr<Core::Calculator> calculator, YAML::Node& settingsCopy) const {
    std::unique_ptr<ElementaryStepOptimization::ElementaryStepOptimizerBase> optimizer;
    if (optimizer_ == "steepestdescent" || optimizer_ == "sd") {
      auto tmp = std::make_unique<ElementaryStepOptimization::ElementaryStepOptimizer<Utils::SteepestDescent>>(
          *calculator, interpolatedProfile);
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
      auto tmp = std::make_unique<ElementaryStepOptimization::ElementaryStepOptimizer<Utils::Lbfgs>>(*calculator,
                                                                                                     interpolatedProfile);
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
    Scine::Utils::nodeToSettings(settings, settingsCopy);
    optimizer->setSettings(settings);

    int cycles = optimizer->optimize();
    int maxiter = settings.getInt("convergence_max_iterations");
    if (cycles >= maxiter)
      throw std::runtime_error("Problem: Path optimization did not converge.");
    printf("  Converged path after %d optimization cycles.\n", cycles);
    const auto& optimizedProfile = optimizer->getReactionProfile();

    return optimizedProfile;
  }

  std::vector<Utils::AtomCollection> extractTransitionStateGuessStructure(ElementaryStepOptimization::ReactionProfile profile,
                                                                          std::shared_ptr<Core::Calculator> calculator,
                                                                          boost::filesystem::path dir) const {
    const auto& molecularSpline = profile.getMolecularSpline();

    // Generate an empty molecule with the correct element types
    Utils::ElementTypeCollection ec = molecularSpline.getElements();
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
    if (extractTangent_) {
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
    if (tangentFileName_.at(0) != '/')
      tangentPath = boost::filesystem::absolute(generalOutputDir / tangentFileName_);
    // create necessary directories for file if necessary
    boost::filesystem::path parent = tangentPath.parent_path();
    if (!parent.empty())
      // does not give error if directories already exist
      boost::filesystem::detail::create_directories(parent);
    std::ofstream fout(tangentPath.string());
    if (!fout.is_open())
      throw std::runtime_error("Problem when opening/creating file: " + tangentPath.string());
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

  Utils::MolecularTrajectory discretizeSpline(const Utils::ElementTypeCollection& elements,
                                              const Utils::BSplines::BSpline& spline, int numberPoints) const {
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

  void writeTrajectory(Utils::MolecularTrajectory trajectory, std::string filepath) const {
    std::ofstream ostream(filepath, std::ofstream::out);
    Utils::MolecularTrajectoryIO::write(Utils::MolecularTrajectoryIO::format::xyz, ostream, trajectory);
    ostream.close();
  }

  Utils::MolecularTrajectory readTrajectory(std::string filepath) const {
    Utils::MolecularTrajectory trajectory;
    std::filebuf fb;
    if (fb.open(filepath, std::istream::in)) {
      std::istream is(&fb);
      trajectory = Utils::MolecularTrajectoryIO::read(Utils::MolecularTrajectoryIO::format::xyz, is);
      fb.close();
    }
    return trajectory;
  }

  void getInitialValues(int& maxEnergyIndex, std::vector<double>& coordinates, std::vector<double>& energies,
                        std::shared_ptr<Core::Calculator> calculator,
                        const ElementaryStepOptimization::ReactionProfile& profile) const {
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

  double pointDistance(const std::vector<double>& coordinates, int maxEnergyIndex) const {
    double diff = coordinates[maxEnergyIndex + 1] - coordinates[maxEnergyIndex - 1];
    return diff / 2;
  }

  std::vector<double> getNewCoordinates(const std::vector<double>& oldCoordinates, int maxEnergyIndex) const {
    std::vector<double> u(5);
    u[0] = oldCoordinates[maxEnergyIndex - 1];
    u[2] = oldCoordinates[maxEnergyIndex];
    u[4] = oldCoordinates[maxEnergyIndex + 1];
    u[1] = 0.5 * (u[0] + u[2]);
    u[3] = 0.5 * (u[2] + u[4]);
    return u;
  }

  std::vector<double> getNewEnergies(const std::vector<double>& oldEnergies, const std::vector<double>& coordinates,
                                     int maxEnergyIndex, std::shared_ptr<Core::Calculator> calculator,
                                     const Utils::BSplines::MolecularSpline& spline) const {
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

  int getIndexForMaxEnergyAndCheckValidity(const std::vector<double>& energies) const {
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
  mutable bool extractTangent_ = false;
  mutable std::string tangentFileName_ = "";
  mutable double coordinateThresholdForMaximumExtraction_ = 1e-3;
};

} // namespace Readuct
} // namespace Scine

#endif // READUCT_BSPLINEINTERPOLATIONTASK_H_
