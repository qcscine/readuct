/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef READUCT_HESSIANTASK_H_
#define READUCT_HESSIANTASK_H_

/* Readuct */
#include "Tasks/Task.h"
/* Scine */
#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/GeometricDerivatives/NormalModeAnalysis.h>
#include <Utils/GeometricDerivatives/NormalModesContainer.h>
#include <Utils/IO/FormattedString.h>
#include <Utils/IO/MolecularTrajectoryIO.h>
#include <Utils/Properties/Thermochemistry/ThermochemistryCalculator.h>
/* External */
#include "boost/exception/diagnostic_information.hpp"
#include <boost/filesystem.hpp>
/* std c++ */
#include <cstdio>
#include <fstream>
#include <iostream>

namespace Scine {
namespace Readuct {

class HessianTask : public Task {
 public:
  /**
   * @brief Construct a new HessianTask.
   * @param input  The input system names for the task.
   * @param output The output system names for the task.
   * @param logger The logger to/through which all text output will be handled.
   */
  HessianTask(std::vector<std::string> input, std::vector<std::string> output, std::shared_ptr<Core::Log> logger = nullptr)
    : Task(std::move(input), std::move(output), std::move(logger)) {
  }

  std::string name() const override {
    return "Hessian Calculation";
  }

  bool run(SystemsMap& systems, Utils::UniversalSettings::ValueCollection taskSettings, bool testRunOnly = false) const final {
    warningIfMultipleInputsGiven();
    warningIfMultipleOutputsGiven();

    // Read and delete special settings
    bool stopOnError = stopOnErrorExtraction(taskSettings);
    if (!taskSettings.empty()) {
      throw std::logic_error(falseTaskSettingsErrorMessage(name()));
    }

    // If no errors encountered until here, the basic settings should be alright
    if (testRunOnly) {
      return true;
    }

    // Note: _input is guaranteed not to be empty by Task constructor
    auto calc = copyCalculator(systems, _input.front(), name());

    // Check system size
    if (calc->getStructure()->size() == 1) {
      throw std::runtime_error("Cannot perform Hessian task for monoatomic systems.");
    }

    // Get/calculate Hessian
    if (!calc->results().has<Utils::Property::Hessian>() || !calc->results().has<Utils::Property::Thermochemistry>() ||
        !calc->results().has<Utils::Property::Energy>()) {
      calc->setRequiredProperties(Utils::Property::Hessian | Utils::Property::Thermochemistry | Utils::Property::Energy);
      try {
        calc->calculate(name());
        if (!calc->results().get<Utils::Property::SuccessfulCalculation>()) {
          throw std::runtime_error(name() + " was not successful");
        }
      }
      catch (...) {
        if (stopOnError) {
          throw;
        }
        _logger->warning
            << "  " + name() + " was not successful with error:\n  " + boost::current_exception_diagnostic_information()
            << Core::Log::endl;
        return false;
      }
    }
    auto hessian = calc->results().get<Utils::Property::Hessian>();

    // Print frequencies
    auto system = calc->getStructure();
    Utils::NormalModesContainer modes =
        Utils::NormalModeAnalysis::calculateNormalModes(hessian, system->getElements(), system->getPositions());
    auto wavenumbers = modes.getWaveNumbers();
    auto cout = _logger->output;
    cout << "  Vib. Frequencies:\n";
    cout << "  [Rot. and trans. vib. removed, imaginary freq. shown as negative.]\n\n";
    cout.printf("  %4s %8s\n", "#", "cm^-1");
    for (unsigned int i = 0; i < wavenumbers.size(); i++) {
      cout.printf("  %4d %+8.1f\n", i + 1, wavenumbers[i]);
    }
    cout << Core::Log::endl;

    const std::string& outputSystem = (!_output.empty() ? _output[0] : _input[0]);

    // Setup output directory if needed
    boost::filesystem::path dir(outputSystem);
    if (wavenumbers[0] < 0.0) {
      boost::filesystem::create_directory(dir);
    }
    // Print Imaginary modes
    for (unsigned int i = 0; i < wavenumbers.size() && wavenumbers[i] < 0.0; ++i) {
      auto trj = modes.getModeAsMolecularTrajectory(i, *system, 1.0);
      std::string filename = Utils::format("%s.vibmode-%05d.trj.xyz", outputSystem.c_str(), i + 1);
      std::ofstream xyz((dir / filename).string(), std::ofstream::out);
      Utils::MolecularTrajectoryIO::write(Utils::MolecularTrajectoryIO::format::xyz, xyz, trj);
      xyz.close();
    }

    // Store result
    systems[outputSystem] = calc;

    // Print thermochemistry results
    auto thermo = calc->results().get<Utils::Property::Thermochemistry>();
    double temp = calc->settings().getDouble("temperature");
    cout << "  Thermochemical Data:\n";
    cout.printf("  %35s %+14.9f K\n", "Temperature (T):", temp);
    cout.printf("  %35s %14d\n", "Molecular symmetry number:", thermo.overall.symmetryNumber);
    cout.printf("  %35s %+14.9f E_h\n", "Zero Point Vibrational Energy:", thermo.overall.zeroPointVibrationalEnergy);
    cout.printf("  %35s %+14.9f E_h/K\n", "Entropy (S):", thermo.overall.entropy);
    cout.printf("  %35s %+14.9f E_h\n", "Entropy (-TS):", -temp * thermo.overall.entropy);
    cout.printf("  %35s %+14.9f E_h\n", "Enthalpy (H):", thermo.overall.enthalpy);
    double electronicEnergy = calc->results().get<Utils::Property::Energy>();
    cout.printf("  %35s %+14.9f E_h\n", "Gibbs Energy Correction (dG):", thermo.overall.gibbsFreeEnergy - electronicEnergy);
    cout.printf("  %35s %+14.9f E_h\n", "Electronic Energy (E):", electronicEnergy);
    cout.printf("  %35s %+14.9f E_h\n", "Gibbs Energy (E + dG):", thermo.overall.gibbsFreeEnergy);
    cout << Core::Log::endl << Core::Log::endl;

    return true;
  }
};

} // namespace Readuct
} // namespace Scine

#endif // READUCT_HESSIANTASK_H_
