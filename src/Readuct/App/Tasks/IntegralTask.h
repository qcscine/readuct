/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef READUCT_INTEGRALSTASK_H
#define READUCT_INTEGRALSTASK_H

/* Readuct */
#include "Tasks/Task.h"
/* Scine */
#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>
/* External */
#include "boost/exception/diagnostic_information.hpp"
#include <boost/filesystem.hpp>
/* std c++ */
#include <iomanip>
#include <string>
#include <vector>

namespace Scine::Readuct {

/**
 * @class IntegralTask IntegralTask.h
 * @brief This tasks writes integrals, such as the one-electron integrals, to file.
 */
class IntegralTask : public Task {
 public:
  /**
   * @brief Construct a new IntegralTask.
   * @param input  The input system names for the task.
   * @param output The output system names for the task.
   * @param logger The logger to/through which all text output will be handled.
   */
  IntegralTask(std::vector<std::string> input, std::vector<std::string> output, std::shared_ptr<Core::Log> logger = nullptr)
    : Task(std::move(input), std::move(output), std::move(logger)) {
  }
  /**
   * @brief Getter for the task's name.
   * @return The task's name.
   */
  std::string name() const override {
    return "Integral Calculation";
  }
  /**
   * @brief Run the task.
   * @param systems       The input systems.
   * @param taskSettings  The task settings.
   * @param testRunOnly   If true, the main part of the task is not executed. Only used for testing.
   * @param observers     Optional observers. This task does not support observers at the moment.
   * @return Returns true if the task terminated successfully. False, otherwise.
   */
  bool run(SystemsMap& systems, Utils::UniversalSettings::ValueCollection taskSettings, bool testRunOnly = false,
           std::vector<std::function<void(const int&, const Utils::AtomCollection&, const Utils::Results&, const std::string&)>>
               observers = {}) const final {
    // Read and delete special settings
    const bool stopOnError = stopOnErrorExtraction(taskSettings);
    const bool requireHCore = taskSettings.extract("require_one_electron_integrals", true);
    const bool silentCalculator = taskSettings.extract("silent_stdout_calculator", true);
    if (!taskSettings.empty()) {
      std::string keyListing = "\n";
      auto keys = taskSettings.getKeys();
      for (const auto& key : keys) {
        keyListing += "'" + key + "'\n";
      }
      throw std::logic_error("Specified one or more task settings that are not available for this task:" + keyListing);
    }
    if (observers.size() > 0) {
      throw std::logic_error("IntegralTask does not feature algorithm accepting observers, yet one was given");
    }
    if (testRunOnly) {
      return true; // leave out rest in case of task chaining
    }

    // Note: _input is guaranteed not to be empty by Task constructor
    auto calc = copyCalculator(systems, _input.front(), name());
    const auto previousResults = calc->results();
    Utils::CalculationRoutines::setLog(*calc, true, true, !silentCalculator);
    // Check for available properties
    const bool hCoreAvailable = calc->possibleProperties().containsSubSet(Utils::Property::OneElectronMatrix);

    if (requireHCore && !hCoreAvailable) {
      throw std::logic_error(
          "One-electron integrals required, but chosen calculator does not provide them.\n"
          "If you do not need one-electron integrals, set 'require_one_electron_integrals' to 'false'"
          " in the task settings");
    }
    Utils::PropertyList requiredProperties;
    if (requireHCore) {
      requiredProperties.addProperty(Utils::Property::OneElectronMatrix);
    }

    try {
      calc->setRequiredProperties(requiredProperties);
      calc->calculate(name());
      if (!calc->results().get<Utils::Property::SuccessfulCalculation>()) {
        throw std::runtime_error(name() + " was not successful");
      }
    }
    catch (...) {
      if (stopOnError) {
        throw;
      }
      _logger->error
          << "  " + name() + " was not successful with error:\n  " + boost::current_exception_diagnostic_information()
          << Core::Log::endl;
      calc->results() = previousResults + calc->results();
      return false;
    }

    if (requireHCore) {
      const std::string& outputSystem = (!_output.empty() ? _output.front() : _input.front());
      const boost::filesystem::path dir(outputSystem);
      boost::filesystem::create_directory(dir);
      const boost::filesystem::path asciiFile(outputSystem + ".hcore.dat");
      const std::string outPath = (dir / asciiFile).string();
      _logger->output << "Writing core Hamiltonian integrals to file " << outPath << Core::Log::endl;
      std::ofstream outFile(outPath, std::ofstream::out);
      const auto hCoreMatrix = calc->results().template get<Utils::Property::OneElectronMatrix>();
      outFile << "AO core Hamiltonian integrals.\n";
      outFile << std::scientific << std::setprecision(12) << hCoreMatrix << "\n";
      outFile.close();

      // Store result
      if (!_output.empty()) {
        systems[_output[0]] = calc;
      }
      else {
        systems[_input[0]] = calc;
      }
    }

    return true;
  }
};

} // namespace Scine::Readuct

#endif // READUCT_INTEGRALSTASK_H
