/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef READUCT_SINGLEPOINTTASK_H_
#define READUCT_SINGLEPOINTTASK_H_

/* Readuct */
#include "Tasks/Task.h"
/* Scine */
#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/IO/FormattedIOUtils.h>
#include <Utils/IO/Yaml.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>
/* External */
#include "boost/exception/diagnostic_information.hpp"
#include <boost/filesystem.hpp>
/* std c++ */
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

namespace Scine {
namespace Readuct {

class SinglePointTask : public Task {
 public:
  /**
   * @brief Construct a new SinglePointTask.
   * @param input  The input system names for the task.
   * @param output The output system names for the task.
   * @param logger The logger to/through which all text output will be handled.
   */
  SinglePointTask(std::vector<std::string> input, std::vector<std::string> output, std::shared_ptr<Core::Log> logger = nullptr)
    : Task(std::move(input), std::move(output), std::move(logger)) {
  }

  std::string name() const override {
    return "Single Point Calculation";
  }

  bool run(SystemsMap& systems, Utils::UniversalSettings::ValueCollection taskSettings, bool testRunOnly = false) const final {
    warningIfMultipleInputsGiven();
    warningIfMultipleOutputsGiven();

    // Read and delete special settings
    bool stopOnError = stopOnErrorExtraction(taskSettings);
    bool requireCharges = taskSettings.extract("require_charges", true);
    bool requireGradients = taskSettings.extract("require_gradients", false);
    bool requireOrbitalEnergies = taskSettings.extract("orbital_energies", false);
    if (!taskSettings.empty()) {
      throw std::logic_error("Specified a task setting that is not available for this task");
    }

    if (testRunOnly) {
      return true; // leave out rest in case of task chaining
    }

    // Note: _input is guaranteed not to be empty by Task constructor
    auto calc = copyCalculator(systems, _input.front(), name());

    // Check for available properties
    bool chargesAvailable = calc->possibleProperties().containsSubSet(Utils::Property::AtomicCharges);
    bool gradientsAvailable = calc->possibleProperties().containsSubSet(Utils::Property::Gradients);
    bool orbitalEnergiesAvailable = calc->possibleProperties().containsSubSet(Utils::Property::OrbitalEnergies);
    if (requireCharges && !chargesAvailable) {
      throw std::logic_error("Charges required, but chosen calculator does not provide them.\n"
                             "If you do not need charges, set 'require_charges' to 'false' in the task settings");
    }
    if (requireGradients && !gradientsAvailable) {
      throw std::logic_error("Gradients required, but chosen calculator does not provide them.");
    }
    if (requireOrbitalEnergies && !orbitalEnergiesAvailable) {
      throw std::logic_error("Orbital energies required, but chosen calculator does not provide them.");
    }

    // Calculate energy, and possibly atomic charges and/or gradients
    Utils::PropertyList requiredProperties = Utils::Property::Energy;
    if (requireCharges) {
      requiredProperties.addProperty(Utils::Property::AtomicCharges);
    }

    if (requireGradients) {
      requiredProperties.addProperty(Utils::Property::Gradients);
    }
    if (requireOrbitalEnergies) {
      requiredProperties.addProperty(Utils::Property::OrbitalEnergies);
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
      _logger->warning
          << "  " + name() + " was not successful with error:\n  " + boost::current_exception_diagnostic_information()
          << Core::Log::endl;
      return false;
    }

    // Store result
    if (!_output.empty()) {
      systems[_output[0]] = calc;
    }
    else {
      systems[_input[0]] = calc;
    }

    auto energy = calc->results().get<Utils::Property::Energy>();

    // Print energy
    auto cout = _logger->output;
    cout.printf("  The (electronic) energy is: %+16.9f hartree\n\n", energy);
    // Print electronic temperature
    if (calc->settings().valueExists(Utils::SettingsNames::electronicTemperature)) {
      auto etemp = calc->settings().getDouble(Utils::SettingsNames::electronicTemperature);
      cout.printf("  The (electronic) temperature was: %+10.3f K\n\n", etemp);
    }
    // Print charges
    if (requireCharges) {
      auto charges = calc->results().get<Utils::Property::AtomicCharges>();
      auto atomColl = calc->getStructure();
      auto elements = atomColl->getElements();
      cout << "  Atomic Partial Charges:\n\n";
      for (size_t i = 0; i < charges.size(); i++) {
        cout.printf("  %5d %-6s: %+11.6f\n", i, Utils::ElementInfo::symbol(elements[i]).c_str(), charges[i]);
      }
    }
    cout << Core::Log::nl;
    // Print gradients
    if (requireGradients) {
      auto gradients = calc->results().get<Utils::Property::Gradients>();
      cout << "  Gradients (hartree / bohr):\n\n";
      cout << [&gradients](std::ostream& os) { Utils::matrixPrettyPrint(os, gradients); };
    }
    // Print orbital energies
    if (requireOrbitalEnergies) {
      Utils::SingleParticleEnergies orbitalEnergies = calc->results().get<Utils::Property::OrbitalEnergies>();

      cout << "  Orbital Energies:\n\n";

      if (orbitalEnergies.isRestricted()) {
        const auto restrictedEnergies = orbitalEnergies.getRestrictedEnergies();

        std::cout << std::setw(20) << "Index" << std::setw(20) << "Energy / Hartree" << std::endl;
        for (unsigned int i = 0; i < restrictedEnergies.size(); ++i) {
          std::cout << std::setprecision(7) << std::scientific << std::right << std::setw(5) << i << std::setw(5)
                    << restrictedEnergies[i] << std::endl;
        }
      }
      else {
        const auto alphaEnergies = orbitalEnergies.getAlphaEnergies();
        const auto betaEnergies = orbitalEnergies.getBetaEnergies();

        cout << "  Alpha spin:\n\n";
        std::cout << std::setw(5) << "  Index" << std::setw(20) << "Energy / Hartree" << std::endl;
        for (unsigned int i = 0; i < alphaEnergies.size(); ++i) {
          std::cout << std::setprecision(7) << std::scientific << std::right << std::setw(5) << i << std::setw(22)
                    << alphaEnergies[i] << std::endl;
        }

        cout << "\n\n  Beta spin:\n\n";
        std::cout << std::setw(5) << "  Index" << std::setw(20) << "Energy / Hartree" << std::endl;
        for (unsigned int i = 0; i < betaEnergies.size(); ++i) {
          std::cout << std::setprecision(7) << std::scientific << std::right << std::setw(5) << i << std::setw(22)
                    << betaEnergies[i] << std::endl;
        }
      }
    }

    cout << Core::Log::nl << Core::Log::endl;
    return true;
  }
};

} // namespace Readuct
} // namespace Scine

#endif // READUCT_SINGLEPOINTTASK_H_
