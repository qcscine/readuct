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
#include <Utils/CalculatorBasics.h>
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

  bool run(SystemsMap& systems, Utils::UniversalSettings::ValueCollection taskSettings, bool testRunOnly = false,
           std::vector<std::function<void(const int&, const Utils::AtomCollection&, const Utils::Results&, const std::string&)>>
               observers = {}) const final {
    warningIfMultipleInputsGiven();
    warningIfMultipleOutputsGiven();

    // Read and delete special settings
    bool stopOnError = stopOnErrorExtraction(taskSettings);
    bool requireCharges = taskSettings.extract("require_charges", true);
    bool requireGradients = taskSettings.extract("require_gradients", false);
    bool requireStressTensor = taskSettings.extract("require_stress_tensor", false);
    bool requireBondOrders = taskSettings.extract("require_bond_orders", false);
    bool requireOrbitalEnergies = taskSettings.extract("orbital_energies", false);
    bool silentCalculator = taskSettings.extract("silent_stdout_calculator", true);
    int spinPropensityCheck = taskSettings.extract("spin_propensity_check", 0);
    if (!taskSettings.empty()) {
      std::string keyListing = "\n";
      auto keys = taskSettings.getKeys();
      for (const auto& key : keys) {
        keyListing += "'" + key + "'\n";
      }
      throw std::logic_error("Specified one or more task settings that are not available for this task:" + keyListing);
    }
    if (observers.size() > 0) {
      throw std::logic_error("SinglePointTask does not feature algorithm accepting observers, yet one was given");
    }

    if (testRunOnly) {
      return true; // leave out rest in case of task chaining
    }

    // Note: _input is guaranteed not to be empty by Task constructor
    auto calc = copyCalculator(systems, _input.front(), name());
    Utils::CalculationRoutines::setLog(*calc, true, true, !silentCalculator);
    // Check for available properties
    bool chargesAvailable = calc->possibleProperties().containsSubSet(Utils::Property::AtomicCharges);
    bool gradientsAvailable = calc->possibleProperties().containsSubSet(Utils::Property::Gradients);
    bool stressTensorAvailable = calc->possibleProperties().containsSubSet(Utils::Property::StressTensor);
    bool orbitalEnergiesAvailable = calc->possibleProperties().containsSubSet(Utils::Property::OrbitalEnergies);
    bool bondOrdersAvailable = calc->possibleProperties().containsSubSet(Utils::Property::BondOrderMatrix);
    if (requireCharges && !chargesAvailable) {
      throw std::logic_error("Charges required, but chosen calculator does not provide them.\n"
                             "If you do not need charges, set 'require_charges' to 'false' in the task settings");
    }
    if (requireGradients && !gradientsAvailable) {
      throw std::logic_error("Gradients required, but chosen calculator does not provide them.");
    }
    if (requireStressTensor && !stressTensorAvailable) {
      throw std::logic_error("Stress tensor required, but chosen calculator does not provide it.");
    }
    if (requireOrbitalEnergies && !orbitalEnergiesAvailable) {
      throw std::logic_error("Orbital energies required, but chosen calculator does not provide them.");
    }
    if (requireBondOrders && !bondOrdersAvailable) {
      throw std::logic_error("Bond orders required, but chosen calculator does not provide them.");
    }
    // Calculate energy, and possibly atomic charges and/or gradients
    Utils::PropertyList requiredProperties = Utils::Property::Energy;
    if (requireCharges) {
      requiredProperties.addProperty(Utils::Property::AtomicCharges);
    }
    if (requireGradients) {
      requiredProperties.addProperty(Utils::Property::Gradients);
    }
    if (requireStressTensor) {
      requiredProperties.addProperty(Utils::Property::StressTensor);
    }
    if (requireOrbitalEnergies) {
      requiredProperties.addProperty(Utils::Property::OrbitalEnergies);
    }
    if (requireBondOrders) {
      requiredProperties.addProperty(Utils::Property::BondOrderMatrix);
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
      return false;
    }

    if (spinPropensityCheck != 0) {
      if (!calc->settings().valueExists(Utils::SettingsNames::spinMultiplicity)) {
        _logger->warning << "Warning: " << calc->name()
                         << " does not allow multiplicity changes, skipping spin propensity check" << Core::Log::endl;
      }
      else {
        calc = Utils::CalculationRoutines::spinPropensity(*calc, *_logger, spinPropensityCheck);
      }
      // some calculators are lazy and don't copy properties, hence we recalculate if necessary
      if (!calc->getRequiredProperties().containsSubSet(requiredProperties)) {
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
          return false;
        }
      }
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
      cout << Core::Log::nl;
    }
    // Print stress tensor
    if (requireStressTensor) {
      Eigen::Matrix3d stressTensor = calc->results().get<Utils::Property::StressTensor>();
      cout << "  Stress tensor (hartree / bohr^3):\n\n";
      cout << [&stressTensor](std::ostream& os) { Utils::matrixPrettyPrint(os, stressTensor); };
      cout << Core::Log::nl;
    }
    // Print bond orders
    if (requireBondOrders) {
      cout << "  Atom#1 Atom#2 Bond Order\n\n";
      auto bos = calc->results().get<Utils::Property::BondOrderMatrix>();
      const auto& mat = bos.getMatrix();
      for (int i = 0; i < mat.outerSize(); i++) {
        for (typename Eigen::SparseMatrix<double>::InnerIterator it(mat, i); it; ++it) {
          if (it.value() > 0.3 && it.row() < it.col()) {
            cout.printf("  %6d %6d %+16.9f\n", it.row(), it.col(), it.value());
          }
        }
      }
      cout << Core::Log::nl << Core::Log::endl;
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
