/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef READUCT_BONDORDERTASK_H_
#define READUCT_BONDORDERTASK_H_

/* Readuct */
#include "Tasks/Task.h"
/* Scine */
#include <Core/Interfaces/Calculator.h>
#include <Utils/CalculatorBasics/Results.h>
/* External */
#include "boost/exception/diagnostic_information.hpp"
#include <boost/filesystem.hpp>
/* std c++ */
#include <cstdio>
#include <fstream>
#include <iostream>

namespace Scine {
namespace Readuct {

class BondOrderTask : public Task {
 public:
  /**
   * @brief Construct a new BondOrderTask.
   * @param input  The input system names for the task.
   * @param output The output system names for the task.
   * @param logger The logger to/through which all text output will be handled.
   */
  BondOrderTask(std::vector<std::string> input, std::vector<std::string> output, std::shared_ptr<Core::Log> logger = nullptr)
    : Task(std::move(input), std::move(output), std::move(logger)) {
  }

  std::string name() const override {
    return "Bond Order Calculation";
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

    // Calculate bond orders and energy if not present in the results yet
    if (!calc->results().has<Utils::Property::BondOrderMatrix>()) {
      calc->setRequiredProperties(Utils::Property::BondOrderMatrix | Utils::Property::Energy);
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

    auto energy = calc->results().get<Utils::Property::Energy>();
    auto bos = calc->results().get<Utils::Property::BondOrderMatrix>();

    // Print
    auto cout = _logger->output;
    cout.printf("  The (electronic) energy is: %+16.9f hartree\n\n\n", energy);
    cout << "  Showing all bond orders >0.3 a.u.:\n\n";
    cout << "  Atom#1 Atom#2 Bond Order\n\n";
    const auto& mat = bos.getMatrix();
    for (int i = 0; i < mat.outerSize(); i++) {
      for (typename Eigen::SparseMatrix<double>::InnerIterator it(mat, i); it; ++it) {
        if (it.value() > 0.3 && it.row() < it.col()) {
          cout.printf("  %6d %6d %+16.9f\n", it.row(), it.col(), it.value());
        }
      }
    }
    cout << Core::Log::nl << Core::Log::endl;

    // Store result
    if (!_output.empty()) {
      systems[_output.front()] = std::move(calc);
    }
    else {
      systems[_input.front()] = std::move(calc);
    }

    return true;
  }
};

} // namespace Readuct
} // namespace Scine

#endif // READUCT_BONDORDERTASK_H_
