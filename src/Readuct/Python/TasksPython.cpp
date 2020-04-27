/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Tasks/AfirOptimizationTask.h"
#include "Tasks/BSplineInterpolationTask.h"
#include "Tasks/BondOrderTask.h"
#include "Tasks/GeometryOptimizationTask.h"
#include "Tasks/HessianTask.h"
#include "Tasks/IrcTask.h"
#include "Tasks/SinglePointTask.h"
#include "Tasks/Task.h"
#include "Tasks/TsOptimizationTask.h"
#include <Core/Interfaces/Calculator.h>
#include <Utils/IO/Logger.h>
#include <Utils/IO/Yaml.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <pybind11/pybind11.h>
#include <iostream>

using namespace Scine;

template<class TaskName>
std::tuple<pybind11::dict, bool> run(pybind11::dict inputSystems, pybind11::list inputNames, pybind11::kwargs kwargs) {
  // Start logger
  std::string loggingVerbosity = "info";
  if (kwargs.contains(Utils::SettingsNames::loggerVerbosity)) {
    loggingVerbosity = kwargs[Utils::SettingsNames::loggerVerbosity].cast<std::string>();
  }
  Utils::Log::startConsoleLogging(loggingVerbosity);

  // Read all systems
  std::map<std::string, std::shared_ptr<Core::Calculator>> systems;
  for (auto item : inputSystems) {
    std::string key = item.first.cast<std::string>();
    std::shared_ptr<Core::Calculator> value = item.second.cast<std::shared_ptr<Core::Calculator>>();
    systems[key] = value;
  }
  // Read all input system names
  std::vector<std::string> input;
  for (auto item : inputNames) {
    std::string name = item.cast<std::string>();
    input.push_back(name);
  }
  // Read all settings
  YAML::Node settingsnode;
  std::vector<std::string> output;
  for (auto item : kwargs) {
    std::string key = item.first.cast<std::string>();
    if (!key.compare("output")) {
      pybind11::list outputNames = item.second.cast<pybind11::list>();
      for (auto item : outputNames) {
        std::string name = item.cast<std::string>();
        output.push_back(name);
      }
    }
    else if (pybind11::isinstance<pybind11::bool_>(item.second)) {
      bool value = item.second.cast<bool>();
      settingsnode[key] = value;
    }
    else if (pybind11::isinstance<pybind11::str>(item.second)) {
      std::string value = item.second.cast<std::string>();
      settingsnode[key] = value;
    }
    else if (pybind11::isinstance<pybind11::int_>(item.second)) {
      int value = item.second.cast<int>();
      settingsnode[key] = value;
    }
    else if (pybind11::isinstance<pybind11::float_>(item.second)) {
      double value = item.second.cast<double>();
      settingsnode[key] = value;
    }
    else if (pybind11::isinstance<pybind11::list>(item.second)) {
      pybind11::list list = item.second.cast<pybind11::list>();
      std::vector<std::string> strings;
      for (auto item : list) {
        std::string str = pybind11::str(item).cast<std::string>();
        strings.push_back(str);
      }
      settingsnode[key] = strings;
    }
    else {
      throw std::runtime_error(key + " could not be converted from Python into C++ a variable, check its type!");
    }
  }
  // Run
  auto task = std::make_unique<TaskName>(input, output);
  auto success = task->run(systems, settingsnode);
  // Update systems
  for (auto item : systems) {
    inputSystems[pybind11::cast(item.first)] = item.second;
  }

  return {inputSystems, success};
}

void init_tasks(pybind11::module& m) {
  m.def("run_single_point_task", &run<Readuct::SinglePointTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"));
  m.def("run_sp_task", &run<Readuct::SinglePointTask>, pybind11::arg("systems"), pybind11::arg("names_of_used_systems"));
  m.def("run_optimization_task", &run<Readuct::GeometryOptimizationTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"));
  m.def("run_bond_order_task", &run<Readuct::BondOrderTask>, pybind11::arg("systems"), pybind11::arg("names_of_used_systems"));
  m.def("run_opt_task", &run<Readuct::GeometryOptimizationTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"));
  m.def("run_hessian_task", &run<Readuct::HessianTask>, pybind11::arg("systems"), pybind11::arg("names_of_used_systems"));
  m.def("run_freq_task", &run<Readuct::HessianTask>, pybind11::arg("systems"), pybind11::arg("names_of_used_systems"));
  m.def("run_transition_state_optimization_task", &run<Readuct::TsOptimizationTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"));
  m.def("run_tsopt_task", &run<Readuct::TsOptimizationTask>, pybind11::arg("systems"), pybind11::arg("names_of_used_systems"));
  m.def("run_irc_task", &run<Readuct::IrcTask>, pybind11::arg("systems"), pybind11::arg("names_of_used_systems"));
  m.def("run_afir_task", &run<Readuct::AfirOptimizationTask>, pybind11::arg("systems"), pybind11::arg("names_of_used_systems"));
  m.def("run_bspline_task", &run<Readuct::BSplineInterpolationTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"));
}
