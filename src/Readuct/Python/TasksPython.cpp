/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Tasks/AfirOptimizationTask.h"
#include "Tasks/BSplineInterpolationTask.h"
#include "Tasks/BondOrderTask.h"
#include "Tasks/GeometryOptimizationTask.h"
#include "Tasks/HessianTask.h"
#include "Tasks/IrcTask.h"
#include "Tasks/NtOptimization2Task.h"
#include "Tasks/NtOptimizationTask.h"
#include "Tasks/SinglePointTask.h"
#include "Tasks/Task.h"
#include "Tasks/TsOptimizationTask.h"
#include <Core/Interfaces/Calculator.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>

using namespace Scine;

namespace {

using SystemsMap = std::map<std::string, std::shared_ptr<Core::Calculator>>;

template<class TaskType>
std::pair<SystemsMap, bool> run(SystemsMap systems, std::vector<std::string> inputNames, bool testRunOnly,
                                const pybind11::kwargs& kwargs) {
  Utils::UniversalSettings::ValueCollection settingValues;

  std::vector<std::string> output;
  std::shared_ptr<Core::Log> logger = nullptr;
  for (auto item : kwargs) {
    std::string key = item.first.cast<std::string>();
    if (key == "output") {
      pybind11::list outputNames = item.second.cast<pybind11::list>();
      for (auto item : outputNames) {
        std::string name = item.cast<std::string>();
        output.push_back(name);
      }
    }
    else if (key == "logger") {
      logger = item.second.cast<std::shared_ptr<Core::Log>>();
    }
    else if (pybind11::isinstance<pybind11::bool_>(item.second)) {
      bool value = item.second.cast<bool>();
      settingValues.addBool(key, value);
    }
    else if (pybind11::isinstance<pybind11::str>(item.second)) {
      std::string value = item.second.cast<std::string>();
      settingValues.addString(key, value);
    }
    else if (pybind11::isinstance<pybind11::int_>(item.second)) {
      int value = item.second.cast<int>();
      settingValues.addInt(key, value);
    }
    else if (pybind11::isinstance<pybind11::float_>(item.second)) {
      double value = item.second.cast<double>();
      settingValues.addDouble(key, value);
    }
    else if (pybind11::isinstance<pybind11::list>(item.second)) {
      pybind11::list list = item.second.cast<pybind11::list>();
      if (!list.empty()) {
        // Figure out type from first item
        const auto& firstItem = *list.begin();
        if (pybind11::isinstance<pybind11::int_>(firstItem)) {
          std::vector<int> ints;
          for (const auto& item : list) {
            ints.push_back(item.cast<int>());
          }
          settingValues.addIntList(key, ints);
        }
        else if (pybind11::isinstance<pybind11::str>(firstItem)) {
          std::vector<std::string> strings;
          for (const auto& item : list) {
            strings.push_back(item.cast<std::string>());
          }
          settingValues.addStringList(key, strings);
        }
        else {
          throw std::runtime_error("CollectionLists are unimplemented here!");
        }
      }
    }
    else {
      throw std::runtime_error(key + " could not be converted from Python into a C++ variable, check its type!");
    }
  }

  // Run the task
  TaskType task(std::move(inputNames), std::move(output), std::move(logger));
  const bool success = task.run(systems, settingValues, testRunOnly);

  return {systems, success};
}

} // namespace

void init_tasks(pybind11::module& m) {
  m.def("run_single_point_task", &run<Readuct::SinglePointTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"), pybind11::arg("test_run_only") = false);
  m.def("run_sp_task", &run<Readuct::SinglePointTask>, pybind11::arg("systems"), pybind11::arg("names_of_used_systems"),
        pybind11::arg("test_run_only") = false);

  m.def("run_optimization_task", &run<Readuct::GeometryOptimizationTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"), pybind11::arg("test_run_only") = false);
  m.def("run_opt_task", &run<Readuct::GeometryOptimizationTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"), pybind11::arg("test_run_only") = false);

  m.def("run_bond_order_task", &run<Readuct::BondOrderTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"), pybind11::arg("test_run_only") = false);

  m.def("run_hessian_task", &run<Readuct::HessianTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"), pybind11::arg("test_run_only") = false);
  m.def("run_freq_task", &run<Readuct::HessianTask>, pybind11::arg("systems"), pybind11::arg("names_of_used_systems"),
        pybind11::arg("test_run_only") = false);

  m.def("run_transition_state_optimization_task", &run<Readuct::TsOptimizationTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"), pybind11::arg("test_run_only") = false);
  m.def("run_tsopt_task", &run<Readuct::TsOptimizationTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"), pybind11::arg("test_run_only") = false);

  m.def("run_irc_task", &run<Readuct::IrcTask>, pybind11::arg("systems"), pybind11::arg("names_of_used_systems"),
        pybind11::arg("test_run_only") = false);

  m.def("run_nt_task", &run<Readuct::NtOptimizationTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"), pybind11::arg("test_run_only") = false);

  m.def("run_nt2_task", &run<Readuct::NtOptimization2Task>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"), pybind11::arg("test_run_only") = false);

  m.def("run_afir_task", &run<Readuct::AfirOptimizationTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"), pybind11::arg("test_run_only") = false);

  m.def("run_bspline_task", &run<Readuct::BSplineInterpolationTask>, pybind11::arg("systems"),
        pybind11::arg("names_of_used_systems"), pybind11::arg("test_run_only") = false);
}
