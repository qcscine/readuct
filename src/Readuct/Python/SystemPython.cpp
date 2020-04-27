/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Core/Interfaces/Calculator.h>
#include <Core/ModuleManager.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/Geometry.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/Yaml.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>

using namespace Scine;

std::shared_ptr<Core::Calculator> getCalculator(std::string method_family, std::string program) {
  // Load module manager
  auto& manager = Core::ModuleManager::getInstance();

  for (auto& x : program)
    x = std::tolower(x);
  program[0] = std::toupper(program[0]);
  for (auto& x : method_family)
    x = std::toupper(x);

  // Generate Calculator
  std::shared_ptr<Core::Calculator> calc;
  try {
    calc = manager.get<Core::Calculator>(Core::Calculator::supports(method_family), program);
  }
  catch (...) {
    if (program.empty()) {
      std::cout << "No SCINE module providing '" << method_family << "' is currently loaded.\n";
      std::cout << "Please add the module to the SCINE_MODULE_PATH in order for it to be accessible.\n";
      throw std::runtime_error("Failed to load method/program.");
    }
    else {
      std::cout << "No SCINE module named '" << program << "' providing '" << method_family << "' is currently loaded.\n";
      std::cout << "Please add the module to the SCINE_MODULE_PATH in order for it to be accessible.\n";
      throw std::runtime_error("Failed to load method/program.");
    }
  }
  // Return Calculator
  return calc;
}

std::shared_ptr<Core::Calculator> loadSystem(std::string path, std::string method_family, pybind11::kwargs kwargs) {
  // Read all optional arguments and turn them into a node
  std::string program = "";
  YAML::Node settingsnode;
  for (auto item : kwargs) {
    std::string key = item.first.cast<std::string>();
    std::string value = pybind11::str(item.second).cast<std::string>();
    if (!key.compare("program")) {
      program = value;
    }
    else {
      settingsnode[key] = value;
    }
  }

  // Load molecule
  auto readResults = Utils::ChemicalFileHandler::read(path);
  // Generate Calculator
  auto calc = getCalculator(method_family, program);
  // Apply settings to Calculator
  nodeToSettings(calc->settings(), settingsnode);
  // Set initial structure
  calc->setStructure(readResults.first);
  return calc;
}

std::vector<std::string> getAvailableSettings(std::string method_family, std::string program) {
  // Get the correct Calculator
  auto calc = getCalculator(method_family, program);

  // Get the names of the available settings from the calculator's settings
  std::vector<std::string> availableSettings;
  Utils::UniversalSettings::DescriptorCollection settings = calc->settings().getDescriptorCollection();
  for (const auto& s : settings) {
    availableSettings.push_back(s.first);
  }
  return availableSettings;
}

Utils::PropertyList getPossiblePropertiesByStrings(std::string method_family, std::string program) {
  // Get the correct Calculator
  auto calc = getCalculator(method_family, program);
  return calc->possibleProperties();
}

void init_system(pybind11::module& m) {
  m.def("load_system", &loadSystem, pybind11::arg("path"), pybind11::arg("method_family"));
  m.def("get_available_settings", &getAvailableSettings, pybind11::arg("method_family"), pybind11::arg("program"));
  m.def("get_possible_properties", &getPossiblePropertiesByStrings, pybind11::arg("method_family"), pybind11::arg("program"));
}
