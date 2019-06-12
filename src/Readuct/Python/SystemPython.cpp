/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Core/Interfaces/Calculator.h>
#include <Core/ModuleManager.h>
#include <Utils/Geometry.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/IO/Yaml.h>
#include <pybind11/pybind11.h>
#include <iostream>
// #include <pybind11/stl.h>

using namespace Scine;

std::shared_ptr<Core::Calculator> loadSystem(std::string path, std::string method, pybind11::kwargs kwargs) {
  // Load module manager
  auto& manager = Core::ModuleManager::getInstance();
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
  std::shared_ptr<Core::Calculator> calc;
  try {
    calc = manager.get<Core::Calculator>(method, program);
  }
  catch (...) {
    if (program.empty()) {
      std::cout << "No SCINE module providing '" << method << "' is currently loaded.\n";
      std::cout << "Please add the module to the SCINE_MODULE_PATH in order for it to be accessible.\n";
      throw std::runtime_error("Failed to load method/program.");
    }
    else {
      std::cout << "No SCINE module named '" << program << "' providing '" << method << "' is currently loaded.\n";
      std::cout << "Please add the module to the SCINE_MODULE_PATH in order for it to be accessible.\n";
      throw std::runtime_error("Failed to load method/program.");
    }
  }
  // Apply settings to Calculator
  nodeToSettings(calc->settings(), settingsnode);
  // Set initial structure
  calc->setStructure(readResults.first);
  return calc;
}

void init_system(pybind11::module& m) {
  m.def("load_system", &loadSystem, pybind11::arg("path"), pybind11::arg("method"));
}
