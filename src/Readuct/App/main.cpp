/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

/* Readuct */
#include "Tasks/TaskFactory.h"
#include "io.h"
/* Scine */
#include <Core/Interfaces/Calculator.h>
#include <Core/Log.h>
#include <Utils/IO/Yaml.h>

/* Boost program arguments */
#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

/* C++ std */
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace Scine;
using namespace Scine::Readuct;

int main(int argc, char* argv[]) {
  namespace po = boost::program_options;

  // Arguments
  po::options_description optionsDescription("Recognized options");
  optionsDescription.add_options()("help,h", "Produce this help message")("config,c", po::value<std::string>(),
                                                                          "YAML input file to read");

  po::variables_map optionsMap;
  po::positional_options_description positionalDescription;
  positionalDescription.add("config", 1);
  po::store(po::command_line_parser(argc, argv)
                .options(optionsDescription)
                .positional(positionalDescription)
                .style(po::command_line_style::unix_style | po::command_line_style::allow_long_disguise)
                .run(),
            optionsMap);
  po::notify(optionsMap);

  // Initialize Logger
  auto logger = std::make_shared<Core::Log>();
  auto cout = logger->output;

  // Handle help
  if (optionsMap.count("help") > 0 || optionsMap.count("config") == 0) {
    cout << optionsDescription << Core::Log::endl;
    return 0;
  }

  // Load input file
  const std::string filename = optionsMap["config"].as<std::string>();
  if (!boost::filesystem::exists(filename)) {
    cout << "Specified config file does not exist!\n";
    return 1;
  }

  // Header
  cout << R"(#=============================================================================#)" << Core::Log::endl;
  cout << R"(|   ______   _______     ___       _______   __    __    ______  __________   |)" << Core::Log::endl;
  cout << R"(|  |   _  \ |   ____|   /   \     |       \ |  |  |  |  /      ||          |  |)" << Core::Log::endl;
  cout << R"(|  |  |_)  ||  |__     /  ^  \    |  .--.  ||  |  |  | |  ,----'`---|  |---`  |)" << Core::Log::endl;
  cout << R"(|  |      / |   __|   /  /_\  \   |  |  |  ||  |  |  | |  |         |  |      |)" << Core::Log::endl;
  cout << R"(|  |  |\  \ |  |____ /  _____  \  |  '--'  ||  `--'  | |  `----.    |  |      |)" << Core::Log::endl;
  cout << R"(|  | _| `._||_______/__/     \__\ |_______/  \______/   \______|    |__|      |)" << Core::Log::endl;
  cout << R"(|                                                                             |)" << Core::Log::endl;
  cout << R"(#=============================================================================#)" << Core::Log::endl;
  cout << Core::Log::endl << Core::Log::endl;

  std::map<std::string, std::shared_ptr<Core::Calculator>> systems;
  std::vector<std::shared_ptr<Task>> tasks;
  std::vector<Utils::UniversalSettings::ValueCollection> tasksettings;
  try {
    auto data = loadYamlFile(filename);
    auto systemsWithMetaInfo = std::get<0>(data);
    for (const auto& x : systemsWithMetaInfo) {
      systems[x.first] = std::get<2>(x.second);
    }
    tasks = std::get<1>(data);
    tasksettings = std::get<2>(data);
  }
  catch (std::logic_error& e) {
    cout << e.what() << Core::Log::endl;
    return 1;
  }

  // Run tasks
  for (size_t i = 0; i < tasks.size(); i++) {
    auto name = tasks[i]->name();
    cout.printf("  #============#    #==%s==#\n", std::string(name.size(), '=').c_str());
    cout.printf("  |  Task %3d  |    |  %s  |\n", static_cast<int>(i + 1), name.c_str());
    cout.printf("  #============#    #==%s==#\n", std::string(name.size(), '=').c_str());
    cout << Core::Log::endl << Core::Log::endl;
    const auto& inpt = tasks[i]->input();
    const auto& outpt = tasks[i]->output();

    cout << "  Input System(s):  " + tasks[i]->input()[0] << Core::Log::endl;
    for (unsigned int j = 1; j < inpt.size(); j++) {
      cout << "                    " + inpt[j] << Core::Log::endl;
    }
    cout << Core::Log::endl;
    if (!outpt.empty()) {
      cout << "  Output System(s): " + outpt[0] << Core::Log::endl;
      for (unsigned int j = 1; j < outpt.size(); j++) {
        cout << "                    " + outpt[j] << Core::Log::endl;
      }
      cout << Core::Log::endl;
    }
    tasks[i]->run(systems, tasksettings[i]);
  }

  cout << "  Completed all tasks, exiting!" << Core::Log::endl;
  return 0;
}
