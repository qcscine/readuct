/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ReaductModule.h"
#include <Core/DerivedModule.h>
#include <Core/Interfaces/Calculator.h>
#include <Core/Interfaces/ReactionPathOptimizer.h>
#include <Core/Interfaces/TransitionStateOptimizer.h>
#include <Readuct/ReaductReactionPathOptimizer.h>
#include <Readuct/ReaductTransitionStateOptimizer.h>

namespace Scine {
namespace Readuct {

ReaductModule::ReaductModule() noexcept
  : transition_state_optimizer{ReaductTransitionStateOptimizer::model}, reaction_path_optimizer{ReaductReactionPathOptimizer::model} {
  Core::DerivedModule::checkInvariants(*this);
}

std::string ReaductModule::name() const noexcept {
  return "ReaDuct";
}

boost::any ReaductModule::get(const std::string& interface, const std::string& model) const {
  using tsOptimizerPtr = std::shared_ptr<Core::TransitionStateOptimizer>;
  using reactionPathOptimizerPtr = std::shared_ptr<Core::ReactionPathOptimizer>;
  if (interface == Core::TransitionStateOptimizer::interface) {
    if (model == ReaductTransitionStateOptimizer::model) {
      return static_cast<tsOptimizerPtr>(std::make_shared<ReaductTransitionStateOptimizer>());
    }
  }
  else if (interface == Core::ReactionPathOptimizer::interface) {
    if (model == ReaductReactionPathOptimizer::model) {
      return static_cast<reactionPathOptimizerPtr>(std::make_shared<ReaductReactionPathOptimizer>());
    }
  }
  // Throw an exception if we cannot match a model
  throw Core::ClassNotImplementedError();
}

bool ReaductModule::has(const std::string& interface, const std::string& model) const noexcept {
  return Core::DerivedModule::has(interface, model, *this);
}

std::vector<std::string> ReaductModule::announceInterfaces() const noexcept {
  return Core::DerivedModule::announceInterfaces(*this);
}

std::vector<std::string> ReaductModule::announceModels(const std::string& concept) const noexcept {
  return Core::DerivedModule::announceModels(concept, *this);
}

std::shared_ptr<Core::Module> ReaductModule::make() {
  return std::make_shared<ReaductModule>();
}

} // namespace Readuct
} // namespace Scine
