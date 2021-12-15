/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef READUCT_READUCTMODULE_H
#define READUCT_READUCTMODULE_H

#include <Core/Module.h>
#include <boost/dll/alias.hpp>
#include <boost/hana/define_struct.hpp>
#include <memory>

namespace Scine {
namespace Readuct {

class ReaductModule : public Scine::Core::Module {
 public:
  BOOST_HANA_DEFINE_STRUCT(ReaductModule, (std::vector<std::string>, transition_state_optimizer),
                           (std::vector<std::string>, reaction_path_optimizer));

  ReaductModule() noexcept;

  std::string name() const noexcept final;

  boost::any get(const std::string& interface, const std::string& model) const final;

  bool has(const std::string& interface, const std::string& model) const noexcept final;

  std::vector<std::string> announceInterfaces() const noexcept final;

  std::vector<std::string> announceModels(const std::string& concept) const noexcept final;

  static std::shared_ptr<Module> make();
};

} // namespace Readuct
} // namespace Scine

// At global namespace, define the entry point for the module.
BOOST_DLL_ALIAS(Scine::Readuct::ReaductModule::make, moduleFactory)

#endif // READUCT_READUCTMODULE_H
