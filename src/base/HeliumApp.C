#include "HeliumApp.h"
#include "HeliumRevision.h"
#include "MooseSyntax.h"
#include "AppFactory.h"

// Modules
#include "ModulesApp.h"

InputParameters
HeliumApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_output_syntax") = false;
  return params;
}

registerKnownLabel("HeliumApp");

HeliumApp::HeliumApp(InputParameters parameters) : MooseApp(parameters)
{
  HeliumApp::registerAll(_factory, _action_factory, _syntax);
}

// External entry point for dynamic application loading
extern "C" void
HeliumApp__registerApps()
{
  HeliumApp::registerApps();
}

void
HeliumApp::registerApps()
{
  registerApp(HeliumApp);
}

// External entry point for dynamic object registration
extern "C" void
HeliumApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  HeliumApp::registerAll(f, af, s);
}

void
HeliumApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  Registry::registerObjectsTo(f, {"HeliumApp"});
  Registry::registerActionsTo(af, {"HeliumApp"});

  ModulesApp::registerAllObjects<HeliumApp>(f, af, s);
}
