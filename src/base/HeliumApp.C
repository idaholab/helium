#include "HeliumApp.h"
#include "HeliumRevision.h"
#include "MooseSyntax.h"
#include "AppFactory.h"

// Modules
#include "FluidPropertiesApp.h"

// Fluid properties
#include "HeliumFluidProperties.h"

template <>
InputParameters
validParams<HeliumApp>()
{
  InputParameters params = validParams<MooseApp>();
  params.set<bool>("use_legacy_output_syntax") = false;
  return params;
}

HeliumApp::HeliumApp(InputParameters parameters) : MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  FluidPropertiesApp::registerObjects(_factory);
  HeliumApp::registerObjects(_factory);

  FluidPropertiesApp::associateSyntax(_syntax, _action_factory);
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
HeliumApp__registerObjects(Factory & factory)
{
  HeliumApp::registerObjects(factory);
}

void
HeliumApp::registerObjects(Factory & factory)
{
  Registry::registerObjectsTo(factory, {"HeliumApp"});
}
