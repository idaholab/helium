#pragma once

#include "MooseApp.h"

class Factory;
class HeliumApp;

template <>
InputParameters validParams<HeliumApp>();

class HeliumApp : public MooseApp
{
public:
  HeliumApp(InputParameters parameters);

  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);

protected:
};
