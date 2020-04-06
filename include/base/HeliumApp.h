#pragma once

#include "MooseApp.h"

class Factory;

class HeliumApp : public MooseApp
{
public:
  HeliumApp(InputParameters parameters);

public:
  static InputParameters validParams();
  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);
};
