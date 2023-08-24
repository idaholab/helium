//* This file is part of helium
//* https://github.com/idaholab/helium
//*
//* All rights reserved, see NOTICE.txt for full restrictions
//* https://github.com/idaholab/helium/blob/master/NOTICE.txt
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

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
