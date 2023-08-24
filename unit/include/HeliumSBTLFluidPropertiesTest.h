//* This file is part of helium
//* https://github.com/idaholab/helium
//*
//* All rights reserved, see NOTICE.txt for full restrictions
//* https://github.com/idaholab/helium/blob/master/NOTICE.txt
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "HeliumSBTLFluidProperties.h"
#include "MooseObjectUnitTest.h"

class HeliumSBTLFluidPropertiesTest : public MooseObjectUnitTest
{
public:
  HeliumSBTLFluidPropertiesTest() : MooseObjectUnitTest("HeliumApp") { buildObjects(); }

protected:
  void buildObjects()
  {
    InputParameters uo_pars = _factory.getValidParams("HeliumSBTLFluidProperties");
    _fe_problem->addUserObject("HeliumSBTLFluidProperties", "fp", uo_pars);
    _fp = &_fe_problem->getUserObject<HeliumSBTLFluidProperties>("fp");
  }

  const HeliumSBTLFluidProperties * _fp;
};
