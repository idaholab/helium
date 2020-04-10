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
