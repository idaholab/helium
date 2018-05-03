#ifndef HELIUMFLUIDPROPERTIESTEST_H
#define HELIUMFLUIDPROPERTIESTEST_H

#include "HeliumFluidProperties.h"
#include "MooseObjectUnitTest.h"

class HeliumFluidPropertiesTest : public MooseObjectUnitTest
{
public:
  HeliumFluidPropertiesTest() : MooseObjectUnitTest("HeliumApp") { buildObjects(); }

protected:
  void buildObjects()
  {
    InputParameters uo_pars = _factory.getValidParams("HeliumFluidProperties");
    _fe_problem->addUserObject("HeliumFluidProperties", "fp", uo_pars);
    _fp = &_fe_problem->getUserObject<HeliumFluidProperties>("fp");
  }

  const HeliumFluidProperties * _fp;
};

#endif
