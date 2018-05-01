#ifndef HELIUMFLUIDPROPERTIESTEST_H
#define HELIUMFLUIDPROPERTIESTEST_H

#include "gtest/gtest.h"

#include "FEProblem.h"
#include "AppFactory.h"
#include "GeneratedMesh.h"
#include "HeliumFluidProperties.h"

class HeliumFluidPropertiesTest : public ::testing::Test
{
public:
  HeliumFluidPropertiesTest()
    : _app(AppFactory::createAppShared("HeliumApp", 0, nullptr)), _factory(_app->getFactory())
  {
  }

protected:
  void SetUp()
  {
    Moose::registerObjects(_factory);
    registerObjects();
    buildObjects();
  }

  void TearDown()
  {
    delete _fe_problem;
    delete _mesh;
  }

  void registerObjects() { }

  void buildObjects()
  {
    InputParameters mesh_params = _factory.getValidParams("GeneratedMesh");
    mesh_params.set<MooseEnum>("dim") = "3";
    mesh_params.set<std::string>("name") = "mesh";
    mesh_params.set<std::string>("_object_name") = "name1";
    _mesh = new GeneratedMesh(mesh_params);

    InputParameters problem_params = _factory.getValidParams("FEProblem");
    problem_params.set<MooseMesh *>("mesh") = _mesh;
    problem_params.set<std::string>("name") = "problem";
    problem_params.set<std::string>("_object_name") = "name2";
    _fe_problem = new FEProblem(problem_params);

    InputParameters uo_pars = _factory.getValidParams("HeliumFluidProperties");
    _fe_problem->addUserObject("HeliumFluidProperties", "fp", uo_pars);
    _fp = &_fe_problem->getUserObject<HeliumFluidProperties>("fp");
  }

  std::shared_ptr<MooseApp> _app;
  Factory & _factory;
  MooseMesh * _mesh;
  FEProblem * _fe_problem;
  const HeliumFluidProperties * _fp;
};
#endif
