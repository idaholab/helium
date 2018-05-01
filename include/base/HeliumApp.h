#ifndef HELIUMAPP_H
#define HELIUMAPP_H

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
  static void registerObjects(Factory & factory);

protected:
};

#endif /* HELIUMAPP_H */
