//* This file is part of helium
//* https://github.com/idaholab/helium
//*
//* All rights reserved, see NOTICE.txt for full restrictions
//* https://github.com/idaholab/helium/blob/master/NOTICE.txt
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeliumApp.h"
// Moose Includes
#include "MooseInit.h"
#include "Moose.h"
#include "MooseApp.h"
#include "AppFactory.h"

// Create a performance log
PerfLog Moose::perf_log("Helium");

// Begin the main program.
int
main(int argc, char * argv[])
{
  // Initialize MPI, solvers and MOOSE
  MooseInit init(argc, argv);

  // Register this application's MooseApp and any it depends on
  HeliumApp::registerApps();

  std::shared_ptr<MooseApp> app = AppFactory::createAppShared("HeliumApp", argc, argv);

  // Execute the application
  app->run();

  return 0;
}
