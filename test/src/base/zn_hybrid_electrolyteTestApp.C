//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "zn_hybrid_electrolyteTestApp.h"
#include "zn_hybrid_electrolyteApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

template <>
InputParameters
validParams<zn_hybrid_electrolyteTestApp>()
{
  InputParameters params = validParams<zn_hybrid_electrolyteApp>();
  return params;
}

zn_hybrid_electrolyteTestApp::zn_hybrid_electrolyteTestApp(InputParameters parameters) : MooseApp(parameters)
{
  zn_hybrid_electrolyteTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

zn_hybrid_electrolyteTestApp::~zn_hybrid_electrolyteTestApp() {}

void
zn_hybrid_electrolyteTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  zn_hybrid_electrolyteApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"zn_hybrid_electrolyteTestApp"});
    Registry::registerActionsTo(af, {"zn_hybrid_electrolyteTestApp"});
  }
}

void
zn_hybrid_electrolyteTestApp::registerApps()
{
  registerApp(zn_hybrid_electrolyteApp);
  registerApp(zn_hybrid_electrolyteTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
zn_hybrid_electrolyteTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  zn_hybrid_electrolyteTestApp::registerAll(f, af, s);
}
extern "C" void
zn_hybrid_electrolyteTestApp__registerApps()
{
  zn_hybrid_electrolyteTestApp::registerApps();
}
