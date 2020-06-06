#include "zn_hybrid_electrolyteApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

template <>
InputParameters
validParams<zn_hybrid_electrolyteApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

zn_hybrid_electrolyteApp::zn_hybrid_electrolyteApp(InputParameters parameters) : MooseApp(parameters)
{
  zn_hybrid_electrolyteApp::registerAll(_factory, _action_factory, _syntax);
}

zn_hybrid_electrolyteApp::~zn_hybrid_electrolyteApp() {}

void
zn_hybrid_electrolyteApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"zn_hybrid_electrolyteApp"});
  Registry::registerActionsTo(af, {"zn_hybrid_electrolyteApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
zn_hybrid_electrolyteApp::registerApps()
{
  registerApp(zn_hybrid_electrolyteApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
zn_hybrid_electrolyteApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  zn_hybrid_electrolyteApp::registerAll(f, af, s);
}
extern "C" void
zn_hybrid_electrolyteApp__registerApps()
{
  zn_hybrid_electrolyteApp::registerApps();
}
