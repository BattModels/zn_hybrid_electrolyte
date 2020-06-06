//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef ZN_HYBRID_ELECTROLYTETESTAPP_H
#define ZN_HYBRID_ELECTROLYTETESTAPP_H

#include "MooseApp.h"

class zn_hybrid_electrolyteTestApp;

template <>
InputParameters validParams<zn_hybrid_electrolyteTestApp>();

class zn_hybrid_electrolyteTestApp : public MooseApp
{
public:
  zn_hybrid_electrolyteTestApp(InputParameters parameters);
  virtual ~zn_hybrid_electrolyteTestApp();

  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs = false);
};

#endif /* ZN_HYBRID_ELECTROLYTETESTAPP_H */
