//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef ZN_HYBRID_ELECTROLYTEAPP_H
#define ZN_HYBRID_ELECTROLYTEAPP_H

#include "MooseApp.h"

class zn_hybrid_electrolyteApp;

template <>
InputParameters validParams<zn_hybrid_electrolyteApp>();

class zn_hybrid_electrolyteApp : public MooseApp
{
public:
  zn_hybrid_electrolyteApp(InputParameters parameters);
  virtual ~zn_hybrid_electrolyteApp();

  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);
};

#endif /* ZN_HYBRID_ELECTROLYTEAPP_H */
