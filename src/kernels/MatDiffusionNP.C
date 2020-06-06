/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "MatDiffusionNP.h"

template <>
InputParameters
validParams<MatDiffusionNP>()
{
  InputParameters params = MatDiffusionBaseNP<Real>::validParams();
  params.addClassDescription(
      "Diffusion equation Kernel that takes an isotropic Diffusivity from a material property");
  return params;
}

MatDiffusionNP::MatDiffusionNP(const InputParameters & parameters) : MatDiffusionBaseNP<Real>(parameters)
{
}
