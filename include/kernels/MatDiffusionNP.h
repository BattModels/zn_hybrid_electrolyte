/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef MATDIFFUSIONNP_H
#define MATDIFFUSIONNP_H

#include "MatDiffusionBaseNP.h"

/**
 * Isotropic diffusion kernel that takes a diffusion coefficient of type
 * Real. All logic is implemnted in the MatDiffusionBase class
 * template.
 */
class MatDiffusionNP : public MatDiffusionBaseNP<Real>
{
public:
  MatDiffusionNP(const InputParameters & parameters);
};

template <>
InputParameters validParams<MatDiffusionNP>();

#endif // MATDIFFUSIONNP_H
