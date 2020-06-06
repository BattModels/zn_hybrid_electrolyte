/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef OVERPOTENTIAL_H
#define OVERPOTENTIAL_H

#include "Kernel.h"
#include "DerivativeMaterialInterface.h"

class Overpotential;

template <>
InputParameters validParams<Overpotential>();

class Overpotential: public DerivativeMaterialInterface<Kernel>
{
public:
  Overpotential(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  unsigned int _odp_var;
  const VariableValue & _odp;
  const VariableGradient & _grad_odp;
  unsigned int _ep_var;
  const VariableValue & _ep;
  unsigned int _miu_var;
  const VariableValue & _miu;
  /// Mobility
  const MaterialProperty<Real> & _K;
  const MaterialProperty<Real> & _T;
  const MaterialProperty<Real> & _dTdpot;
  const MaterialProperty<Real> & _dTdeta;
  const MaterialProperty<Real> & _dTdmiu;
  /// Interfacial parameter
};

#endif // OVERPOTENTIAL_H
