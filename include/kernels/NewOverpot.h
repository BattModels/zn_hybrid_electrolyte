/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef NewOverpot_H
#define NewOverpot_H

#include "Kernel.h"
#include "DerivativeMaterialInterface.h"

class NewOverpot;

template <>
InputParameters validParams<NewOverpot>();

class NewOverpot: public DerivativeMaterialInterface<Kernel>
{
public:
  NewOverpot(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  unsigned int _odp_var;
  const VariableValue & _odp;
  const VariableGradient & _grad_odp;
  unsigned int _ep_var;
  const VariableValue & _ep;
  const VariableGradient & _grad_ep;
  unsigned int _miu_var;
  const VariableValue & _miu;
  /// Mobility
  const MaterialProperty<Real> & _K;
  const MaterialProperty<Real> & _E;
  const MaterialProperty<Real> & _T;
  const MaterialProperty<Real> & _dEdeta;
  const MaterialProperty<Real> & _dTdeta;
  const MaterialProperty<Real> & _dTdmiu;
  /// Interfacial parameter
};

#endif // NewOverpot_H
