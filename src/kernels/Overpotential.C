/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "Overpotential.h"

template <>
InputParameters
validParams<Overpotential>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription(
      "Add in Overpotential");
  params.addParam<MaterialPropertyName>("K_name","Gradient coefficient, use kappa in the input file");
  params.addRequiredCoupledVar(
      "odp", "orderparameter");
  params.addRequiredCoupledVar(
      "ep", "coupled variable 2,the electric potential");
  params.addRequiredCoupledVar(
        "miu", "coupled variable 3,the chemical potential");
  params.addRequiredParam<MaterialPropertyName>(
        "t_name", "Corresponding to overpotential0 in the main input");
  return params;
}

Overpotential::Overpotential(const InputParameters & parameters)
: DerivativeMaterialInterface<Kernel>(parameters),
 _odp_var(coupled("odp")),
 _odp(coupledValue("odp")),
 _grad_odp(coupledGradient("odp")),
  _ep_var(coupled("ep")),
  _ep(coupledValue("ep")),
  _miu_var(coupled("miu")),
  _miu(coupledValue("miu")),
  _K(getMaterialProperty<Real>("K_name")),
  _T(getMaterialProperty<Real>("t_name")),
  _dTdpot(getMaterialPropertyDerivative<Real>("t_name", getVar("ep", 0)->name())),
  _dTdeta(getMaterialPropertyDerivative<Real>("t_name", getVar("odp", 0)->name())),
  _dTdmiu(getMaterialPropertyDerivative<Real>("t_name", getVar("miu", 0)->name()))
{
}

Real
Overpotential::computeQpResidual()
{
  return (_u[_qp]+_T[_qp])*_test[_i][_qp]-_K[_qp]*_grad_odp[_qp]*_grad_test[_i][_qp];
}

Real
Overpotential::computeQpJacobian()
{
  return _phi[_j][_qp]*_test[_i][_qp];
}
Real
Overpotential::computeQpOffDiagJacobian(unsigned int jvar)
{
   
   if (jvar == _odp_var)
	   return -_K[_qp]*_grad_phi[_j][_qp]*_grad_test[_i][_qp]+_dTdeta[_qp]*_phi[_j][_qp]*_test[_i][_qp];
  else if (jvar == _ep_var)
     return _dTdpot[_qp]*_phi[_j][_qp]*_test[_i][_qp];
 else if (jvar == _miu_var)
     return _dTdmiu[_qp]*_phi[_j][_qp]*_test[_i][_qp];
  else
     return 0;
}
