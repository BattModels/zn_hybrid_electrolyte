/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "DiffNP.h"

template <>
InputParameters
validParams<DiffNP>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription(
      "Add in additional terms grad(w)*grad(w)+Grad(w)*grad(phi)");
  params.addParam<MaterialPropertyName>("Q_name","The mobility for DiffNP used with the kernel for phi");
  params.addParam<MaterialPropertyName>("QM_name","The mobility for DiffNP used with th kernel for c");
  params.addRequiredCoupledVar(
      "cp", "applied potential");
  params.addRequiredCoupledVar(
      "cv", "coupled variable");
  return params;
}

DiffNP::DiffNP(const InputParameters & parameters)
: DerivativeMaterialInterface<Kernel>(parameters),
	_cp_var(coupled("cp")),
	_cp(coupledValue("cp")),
  _cv_var(coupled("cv")),
  _cv(coupledValue("cv")),
	_grad_cp(coupledGradient("cp")),
	_Q(getMaterialProperty<Real>("Q_name")),
	_QM(getMaterialProperty<Real>("QM_name")),
	_dQ(getMaterialPropertyDerivative<Real>("Q_name", _var.name())),
  _dQv(getMaterialPropertyDerivative<Real>("Q_name", getVar("cv", 0)->name())),
  _dQMv(getMaterialPropertyDerivative<Real>("QM_name", getVar("cv", 0)->name())),
  _dQM(getMaterialPropertyDerivative<Real>("QM_name", _var.name()))
{
}

Real
DiffNP::computeQpResidual()
{
    return  _Q[_qp]*_grad_u[_qp]*_grad_u[_qp]*_test[_i][_qp]
           +_QM[_qp]*_grad_cp[_qp]*_grad_u[_qp]*_test[_i][_qp];
}

Real
DiffNP::computeQpJacobian()
{
  return  _Q[_qp]*_grad_phi[_j][_qp]*_grad_u[_qp]*_test[_i][_qp]
          +_QM[_qp]*_grad_phi[_j][_qp]*_grad_cp[_qp]*_test[_i][_qp];
}
Real
DiffNP::computeQpOffDiagJacobian(unsigned int jvar)
{
   if (jvar == _cp_var)
	   return _QM[_qp]*_grad_phi[_j][_qp]*_grad_u[_qp]*_test[_i][_qp];
   else  if (jvar == _cv_var)
    return  _dQv[_qp]*_grad_u[_qp]*_grad_u[_qp]*_test[_i][_qp]*_phi[_j][_qp]
           +_dQMv[_qp]*_grad_cp[_qp]*_grad_u[_qp]*_test[_i][_qp]*_phi[_j][_qp];
   else
        return 0;
}
