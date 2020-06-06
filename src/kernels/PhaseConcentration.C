/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "PhaseConcentration.h"

template <>
InputParameters
validParams<PhaseConcentration>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("KKS model kernel to enforce the decomposition of concentration into "
                             "phase concentration  (1-h(eta))*ca + h(eta)*cb - c = 0. The "
                             "non-linear variable of this kernel is c");
  params.addRequiredCoupledVar("w", "chemical potential");
  params.addRequiredCoupledVar("eta", "Phase a/b order parameter");
  params.addParam<MaterialPropertyName>(
      "cl_name", "negative of liquid concentration");
  params.addParam<MaterialPropertyName>(
      "cs_name", "negative of solid concentration");
  params.addParam<MaterialPropertyName>(
      "h_name", "h", "Base name for the switching function h(eta)"); // TODO: everywhere else this
                                                                     // is called just "h"
  return params;
}

// Phase interpolation func
PhaseConcentration::PhaseConcentration(const InputParameters & parameters)
  : DerivativeMaterialInterface<Kernel>(parameters),
    _w(coupledValue("w")),
    _w_var(coupled("w")),
    _eta(coupledValue("eta")),
    _eta_var(coupled("eta")),
    _cl(getMaterialProperty<Real>("cl_name")),
    _cs(getMaterialProperty<Real>("cs_name")),
    _dcl(getMaterialPropertyDerivative<Real>("cl_name", getVar("w", 0)->name())),
    _dcs(getMaterialPropertyDerivative<Real>("cs_name", getVar("w", 0)->name())),
    _prop_h(getMaterialProperty<Real>("h_name")),
    _prop_dh(getMaterialPropertyDerivative<Real>("h_name", getVar("eta", 0)->name()))
{
}

Real
PhaseConcentration::computeQpResidual()
{
  // R = -(1-h(eta))*ca - h(eta)*cb + c
  return _test[_i][_qp] * ((1.0 - _prop_h[_qp]) * _cl[_qp]+_prop_h[_qp] * _cs[_qp] + _u[_qp]);
}

Real
PhaseConcentration::computeQpJacobian()
{
  return _test[_i][_qp] * _phi[_j][_qp];
}

Real
PhaseConcentration::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _w_var)
    return _test[_i][_qp] * ((1.0 - _prop_h[_qp]) * _dcl[_qp]+_prop_h[_qp] * _dcs[_qp]) * _phi[_j][_qp];

  else if (jvar == _eta_var)
      return _test[_i][_qp] * (_cs[_qp] - _cl[_qp]) * _prop_dh[_qp] * _phi[_j][_qp];
  else
      return 0.0;
}
