/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef MATDIFFUSIONBASENP_H
#define MATDIFFUSIONBASENP_H

#include "Kernel.h"
#include "JvarMapInterface.h"
#include "DerivativeMaterialInterface.h"

/**
 * This class template implements a diffusion kernel with a mobility that can vary
 * spatially and can depend on variables in the simulation. Two classes are derived from
 * this template, MatDiffusion for isotropic diffusion and MatAnisoDiffusion for
 * anisotropic diffusion.
 *
 * \tparam T Type of the diffusion coefficient parameter. This can be Real for
 *           isotropic diffusion or RealTensorValue for the general anisotropic case.
 */
template <typename T>
class MatDiffusionBaseNP : public DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>
{
public:
  MatDiffusionBaseNP(const InputParameters & parameters);

  virtual void initialSetup();

  /// in class templates this function has to be a static member
  static InputParameters validParams();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  virtual Real computeQpCJacobian();

  /// diffusion coefficient
  const MaterialProperty<T> & _D;
  const MaterialProperty<T> & _DN;
  /// diffusion coefficient derivative w.r.t. the kernel variable
  const MaterialProperty<T> & _dDdc;
  const MaterialProperty<T> & _dDNdc;

  /// diffusion coefficient derivatives w.r.t. coupled variables
  std::vector<const MaterialProperty<T> *> _dDdarg;

  /// is the kernel used in a coupled form?
  const bool _is_coupled;

  /// int label for the Concentration
  unsigned int _conc_var;

  /// Gradient of the concentration
  const VariableGradient & _grad_conc;
    /// int label for the Concentration
  unsigned int _pot_var;
  const VariableValue & _pot;
  unsigned int _cve_var;
  const VariableValue & _cve;
  const MaterialProperty<Real> & _dDdcve;
  const MaterialProperty<Real> & _dDNdcve;
  /// Gradient of the concentration
  const VariableGradient & _grad_pot;
};

template <typename T>
InputParameters
MatDiffusionBaseNP<T>::validParams()
{
  InputParameters params = ::validParams<Kernel>();
  params.addParam<MaterialPropertyName>("D_name", "D", "The name of the diffusivity");
  params.addParam<MaterialPropertyName>("DN_name", "DN", "The name of the Nerst-Plank diffusivity");
  params.addCoupledVar("args", "Vector of arguments of the diffusivity");
  params.addCoupledVar("conc",
                       "Coupled concentration variable for kernel to operate on; if this "
                       "is not specified, the kernel's nonlinear variable will be used as "
                       "usual");
   params.addRequiredCoupledVar("pot", "electric potential");
   params.addRequiredCoupledVar("cve", "coupled variable eta");
  return params;
}

template <typename T>
MatDiffusionBaseNP<T>::MatDiffusionBaseNP(const InputParameters & parameters)
  : DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>(parameters),
    _D(getMaterialProperty<T>("D_name")),
	_DN(getMaterialProperty<T>("DN_name")),
    _dDdc(getMaterialPropertyDerivative<T>("D_name", _var.name())),
	_dDNdc(getMaterialPropertyDerivative<T>("DN_name", _var.name())),
    _dDdarg(_coupled_moose_vars.size()),
    _is_coupled(isCoupled("conc")),
    _conc_var(_is_coupled ? coupled("conc") : _var.number()),
    _grad_conc(_is_coupled ? coupledGradient("conc") : _grad_u),
	_pot_var(coupled("pot")),
	_pot(coupledValue("pot")),
  _cve_var(coupled("cve")),
	_cve(coupledValue("cve")),
  _dDdcve(getMaterialPropertyDerivative<Real>("D_name", getVar("cve", 0)->name())),
	_dDNdcve(getMaterialPropertyDerivative<Real>("DN_name", getVar("cve", 0)->name())),
	_grad_pot(coupledGradient("pot"))
{
  // fetch derivatives
  for (unsigned int i = 0; i < _dDdarg.size(); ++i)
    _dDdarg[i] = &getMaterialPropertyDerivative<T>("D_name", _coupled_moose_vars[i]->name());
}

template <typename T>
void
MatDiffusionBaseNP<T>::initialSetup()
{
  validateNonlinearCoupling<Real>("D_name");
  validateNonlinearCoupling<Real>("DN_name");
}

template <typename T>
Real
MatDiffusionBaseNP<T>::computeQpResidual()
{
  return _D[_qp] * _grad_conc[_qp] * _grad_test[_i][_qp]+_DN[_qp] * _grad_pot[_qp] * _grad_test[_i][_qp];
}

template <typename T>
Real
MatDiffusionBaseNP<T>::computeQpJacobian()
{
  Real sum = _phi[_j][_qp] * _dDdc[_qp] * _grad_conc[_qp] * _grad_test[_i][_qp];
  sum +=_phi[_j][_qp] * _dDNdc[_qp] * _grad_pot[_qp] * _grad_test[_i][_qp];
  if (!_is_coupled)
    sum += computeQpCJacobian();

  return sum;
}

template <typename T>
Real
MatDiffusionBaseNP<T>::computeQpOffDiagJacobian(unsigned int jvar)
{
  // get the coupled variable jvar is referring to
  const unsigned int cvar = mapJvarToCvar(jvar);

  Real sum = (*_dDdarg[cvar])[_qp] * _phi[_j][_qp] * _grad_conc[_qp] * _grad_test[_i][_qp];
  if (jvar==_conc_var)
    sum += computeQpCJacobian();

  return sum;
  if (jvar==_pot_var)
 return _DN[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
 if (jvar==_cve_var)
 return _dDNdcve[_qp] * _grad_pot[_qp] * _grad_test[_i][_qp]*_phi[_j][_qp]
 +_dDdcve[_qp] * _grad_conc[_qp] * _grad_test[_i][_qp]*_phi[_j][_qp];
}

template <typename T>
Real
MatDiffusionBaseNP<T>::computeQpCJacobian()
{
  return _D[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}
#endif // MATDIFFUSIONBASENP_H
