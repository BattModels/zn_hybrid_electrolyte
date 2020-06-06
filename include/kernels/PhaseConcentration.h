/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef PhaseConcentration_H
#define PhaseConcentration_H

#include "Kernel.h"
#include "DerivativeMaterialInterface.h"

// Forward Declarations
class PhaseConcentration;

template <>
InputParameters validParams<PhaseConcentration>();

/**
 * Enforce sum of phase concentrations to be the real concentration.
 *
 * \f$ c=h(\eta)c_a+\left(1-h(\eta)\right)c_b\f$
 *
 * The non-linear variable for this Kernel is the concentration \f$ c_b \f$, while
 * \f$ c_a \f$ and \f$ c \f$ are supplied as coupled variables.
 * (compare this to KKSPhaseChemicalPotential, where the non-linear variable is
 * the other phase concentration \f$ c_a \f$!)
 *
 * \see KKSPhaseChemicalPotential
 * \see KKSHEtaPolyMaterial
 */
class PhaseConcentration : public DerivativeMaterialInterface<Kernel>
{
public:
  PhaseConcentration(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const VariableValue & _w;
  unsigned int _w_var;

  const VariableValue & _eta;
  unsigned int _eta_var;
  const MaterialProperty<Real> & _cl;
  const MaterialProperty<Real> & _cs;
  const MaterialProperty<Real> & _dcl;
  const MaterialProperty<Real> & _dcs;
  /// Switching function \f$ h(\eta) \f$
  const MaterialProperty<Real> & _prop_h;

  /// Derivative of the switching function \f$ \frac d{d\eta} h(\eta) \f$
  const MaterialProperty<Real> & _prop_dh;
};

#endif // PhaseConcentration_H
