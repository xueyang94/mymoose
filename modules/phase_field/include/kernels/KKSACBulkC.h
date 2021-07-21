//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// #include "KKSACBulkBase.h"
#include "Kernel.h"

// Forward Declarations

/**
 * KKSACBulkBase child class for the phase concentration difference term
 * \f$ \frac{dh}{d\eta}\frac{dF_a}{dc_a}(c_a-c_b) \f$
 * in the the Allen-Cahn bulk residual.
 *
 * The non-linear variable for this Kernel is the order parameter 'eta'.
 */
// class KKSACBulkC : public KKSACBulkBase
class KKSACBulkC : public Kernel
{
public:
  static InputParameters validParams();

  KKSACBulkC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const MaterialProperty<Real> & _c1;
  const MaterialProperty<Real> & _c2;
  const MaterialProperty<Real> & _dc1dc;
  const MaterialProperty<Real> & _dc1deta;
  const MaterialProperty<Real> & _dc2dc;
  const MaterialProperty<Real> & _dc2deta;

  const MaterialProperty<Real> & _L;

  // Chemical potential
  unsigned int _w_var;
  const VariableValue & _w;
};
