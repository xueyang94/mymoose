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
 * KKSACBulkBase child class for the free energy difference term
 * \f$ -\frac{dh}{d\eta}(F_a-F_b)+w\frac{dg}{d\eta} \f$
 * in the the Allen-Cahn bulk residual.
 *
 * The non-linear variable for this Kernel is the order parameter 'eta'.
 */
class KKSACBulkF : public Kernel
{
public:
  static InputParameters validParams();

  KKSACBulkF(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /// double well height parameter
  Real _m;

  const MaterialProperty<Real> & _c1;
  const MaterialProperty<Real> & _c2;

  const MaterialProperty<Real> & _dc1dc;
  const MaterialProperty<Real> & _dc1deta;
  const MaterialProperty<Real> & _dc2dc;
  const MaterialProperty<Real> & _dc2deta;

  const MaterialProperty<Real> & _L;
  const MaterialProperty<Real> & _f1;
  const MaterialProperty<Real> & _f2;

  // Chemical potential
  unsigned int _w_var;
};
