//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"
// #include "SplitCHBase.h"
#include "JvarMapInterface.h"
#include "DerivativeMaterialInterface.h"

// Forward Declarations

/**
 * SplitCHBulk child class that takes all the necessary data from a
 * KKSBaseMaterial.
 * We calculate \f$ \frac{\partial F_a}{\partial c_a} \f$.
 * This takes advantage of the KKS identity
 *
 * \f$ dF/dc = dF_a/dc_a (= dF_b/dc_b) \f$
 *
 * The non-linear variable for this Kernel is the concentration 'c'.
 * The user picks one phase free energy \f$ F_a \f$ (f_base) and its corresponding
 * phase concentration \f$ c_a \f$
 */

class KKSSplitCHCRes : public DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>
// class KKSSplitCHCRes : public DerivativeMaterialInterface<JvarMapKernelInterface<SplitCHBase>>
{
public:
  static InputParameters validParams();

  KKSSplitCHCRes(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  unsigned int _num_j;
  const JvarMap & _c_map;
  unsigned int _num_c;
  int _o;
  /// Chemical potential
  unsigned int _w_var;
  const VariableValue & _w;
  unsigned int _eta_var;
  std::vector<MaterialPropertyName> _c1_names;
  std::vector<MaterialPropertyName> _dc1db_names;
  std::vector<std::vector<const MaterialProperty<Real> *>> _prop_dc1db;

  std::vector<MaterialPropertyName> _dc1deta_names;
  std::vector<const MaterialProperty<Real> *> _prop_dc1deta;

  const MaterialPropertyName _F1_name;
  std::vector<const MaterialProperty<Real> *> _prop_dF1dc1;
  std::vector<const MaterialProperty<Real> *> _prop_d2F1dc1db1;
};
