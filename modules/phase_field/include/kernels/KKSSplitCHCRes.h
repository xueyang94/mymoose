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

  const MaterialProperty<Real> & _dc1dc;

  const MaterialPropertyName _c1_name;

  const MaterialProperty<Real> & _first_df1;
  const MaterialProperty<Real> & _second_df1;

  /// Chemical potential
  unsigned int _w_var;
  const VariableValue & _w;

  std::vector<VariableName> _eta_names;
  const JvarMap & _eta_map;
  unsigned int _num_j;
  std::vector<MaterialPropertyName> _dc1detaj_names;
  std::vector<const MaterialProperty<Real> *> _prop_dc1detaj;

  const MaterialProperty<Real> & _dc1db;
  std::vector<VariableName> _coupled_c_names; //
  const JvarMap & _coupled_c_map;             //
  unsigned int _num_coupled_c;                //
};
