//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KKSMultiACBulkBase.h"

// Forward Declarations

/**
 * KKSMultiACBulkBase child class for the free energy term
 * \f$ \sum_j \frac{\partial h_j}{\partial \eta_i} F_j + w_i \frac{dg}{d\eta_i} \f$
 * in the the Allen-Cahn bulk residual.
 *
 * The non-linear variable for this Kernel is the order parameter \f$ eta_i \f$.
 */
class KKSMultiACBulkF : public KKSMultiACBulkBase
{
public:
  static InputParameters validParams();

  KKSMultiACBulkF(const InputParameters & parameters);

protected:
  virtual Real computeDFDOP(PFFunctionType type);
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  std::vector<MaterialPropertyName> _ci_names;
  std::vector<VariableName> _eta_names;
  const JvarMap & _eta_map;

  /// Position of the nonlinear variable in the list of cj's
  int _k;

  std::vector<std::vector<const MaterialProperty<Real> *>> _prop_d2hjdetapdetai;

  std::vector<MaterialPropertyName> _dcidc_names;
  std::vector<const MaterialProperty<Real> *> _prop_dcidc;

  std::vector<MaterialPropertyName> _dcidetaj_names;
  std::vector<std::vector<const MaterialProperty<Real> *>> _prop_dcidetaj;

  std::vector<const MaterialProperty<Real> *> _prop_dFidci;

  /// double well height parameter
  Real _wi;

  MaterialPropertyName _gi_name;
  const MaterialProperty<Real> & _prop_dgi;
  // std::vector<const MaterialProperty<Real> *> _prop_d2gpdetapdetai;
  const MaterialProperty<Real> & _prop_d2gi;

  unsigned int _c_var;
};
