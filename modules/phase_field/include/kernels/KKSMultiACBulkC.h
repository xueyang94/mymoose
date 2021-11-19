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
 * KKSACBulkBase child class for the phase concentration term
 * \f$ - \sum_j \frac{dF_1}{dc_1} \frac{dh_j}{d\eta_i} (c_j) \f$
 * in the the Allen-Cahn bulk residual.
 *
 * The non-linear variable for this Kernel is the order parameter 'eta_i'.
 */
class KKSMultiACBulkC : public KKSMultiACBulkBase
{
public:
  static InputParameters validParams();

  KKSMultiACBulkC(const InputParameters & parameters);

protected:
  virtual Real computeDFDOP(PFFunctionType type);
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  std::vector<VariableName> _eta_names;
  const JvarMap & _eta_map;
  /// Position of the nonlinear variable in the list of cj's
  int _k;
  unsigned int _c_var;
  std::vector<MaterialPropertyName> _ci_names;
  std::vector<const MaterialProperty<Real> *> _prop_ci;
  const unsigned int _num_other_b;
  const JvarMap & _other_b_map;
  std::vector<MaterialPropertyName> _other_b1_names;

  std::vector<MaterialPropertyName> _dcidc_names;
  std::vector<const MaterialProperty<Real> *> _prop_dcidc;
  std::vector<MaterialPropertyName> _dcidb_names;
  std::vector<std::vector<const MaterialProperty<Real> *>> _prop_dcidb;
  std::vector<MaterialPropertyName> _db1dc_names;
  std::vector<const MaterialProperty<Real> *> _prop_db1dc;
  std::vector<MaterialPropertyName> _db1db_names;
  std::vector<const MaterialProperty<Real> *> _prop_db1db;
  std::vector<MaterialPropertyName> _dcidetaj_names;
  std::vector<std::vector<const MaterialProperty<Real> *>> _prop_dcidetaj;
  std::vector<MaterialPropertyName> _db1detaj_names;
  std::vector<std::vector<const MaterialProperty<Real> *>> _prop_db1detaj;

  std::vector<std::vector<const MaterialProperty<Real> *>> _prop_d2hjdetaidetap;
  const MaterialProperty<Real> & _first_df1;
  const MaterialProperty<Real> & _second_df1;
  std::vector<const MaterialProperty<Real> *> _d2F1dc1db1;
};
