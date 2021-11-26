//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KKSACBulkBase.h"

// Forward Declarations

/**
 * KKSACBulkBase child class for the phase concentration difference term
 * \f$ \frac{dh}{d\eta}\frac{dF_a}{dc_a}(c_a-c_b) \f$
 * in the the Allen-Cahn bulk residual.
 *
 * The non-linear variable for this Kernel is the order parameter 'eta'.
 */

class KKSACBulkC : public KKSACBulkBase
{
public:
  static InputParameters validParams();

  KKSACBulkC(const InputParameters & parameters);

protected:
  virtual Real computeDFDOP(PFFunctionType type);
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  std::vector<VariableName> _c_names;
  const JvarMap & _c_map;
  const unsigned int _num_c;
  const unsigned int _num_j;
  std::vector<MaterialPropertyName> _ci_names;
  std::vector<std::vector<MaterialPropertyName>> _ci_name_matrix;
  std::vector<std::vector<const MaterialProperty<Real> *>> _prop_ci;
  std::vector<MaterialPropertyName> _dcidb_names;
  std::vector<std::vector<std::vector<const MaterialProperty<Real> *>>> _prop_dcidb;
  std::vector<MaterialPropertyName> _dcideta_names;
  std::vector<std::vector<const MaterialProperty<Real> *>> _prop_dcideta;
  const MaterialProperty<Real> & _L;

  const MaterialPropertyName _Fa_name;
  std::vector<const MaterialProperty<Real> *> _first_dFa;
  std::vector<std::vector<const MaterialProperty<Real> *>> _d2F1dc1db1;
};
