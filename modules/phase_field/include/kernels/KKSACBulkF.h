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
 * KKSACBulkBase child class for the free energy difference term
 * \f$ -\frac{dh}{d\eta}(F_a-F_b)+w\frac{dg}{d\eta} \f$
 * in the the Allen-Cahn bulk residual.
 *
 * The non-linear variable for this Kernel is the order parameter 'eta'.
 */
// class KKSACBulkF : public Kernel
class KKSACBulkF : public KKSACBulkBase
{
public:
  static InputParameters validParams();

  KKSACBulkF(const InputParameters & parameters);

protected:
  virtual Real computeDFDOP(PFFunctionType type);
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const JvarMap & _c_map;
  const unsigned int _num_c;

  std::vector<MaterialPropertyName> _ci_names;
  std::vector<MaterialPropertyName> _dcidb_names;
  std::vector<std::vector<std::vector<const MaterialProperty<Real> *>>> _prop_dcidb;
  std::vector<MaterialPropertyName> _dcideta_names;
  std::vector<std::vector<const MaterialProperty<Real> *>> _prop_dcideta;
  const MaterialPropertyName _Fa_name;
  std::vector<const MaterialProperty<Real> *> _first_dFa;
  const MaterialPropertyName _Fb_name;
  const MaterialProperty<Real> & _prop_Fb;
  std::vector<const MaterialProperty<Real> *> _first_dFb;
  const MaterialProperty<Real> & _prop_dg;
  const MaterialProperty<Real> & _prop_d2g;
  const MaterialProperty<Real> & _L;
  /// double well height parameter
  Real _w;
};
