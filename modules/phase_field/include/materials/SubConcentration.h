//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "DerivativeMaterialInterface.h"

class SubConcentration : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  SubConcentration(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  const VariableValue & _c;

  const unsigned int _num_eta;
  const std::vector<VariableName> _eta_names;

  std::vector<MaterialPropertyName> _hj_names;
  std::vector<const MaterialProperty<Real> *> _prop_hj;
  std::vector<std::vector<const MaterialProperty<Real> *>> _prop_dhjdetai;

  // MaterialProperty<Real> & _c1;
  // MaterialProperty<Real> & _c2;
  // MaterialProperty<Real> & _c3;

  std::vector<MaterialPropertyName> _ci_names;
  std::vector<MaterialProperty<Real> *> _ci_prop;

  const std::vector<Real> _ci_IC;

  MaterialProperty<Real> & _dc1dc;
  MaterialProperty<Real> & _dc2dc;
  MaterialProperty<Real> & _dc3dc;

  MaterialProperty<Real> & _dc1deta1;
  MaterialProperty<Real> & _dc1deta2;
  MaterialProperty<Real> & _dc1deta3;
  MaterialProperty<Real> & _dc2deta1;
  MaterialProperty<Real> & _dc2deta2;
  MaterialProperty<Real> & _dc2deta3;
  MaterialProperty<Real> & _dc3deta1;
  MaterialProperty<Real> & _dc3deta2;
  MaterialProperty<Real> & _dc3deta3;

  const SymbolName _c1_name;
  const SymbolName _c2_name;
  const SymbolName _c3_name;

  MaterialBase & _f1;
  MaterialBase & _f2;
  MaterialBase & _f3;

  const MaterialProperty<Real> & _first_df1;
  const MaterialProperty<Real> & _first_df2;
  const MaterialProperty<Real> & _first_df3;

  const MaterialProperty<Real> & _second_df1;
  const MaterialProperty<Real> & _second_df2;
  const MaterialProperty<Real> & _second_df3;

  const Real _abs_tol;
  const Real _rel_tol;
  const unsigned int _maxiter;
};
