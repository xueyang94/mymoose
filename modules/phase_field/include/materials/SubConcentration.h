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
#include "NestedSolve.h"

class SubConcentration : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  SubConcentration(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  const VariableValue & _c;

  const unsigned int _num_eta;
  const std::vector<VariableName> _eta_names;

  std::vector<MaterialPropertyName> _hj_names;
  std::vector<const MaterialProperty<Real> *> _prop_hj;
  std::vector<std::vector<const MaterialProperty<Real> *>> _prop_dhjdetai;

  std::vector<MaterialPropertyName> _ci_names;
  std::vector<MaterialProperty<Real> *> _ci_prop;

  const MaterialProperty<Real> & _c1_old;
  const MaterialProperty<Real> & _c2_old;
  const MaterialProperty<Real> & _c3_old;
  // const std::vector<MaterialProperty<Real> *> _ci_old;

  const std::vector<Real> _ci_IC;

  std::vector<MaterialPropertyName> _dcidc_names;
  std::vector<MaterialProperty<Real> *> _prop_dcidc;

  std::vector<MaterialPropertyName> _dcidetaj_names;
  std::vector<std::vector<MaterialProperty<Real> *>> _prop_dcidetaj;

  MaterialBase & _f1;
  MaterialBase & _f2;
  MaterialBase & _f3;

  const SymbolName _c1_name;
  const SymbolName _c2_name;
  const SymbolName _c3_name;

  const MaterialProperty<Real> & _first_df1;
  const MaterialProperty<Real> & _first_df2;
  const MaterialProperty<Real> & _first_df3;

  const MaterialProperty<Real> & _second_df1;
  const MaterialProperty<Real> & _second_df2;
  const MaterialProperty<Real> & _second_df3;

  MaterialProperty<Real> & _iter;
  const Real _abs_tol;
  const Real _rel_tol;
  NestedSolve _nested_solve;
};
