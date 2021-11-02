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
  // const std::vector<MaterialPropertyName> _ci_old;
  // const std::vector<MaterialProperty<Real> *> _ci_old_prop;

  const std::vector<Real> _ci_IC;

  std::vector<MaterialPropertyName> _dcidc_names;
  std::vector<MaterialProperty<Real> *> _prop_dcidc;
  std::vector<MaterialPropertyName> _coupled_dcidb_names;
  std::vector<MaterialProperty<Real> *> _prop_coupled_dcidb;

  std::vector<MaterialPropertyName> _dcidetaj_names;
  std::vector<std::vector<MaterialProperty<Real> *>> _prop_dcidetaj;

  std::vector<MaterialName> _Fi_material;
  // std::vector<MaterialBase *> _fi_material;
  MaterialBase & _f1;
  MaterialBase & _f2;
  MaterialBase & _f3;

  std::vector<MaterialPropertyName> _Fi_names;
  std::vector<const MaterialProperty<Real> *> _first_dFi;
  std::vector<const MaterialProperty<Real> *> _second_dFi;

  MaterialProperty<Real> & _iter;
  const Real _abs_tol;
  const Real _rel_tol;
  NestedSolve _nested_solve;
};
