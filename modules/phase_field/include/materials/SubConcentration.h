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
  const VariableValue & _eta;
  const MaterialProperty<Real> & _prop_h;

  std::vector<MaterialPropertyName> _ci_names;
  std::vector<MaterialProperty<Real> *> _prop_ci;

  const MaterialProperty<Real> & _c1_old;
  const MaterialProperty<Real> & _c2_old;
  const std::vector<Real> _ci_IC;

  std::vector<MaterialName> _Fi_material;
  MaterialBase & _f1;
  MaterialBase & _f2;

  std::vector<MaterialPropertyName> _Fi_names;
  std::vector<const MaterialProperty<Real> *> _first_dFi;
  std::vector<const MaterialProperty<Real> *> _second_dFi;

  MaterialProperty<Real> & _iter;
  const Real _abs_tol;
  const Real _rel_tol;
  NestedSolve _nested_solve;
};
