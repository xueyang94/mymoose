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

class KKSPhaseConcentrationMaterial : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  KKSPhaseConcentrationMaterial(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  const std::vector<const VariableValue *> _prop_c;
  const unsigned int _num_c;
  const MaterialProperty<Real> & _prop_h;
  std::vector<MaterialPropertyName> _ci_names;
  std::vector<MaterialProperty<Real> *> _prop_ci;
  std::vector<const MaterialProperty<Real> *> _ci_old;
  std::vector<Real> _ci_IC;

  MaterialBase & _f1;
  MaterialBase & _f2;

  std::vector<MaterialPropertyName> _Fi_names;
  std::vector<std::vector<const MaterialProperty<Real> *>> _first_dFi;
  std::vector<std::vector<std::vector<const MaterialProperty<Real> *>>> _second_dFi;

  MaterialProperty<Real> & _iter;
  const Real _abs_tol;
  const Real _rel_tol;
  NestedSolve _nested_solve;
};
