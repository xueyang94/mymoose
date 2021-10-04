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
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  const VariableValue & _c;

  const VariableValue & _eta;

  const MaterialProperty<Real> & _h;

  MaterialProperty<Real> & _c1;
  MaterialProperty<Real> & _c2;
  const MaterialProperty<Real> & _c1_old;
  const MaterialProperty<Real> & _c2_old;

  const Real _c1_initial;
  const Real _c2_initial;

  const Real _abs_tol;
  const Real _rel_tol;

  MaterialProperty<Real> & _dc1dc;
  MaterialProperty<Real> & _dc2dc;
  MaterialProperty<Real> & _dc1deta;
  MaterialProperty<Real> & _dc2deta;

  const SymbolName _c1_name;
  const SymbolName _c2_name;

  MaterialBase & _f1;
  MaterialBase & _f2;

  const MaterialProperty<Real> & _first_df1;
  const MaterialProperty<Real> & _first_df2;

  const MaterialProperty<Real> & _second_df1;
  const MaterialProperty<Real> & _second_df2;

  MaterialProperty<Real> & _iter;
  const unsigned int _min_iter;
  const unsigned int _max_iter;
};
