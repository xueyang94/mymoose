//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// #include "Material.h"
#include "DerivativeMaterialInterface.h"
// #include "DerivativeMaterialPropertyNameInterface.h"

// class SubConcentration : public Material
class SubConcentration : public DerivativeMaterialInterface<Material>
// class SubConcentration : public DerivativeMaterialPropertyNameInterface<Material>
{
public:
  static InputParameters validParams();

  SubConcentration(const InputParameters & parameters);

protected:
  // virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  const VariableValue & _c;
  // std::string _c_name;

  const VariableValue & _eta;
  // std::string _eta_name;

  const MaterialProperty<Real> & _h;

  MaterialProperty<Real> & _c1;
  MaterialProperty<Real> & _c2;
  // const MaterialProperty<Real> & _c1_old;
  // const MaterialProperty<Real> & _c2_old;

  const Real _abs_tol;

  const Real _rel_tol;

  const unsigned int _maxiter;

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
};
