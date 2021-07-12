//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

class SubConcentration : public Material
{
public:
  static InputParameters validParams();

  SubConcentration(const InputParameters & parameters);

protected:
  // virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  const VariableValue & _c;
  const MaterialProperty<Real> & _c_mat;
  const MaterialProperty<Real> & _c_mat_old;
  const MaterialProperty<Real> & _eta_mat;
  const MaterialProperty<Real> & _eta_mat_old;
  const MaterialProperty<Real> & _h;

  MaterialProperty<Real> & _c1;
  MaterialProperty<Real> & _c2;
  const MaterialProperty<Real> & _c1_old;
  const MaterialProperty<Real> & _c2_old;

  const Real _abs_tol;

  const Real _rel_tol;

  const unsigned int _maxiter;

  MaterialProperty<Real> & _dc1dc;
  MaterialProperty<Real> & _dc1deta;
  MaterialProperty<Real> & _dc1dc;
  MaterialProperty<Real> & _dc1deta;
};
