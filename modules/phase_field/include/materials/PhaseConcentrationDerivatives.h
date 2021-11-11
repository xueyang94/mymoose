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

class PhaseConcentrationDerivatives : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  PhaseConcentrationDerivatives(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  const unsigned int _num_c;
  const std::vector<VariableName> _c_names;
  std::vector<std::vector<const MaterialProperty<Real> *>> _prop_ci;

  const std::vector<VariableName> _eta_names;
  const unsigned int _num_eta;
  std::vector<MaterialPropertyName> _Fj_names;
  std::vector<std::vector<std::vector<const MaterialProperty<Real> *>>> _prop_d2Fjdcjdbj;
  std::vector<MaterialPropertyName> _hj_names;
  std::vector<const MaterialProperty<Real> *> _prop_hj;
  std::vector<std::vector<const MaterialProperty<Real> *>> _prop_dhjdetai;
  std::vector<MaterialPropertyName> _ci_names;

  std::vector<std::vector<MaterialPropertyName>> _ci_name_matrix;

  std::vector<MaterialPropertyName> _dcidb_names;
  std::vector<std::vector<std::vector<MaterialProperty<Real> *>>> _prop_dcidb;
  std::vector<MaterialPropertyName> _dcidetaj_names;
  std::vector<std::vector<std::vector<MaterialProperty<Real> *>>> _prop_dcidetaj;
};
