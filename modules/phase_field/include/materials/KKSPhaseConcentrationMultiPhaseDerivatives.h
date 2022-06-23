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

class KKSPhaseConcentrationMultiPhaseDerivatives : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  KKSPhaseConcentrationMultiPhaseDerivatives(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  const unsigned int _num_c;
  const std::vector<VariableName> _eta_names;
  const unsigned int _num_j;
  std::vector<const MaterialProperty<Real> *> _prop_ci;
  std::vector<MaterialPropertyName> _ci_names;
  std::vector<MaterialPropertyName> _dcidb_names;
  std::vector<std::vector<std::vector<MaterialProperty<Real> *>>> _prop_dcidb;
  std::vector<MaterialPropertyName> _dcidetaj_names;
  std::vector<std::vector<std::vector<MaterialProperty<Real> *>>> _prop_dcidetaj;
  std::vector<MaterialPropertyName> _hj_names;
  std::vector<const MaterialProperty<Real> *> _prop_hj;
  std::vector<std::vector<const MaterialProperty<Real> *>> _dhjdetap;
  std::vector<MaterialPropertyName> _Fj_names;
  std::vector<std::vector<std::vector<const MaterialProperty<Real> *>>> _d2Fjdcjdbj;
};
