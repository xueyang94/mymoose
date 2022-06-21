//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADMortarConstraint.h"

class GapHeatConductanceTest : public ADMortarConstraint
{
public:
  static InputParameters validParams();

  GapHeatConductanceTest(const InputParameters & parameters);

  using ADMortarConstraint::computeJacobian;
  void computeJacobian(Moose::MortarType mortar_type) override;

protected:
  ADReal computeQpResidual(Moose::MortarType mortar_type) final;

  const ADMaterialProperty<Real> & _secondary_gap_conductance;
  const ADMaterialProperty<Real> & _primary_gap_conductance;
  const FunctionTempl<Real> & _secondary_mms_function;
  const FunctionTempl<Real> & _primary_mms_function;
};
