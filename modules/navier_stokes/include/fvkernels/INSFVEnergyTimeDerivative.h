//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVTimeKernel.h"

class INSFVEnergyTimeDerivative : public FVTimeKernel
{
public:
  static InputParameters validParams();
  INSFVEnergyTimeDerivative(const InputParameters & params);

protected:
  ADReal computeQpResidual() override;

  /// the density
  const Moose::Functor<ADReal> & _rho;
  /// the heat conductivity
  const Moose::Functor<ADReal> & _cp;
  /// the time derivative of the heat conductivity
  const Moose::Functor<ADReal> & _cp_dot;
};
