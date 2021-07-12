//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KKSSplitCHCRes.h"

registerMooseObject("PhaseFieldApp", KKSSplitCHCRes);

InputParameters
KKSSplitCHCRes::validParams()
{
  InputParameters params = SplitCHBase::validParams();
  params.addClassDescription(
      "KKS model kernel for the split Bulk Cahn-Hilliard term. This kernel operates on the "
      "physical concentration 'c' as the non-linear variable");
  params.addRequiredParam<MaterialPropertyName>("c1_name", "The name of c1");
  params.addRequiredParam<MaterialPropertyName>("dc1dc_name", "The name of dc1/dc");
  params.addRequiredParam<MaterialPropertyName>("dc1deta_name", "The name of dc1/deta");
  params.addRequiredCoupledVar("w",
                               "Chemical potential non-linear helper variable for the split solve");
  return params;
}

KKSSplitCHCRes::KKSSplitCHCRes(const InputParameters & parameters)
  : DerivativeMaterialInterface<JvarMapKernelInterface<SplitCHBase>>(parameters),
    _c1(getMaterialProperty<Real>("c1_name")),
    _dc1dc(getMaterialProperty<Real>("dc1dc_name")),
    _dc1deta(getMaterialProperty<Real>("dc1deta_name")),
    _w_var(coupled("w")),
    _w(coupledValue("w"))

{
}

Real
KKSSplitCHCRes::computeQpResidual()
{
  return (400 * log(_c1[_qp]) - 400 * log(1 - _c1[_qp]) - 28 - _w[_qp]) * _test[_i][_qp];
}

Real
KKSSplitCHCRes::computeQpJacobian()
{
  return _dc1dc * (400 / _c1[_qp] + 400 / (1 - _c1[_qp])) * _phi[_j][_qp] * _test[_i][_qp];
}

Real
KKSSplitCHCRes::computeQpOffDiagJacobian(unsigned int jvar)
{
  // treat w variable explicitly
  if (jvar == _w_var)
    return -_phi[_j][_qp] * _test[_i][_qp];

  // eta is the couple variable
  return _dc1deta[_qp] * (400 / _c1[_qp] + 400 / (1 - _c1[_qp])) * _phi[_j][_qp] * _test[_i][_qp];
}
