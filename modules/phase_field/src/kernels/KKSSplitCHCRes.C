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
  InputParameters params = Kernel::validParams();
  params.addClassDescription(
      "KKS model kernel for the split Bulk Cahn-Hilliard term. This kernel operates on the "
      "physical concentration 'c' as the non-linear variable");
  // params.addRequiredCoupledVar("global_c", "The interpolated concentration");
  // params.addRequiredCoupledVar("eta", "The order parameter");
  params.addRequiredParam<MaterialPropertyName>("dc1dc_name", "The name of dc1/dc");
  params.addRequiredParam<MaterialPropertyName>("dc1deta_name", "The name of dc1/deta");
  params.addRequiredParam<MaterialPropertyName>("df1dc1_name",
                                                "The name of the first derivative of f1 w.r.t. c1");
  params.addRequiredParam<MaterialPropertyName>(
      "d2f1dc1_name", "The name of the second derivative of f1 w.r.t. c1");
  params.addRequiredCoupledVar("w",
                               "Chemical potential non-linear helper variable for the split solve");
  return params;
}

KKSSplitCHCRes::KKSSplitCHCRes(const InputParameters & parameters)
  : Kernel(parameters),
    // : DerivativeMaterialInterface<Kernel>(parameters),
    // _c(coupledValue("global_c")),
    // _c_name(getVar("global_c", 0)->name()),
    // _eta(coupledValue("eta")),
    // _eta_name(getVar("eta", 0)->name()),
    // _dc1dc(getMaterialPropertyDerivative<Real>("c1_name", _c_name)),
    // _dc1deta(getMaterialPropertyDerivative<Real>("c1_name", _eta_name)),
    _dc1dc(getMaterialProperty<Real>("dc1dc_name")),
    _dc1deta(getMaterialProperty<Real>("dc1deta_name")),
    _first_df1(getMaterialProperty<Real>("df1dc1_name")),
    _second_df1(getMaterialProperty<Real>("d2f1dc1_name")),
    _w_var(coupled("w")),
    _w(coupledValue("w"))

{
}

Real
KKSSplitCHCRes::computeQpResidual()
{
  return (_first_df1[_qp] - _w[_qp]) * _test[_i][_qp];
}

Real
KKSSplitCHCRes::computeQpJacobian()
{
  return _second_df1[_qp] * _dc1dc[_qp] * _phi[_j][_qp] * _test[_i][_qp];
}

Real
KKSSplitCHCRes::computeQpOffDiagJacobian(unsigned int jvar)
{
  // treat w variable explicitly
  if (jvar == _w_var)
    return -_phi[_j][_qp] * _test[_i][_qp];

  // eta is the coupled variable
  return _second_df1[_qp] * _dc1deta[_qp] * _phi[_j][_qp] * _test[_i][_qp];
}
