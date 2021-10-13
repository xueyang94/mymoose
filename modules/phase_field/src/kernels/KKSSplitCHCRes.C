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
  // InputParameters params = DerivativeMaterialInterface<Kernel>::validParams();
  InputParameters params = SplitCHBase::validParams();
  params.addClassDescription(
      "KKS model kernel for the split Bulk Cahn-Hilliard term. This kernel operates on the "
      "physical concentration 'c' as the non-linear variable");
  params.addRequiredParam<MaterialPropertyName>("dc1dc_name", "The name of dc1/dc");
  params.addRequiredParam<MaterialPropertyName>("F1_name", "F1");
  params.addRequiredCoupledVar("w",
                               "Chemical potential non-linear helper variable for the split solve");
  params.addCoupledVar("etas", "Order parameters for all phases.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "dc1detaj_names", "The names of dc1/detaj in the order of dc1deta1, dc2deta1, dc3deta1, etc");
  params.addCoupledVar("args_a", "Vector of additional arguments to Fa");
  params.addRequiredCoupledVar("global_c", "The interpolated concentration");
  return params;
}

KKSSplitCHCRes::KKSSplitCHCRes(const InputParameters & parameters)
  // : DerivativeMaterialInterface<Kernel>(parameters),
  : DerivativeMaterialInterface<JvarMapKernelInterface<SplitCHBase>>(parameters),
    _dc1dc(getMaterialProperty<Real>("dc1dc_name")),
    _c1_name("c1"),
    _c2_name("c2"),
    _first_df1(getMaterialPropertyDerivative<Real>("F1_name", _c1_name)),
    _second_df1(getMaterialPropertyDerivative<Real>("F1_name", _c1_name, _c1_name)),
    _w_var(coupled("w")),
    _w(coupledValue("w")),
    _eta_names(coupledComponents("etas")),
    _eta_map(getParameterJvarMap("etas")),
    _num_j(_eta_names.size()),
    _dc1detaj_names(getParam<std::vector<MaterialPropertyName>>("dc1detaj_names")),
    _prop_dc1detaj(_num_j),
    _c_var(coupled("global_c"))
{
  for (unsigned int i = 0; i < _num_j; ++i)
    _prop_dc1detaj[i] = &getMaterialPropertyByName<Real>(_dc1detaj_names[i]);
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
  Real res = 0;

  // treat w variable explicitly
  if (jvar == _w_var)
    return -_phi[_j][_qp] * _test[_i][_qp];

  // if c is the coupled variable
  if (jvar == _c_var)
    return _phi[_j][_qp] * _test[_i][_qp] * _second_df1[_qp] * _dc1dc[_qp];

  // if order parameters are the coupled variables
  auto etavar = mapJvarToCvar(jvar, _eta_map);
  if (etavar >= 0)
    res = _phi[_j][_qp] * _test[_i][_qp] * _second_df1[_qp] * (*_prop_dc1detaj[etavar])[_qp];

  return res;
}
