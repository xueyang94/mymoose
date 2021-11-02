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
  InputParameters params = JvarMapKernelInterface<Kernel>::validParams();
  // InputParameters params = SplitCHBase::validParams();
  params.addClassDescription(
      "KKS model kernel for the split Bulk Cahn-Hilliard term. This kernel operates on the "
      "physical concentration 'c' as the non-linear variable");
  params.addRequiredParam<MaterialPropertyName>("dc1dc_name", "The name of dc1/dc");
  params.addRequiredParam<MaterialPropertyName>("c1_name", "c1");
  params.addRequiredParam<MaterialPropertyName>("F1_name", "F1");
  params.addRequiredCoupledVar("w",
                               "Chemical potential non-linear helper variable for the split solve");
  params.addCoupledVar("etas", "Order parameters for all phases.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "dc1detaj_names",
      "The names of dc1/detaj must be in the same order as etas, for exemple, dc1deta1, dc2deta1, "
      "dc3deta1, etc");
  params.addParam<MaterialPropertyName>(
      "dc1db_name",
      "The derivative of the phase concentration c1 (corresponding to the kernel global c) wrt a "
      "coupled global concentration b");
  params.addRequiredCoupledVar("coupled_global_b",
                               "The coupled interpolated concentration s, l, etc");
  return params;
}

KKSSplitCHCRes::KKSSplitCHCRes(const InputParameters & parameters)
  : DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>(parameters),
    // : DerivativeMaterialInterface<JvarMapKernelInterface<SplitCHBase>>(parameters),
    _dc1dc(getMaterialProperty<Real>("dc1dc_name")),
    _c1_name(getParam<MaterialPropertyName>("c1_name")),
    _first_df1(getMaterialPropertyDerivative<Real>("F1_name", _c1_name)),
    _second_df1(getMaterialPropertyDerivative<Real>("F1_name", _c1_name, _c1_name)),
    _w_var(coupled("w")),
    _w(coupledValue("w")),
    _eta_names(coupledComponents("etas")),
    _eta_map(getParameterJvarMap("etas")),
    _num_j(_eta_names.size()),
    _dc1detaj_names(getParam<std::vector<MaterialPropertyName>>("dc1detaj_names")),
    _prop_dc1detaj(_num_j),

    _dc1db(getMaterialProperty<Real>("dc1db_name")),
    _coupled_c_names(coupledComponents("coupled_global_b")),
    _coupled_c_map(getParameterJvarMap("coupled_global_b")),
    _num_coupled_c(_coupled_c_names.size())

{
  for (unsigned int i = 0; i < _num_j; ++i)
    _prop_dc1detaj[i] = &getMaterialPropertyByName<Real>(_dc1detaj_names[i]);
}

Real
KKSSplitCHCRes::computeQpResidual()
{
  // std::cout << "_first_df1 is " << _first_df1[_qp] << std::endl;
  // std::cout << "_second_df1 is " << _second_df1[_qp] << std::endl;

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
  // std::cout << "KKSSplitCHCRes dc1db or db1dc is " << _dc1db[_qp] << std::endl;

  // treat w variable explicitly
  if (jvar == _w_var)
    return -_phi[_j][_qp] * _test[_i][_qp];

  // if b is the coupled variable
  auto compvar = mapJvarToCvar(jvar, _coupled_c_map);
  if (compvar >= 0)
  {
    Real sum = 0.0;

    for (unsigned int n = 0; n < _num_coupled_c; ++n)
      sum += _second_df1[_qp] * _dc1db[_qp];

    return sum * _phi[_j][_qp] * _test[_i][_qp];
  }

  // if order parameters are the coupled variables
  auto etavar = mapJvarToCvar(jvar, _eta_map);
  if (etavar >= 0)
    return _second_df1[_qp] * (*_prop_dc1detaj[etavar])[_qp] * _phi[_j][_qp] * _test[_i][_qp];

  return 0.0;
}
