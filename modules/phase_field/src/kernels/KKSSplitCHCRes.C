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
  params.addCoupledVar("etas", "Order parameters for all phases.");
  params.addRequiredCoupledVar("coupled_global_b",
                               "The coupled interpolated concentration s, l, etc");
  params.addRequiredCoupledVar("w",
                               "Chemical potential non-linear helper variable for the split solve");
  params.addRequiredParam<MaterialPropertyName>("c1_name", "c1");
  params.addParam<std::vector<MaterialPropertyName>>(
      "coupled_b1",
      "The phase concentration of coupled_global_b in F1. The order must match coupled_global_b.");

  params.addRequiredParam<MaterialPropertyName>("dc1dc_name", "The name of dc1/dc");
  params.addParam<std::vector<MaterialPropertyName>>(
      "dc1db_names",
      "The derivatives of the phase concentration c1 (corresponding to the kernel global c) wrt "
      "coupled global concentrations. The order must match coupled_global_b.");
  params.addParam<std::vector<MaterialPropertyName>>(
      "db1dc_names",
      "The phase concentration of other_global_b in the frist phase of Fj_names taken derivative"
      "wrt kernel variable. The b1 order must match other_global_b.");
  params.addParam<std::vector<MaterialPropertyName>>(
      "db1db_names",
      "The phase concentration of coupled_global_b in F1 taken derivative of the corresponding "
      "coupled_global_b. The order must match coupled_global_b.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "dc1detaj_names",
      "The names of dc1/detaj must be in the same order as etas, for exemple, dc1deta1, dc2deta1, "
      "dc3deta1, etc");
  params.addParam<std::vector<MaterialPropertyName>>(
      "db1detaj_names",
      "The phase concentration of other_global_b in the first phase of Fj_names taken derivative "
      "wrt eta. b1 must match the order of other_global_b. etaj must match the order of etas.");
  params.addRequiredParam<MaterialPropertyName>("F1_name", "F1");
  return params;
}

KKSSplitCHCRes::KKSSplitCHCRes(const InputParameters & parameters)
  : DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>(parameters),
    // : DerivativeMaterialInterface<JvarMapKernelInterface<SplitCHBase>>(parameters),
    _eta_names(coupledComponents("etas")),
    _eta_map(getParameterJvarMap("etas")),
    _num_j(_eta_names.size()),
    _coupled_b_map(getParameterJvarMap("coupled_global_b")),
    _num_coupled_b(coupledComponents("coupled_global_b")),
    _w_var(coupled("w")),
    _w(coupledValue("w")),
    _c1_name(getParam<MaterialPropertyName>("c1_name")),
    _coupled_b1_names(getParam<std::vector<MaterialPropertyName>>("coupled_b1")),

    _dc1dc(getMaterialProperty<Real>("dc1dc_name")),
    _dc1db_names(getParam<std::vector<MaterialPropertyName>>("dc1db_names")),
    _prop_dc1db(_num_coupled_b),
    _db1dc_names(getParam<std::vector<MaterialPropertyName>>("db1dc_names")),
    _prop_db1dc(_num_coupled_b),
    _db1db_names(getParam<std::vector<MaterialPropertyName>>("db1db_names")),
    _prop_db1db(_num_coupled_b),
    _dc1detaj_names(getParam<std::vector<MaterialPropertyName>>("dc1detaj_names")),
    _prop_dc1detaj(_num_j),
    _db1detaj_names(getParam<std::vector<MaterialPropertyName>>("db1detaj_names")),
    _prop_db1detaj(_num_coupled_b),

    _F1_name(getParam<MaterialPropertyName>("F1_name")),
    _first_df1(getMaterialPropertyDerivative<Real>("F1_name", _c1_name)),
    _second_df1(getMaterialPropertyDerivative<Real>("F1_name", _c1_name, _c1_name)),
    _prop_d2F1dc1db1(_num_coupled_b)

{
  // initialize _prop_dc1db
  for (unsigned int i = 0; i < _num_coupled_b; ++i)
    _prop_dc1db[i] = &getMaterialPropertyByName<Real>(_dc1db_names[i]);

  // initialize _prop_db1dc
  for (unsigned int i = 0; i < _num_coupled_b; ++i)
    _prop_db1dc[i] = &getMaterialPropertyByName<Real>(_db1dc_names[i]);

  // initialize _prop_db1db
  for (unsigned int i = 0; i < _num_coupled_b; ++i)
    _prop_db1db[i] = &getMaterialPropertyByName<Real>(_db1db_names[i]);

  // initialize _prop_dc1detaj
  for (unsigned int i = 0; i < _num_j; ++i)
    _prop_dc1detaj[i] = &getMaterialPropertyByName<Real>(_dc1detaj_names[i]);

  // initialize _prop_db1detaj
  for (unsigned int m = 0; m < _num_coupled_b; ++m)
  {
    _prop_db1detaj[m].resize(_num_j);
    for (unsigned int n = 0; n < _num_j; ++n)
      _prop_db1detaj[m][n] = &getMaterialPropertyByName<Real>(_db1detaj_names[m * _num_j + n]);
  }

  // initialize _prop_d2F1dc1db1
  for (unsigned int i = 0; i < _num_coupled_b; ++i)
    _prop_d2F1dc1db1[i] =
        &getMaterialPropertyDerivative<Real>(_F1_name, _c1_name, _coupled_b1_names[i]);
}

Real
KKSSplitCHCRes::computeQpResidual()
{
  return (_first_df1[_qp] - _w[_qp]) * _test[_i][_qp];
}

Real
KKSSplitCHCRes::computeQpJacobian()
{
  Real factor = 0.0;
  factor = _second_df1[_qp] * _dc1dc[_qp];

  for (unsigned int n = 0; n < _num_coupled_b; ++n)
    factor += (*_prop_d2F1dc1db1[n])[_qp] * (*_prop_db1dc[n])[_qp];

  return factor * _phi[_j][_qp] * _test[_i][_qp];
}

Real
KKSSplitCHCRes::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real factor = 0.0;

  // treat w variable explicitly
  if (jvar == _w_var)
    return -_phi[_j][_qp] * _test[_i][_qp];

  // if b is the coupled variable
  auto coupled_b_var = mapJvarToCvar(jvar, _coupled_b_map);
  if (coupled_b_var >= 0)
  {
    factor = _second_df1[_qp] * (*_prop_dc1db[coupled_b_var])[_qp];

    for (unsigned int n = 0; n < _num_coupled_b; ++n)
      factor += (*_prop_d2F1dc1db1[n])[_qp] * (*_prop_db1db[n])[_qp];

    return factor * _phi[_j][_qp] * _test[_i][_qp];
  }

  // if order parameters are the coupled variables
  auto etavar = mapJvarToCvar(jvar, _eta_map);
  if (etavar >= 0)
  {
    factor = _second_df1[_qp] * (*_prop_dc1detaj[etavar])[_qp];

    for (unsigned int n = 0; n < _num_coupled_b; ++n)
      factor += (*_prop_d2F1dc1db1[n])[_qp] * (*_prop_db1detaj[n][etavar])[_qp];

    return factor * _phi[_j][_qp] * _test[_i][_qp];
  }

  return 0.0;
}
