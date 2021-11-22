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
  params.addCoupledVar("global_cs", "Global concentrations c, b, etc.");
  params.addCoupledVar("w", "Chemical potential non-linear helper variable for the split solve");
  params.addParam<std::vector<MaterialPropertyName>>(
      "c1_names",
      "Phase concentrations in the frist phase of etas. The order must match global_c, for "
      "example, c1, b1, etc.");
  params.addParam<std::vector<MaterialPropertyName>>(
      "dc1db_names",
      "The phase concentrations of the frist phase in Fj_names taken derivatives wrt global "
      "concentrations.  c1 and b must match the order of global_c. First keep the same c1 and loop "
      "through b for one species, for example, dc1dc, dc1db, db1dc, db1db");
  params.addParam<std::vector<MaterialPropertyName>>(
      "dc1detaj_names",
      "The phase concentrations of the first phase in Fj_names taken derivatives wrt etas. c1 must "
      "match the order in global_c, and etaj must match the order in etas. First keep the same c1 "
      "and loop through eta, for example, dc1deta1, "
      "dc1deta2, db1deta1, db1deta2.");
  params.addParam<MaterialPropertyName>("F1_name", "F1");
  return params;
}

KKSSplitCHCRes::KKSSplitCHCRes(const InputParameters & parameters)
  : DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>(parameters),
    // : DerivativeMaterialInterface<JvarMapKernelInterface<SplitCHBase>>(parameters),
    _eta_names(coupledComponents("etas")),
    _eta_map(getParameterJvarMap("etas")),
    _num_j(_eta_names.size()),
    _c_names(coupledComponents("global_cs")),
    _c_map(getParameterJvarMap("global_cs")),
    _num_c(coupledComponents("global_cs")),
    _o(-1), // position of nonlinear variable c in the list of global_cs
    _w_var(coupled("w")),
    _w(coupledValue("w")),
    _c1_names(getParam<std::vector<MaterialPropertyName>>("c1_names")),
    _dc1db_names(getParam<std::vector<MaterialPropertyName>>("dc1db_names")),
    _prop_dc1db(_num_c),

    _dc1detaj_names(getParam<std::vector<MaterialPropertyName>>("dc1detaj_names")),
    _prop_dc1detaj(_num_c),

    _F1_name(getParam<MaterialPropertyName>("F1_name")),
    _first_df1(getMaterialPropertyDerivative<Real>("F1_name", _c1_names[0])),
    _prop_d2F1dc1db1(_num_c)

{
  for (unsigned int i = 0; i < _num_c; ++i)
  {
    // get c names and variable indices
    _c_names[i] = getVar("global_cs", i)->name();

    // Set _o to the position of the nonlinear variable in the list of global_cs
    if (coupled("global_cs", i) == _var.number())
      _o = i;
  }

  // initialize _prop_dc1db
  for (unsigned int m = 0; m < _num_c; ++m)
  {
    _prop_dc1db[m].resize(_num_c);

    for (unsigned int n = 0; n < _num_c; ++n)
    {
      _prop_dc1db[m][n] = &getMaterialPropertyByName<Real>(_dc1db_names[m * _num_c + n]);
    }
  }

  // initialize _prop_dc1detaj
  for (unsigned int m = 0; m < _num_c; ++m)
  {
    _prop_dc1detaj[m].resize(_num_j);

    for (unsigned int n = 0; n < _num_j; ++n)
      _prop_dc1detaj[m][n] = &getMaterialPropertyByName<Real>(_dc1detaj_names[m * _num_j + n]);
  }

  // initialize _prop_d2F1dc1db1
  for (unsigned int i = 0; i < _num_c; ++i)
    _prop_d2F1dc1db1[i] =
        &getMaterialPropertyDerivative<Real>(_F1_name, _c1_names[0], _c1_names[i]);
}

Real
KKSSplitCHCRes::computeQpResidual()
{
  return (_first_df1[_qp] - _w[_qp]) * _test[_i][_qp];
}

Real
KKSSplitCHCRes::computeQpJacobian()
{
  Real sum = 0.0;

  for (unsigned int m = 0; m < _num_c; ++m)
    sum += (*_prop_d2F1dc1db1[m])[_qp] * (*_prop_dc1db[m][_o])[_qp];

  return sum * _phi[_j][_qp] * _test[_i][_qp];
}

Real
KKSSplitCHCRes::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real sum = 0.0;

  // treat w variable explicitly
  if (jvar == _w_var)
    return -_phi[_j][_qp] * _test[_i][_qp];

  // if b is the coupled variable
  auto compvar = mapJvarToCvar(jvar, _c_map);
  if (compvar >= 0)
  {
    // This can be further improve by not looping the nonlinear variable c because dRdc is
    // on-diagonal
    for (unsigned int m = 0; m < _num_c; ++m)
      sum += (*_prop_d2F1dc1db1[m])[_qp] * (*_prop_dc1db[m][compvar])[_qp];

    return sum * _phi[_j][_qp] * _test[_i][_qp];
  }

  // if order parameters are the coupled variables
  auto etavar = mapJvarToCvar(jvar, _eta_map);
  if (etavar >= 0)
  {
    for (unsigned int m = 0; m < _num_c; ++m)
      sum += (*_prop_d2F1dc1db1[m])[_qp] * (*_prop_dc1detaj[m][etavar])[_qp];

    return sum * _phi[_j][_qp] * _test[_i][_qp];
  }

  return 0.0;
}
