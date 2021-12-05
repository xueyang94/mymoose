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
  params.addCoupledVar("global_cs", "Global concentrations c, b, etc.");
  params.addCoupledVar("w", "Chemical potential non-linear helper variable for the split solve");
  params.addCoupledVar("eta", "Order parameter");
  params.addParam<std::vector<MaterialPropertyName>>(
      "c1_names",
      "Phase concentrations in F1. The order must match global_c, for "
      "example, c1, b1, etc.");
  params.addParam<std::vector<MaterialPropertyName>>(
      "dc1db_names",
      "The phase concentrations in F1 taken derivatives wrt global "
      "concentrations.  c1 and b must match the order of global_c. First keep the same c1 and loop "
      "through b for one species, for example, dc1dc, dc1db, db1dc, db1db");
  params.addParam<std::vector<MaterialPropertyName>>(
      "dc1deta_names",
      "The phase concentrations in F1 taken derivatives wrt eta. The order must match c1_names, "
      "for example, dc1deta, db1deta, etc.");
  params.addParam<MaterialPropertyName>("F1_name", "F1");
  return params;
}

KKSSplitCHCRes::KKSSplitCHCRes(const InputParameters & parameters)
  : DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>(parameters),
    // : DerivativeMaterialInterface<JvarMapKernelInterface<SplitCHBase>>(parameters),
    _c_map(getParameterJvarMap("global_cs")),
    _num_c(coupledComponents("global_cs")),
    _o(-1), // position of nonlinear variable c in the list of global_cs
    _w_var(coupled("w")),
    _w(coupledValue("w")),
    _eta_var(coupled("eta")),
    _c1_names(getParam<std::vector<MaterialPropertyName>>("c1_names")),
    _dc1db_names(getParam<std::vector<MaterialPropertyName>>("dc1db_names")),
    _prop_dc1db(_num_c),

    _dc1deta_names(getParam<std::vector<MaterialPropertyName>>("dc1deta_names")),
    _prop_dc1deta(_num_c),

    _F1_name(getParam<MaterialPropertyName>("F1_name")),
    _prop_dF1dc1(_num_c),
    _prop_d2F1dc1db1(_num_c)

{
  for (unsigned int i = 0; i < _num_c; ++i)
  {
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

  // initialize _prop_dc1deta
  for (unsigned int m = 0; m < _num_c; ++m)
    _prop_dc1deta[m] = &getMaterialPropertyByName<Real>(_dc1deta_names[m]);

  for (unsigned int i = 0; i < _num_c; ++i)
  {
    // initialize _prop_dF1dc1
    _prop_dF1dc1[i] = &getMaterialPropertyDerivative<Real>(_F1_name, _c1_names[i]);

    // initialize _prop_d2F1dc1db1
    _prop_d2F1dc1db1[i] =
        &getMaterialPropertyDerivative<Real>(_F1_name, _c1_names[_o], _c1_names[i]);
  }
}

Real
KKSSplitCHCRes::computeQpResidual()
{
  return ((*_prop_dF1dc1[_o])[_qp] - _w[_qp]) * _test[_i][_qp];
}

Real
KKSSplitCHCRes::computeQpJacobian()
{
  Real sum = 0.0;

  // for (unsigned int m = 0; m < _num_c; ++m)
  //   sum += (*_prop_d2F1dc1db1[m])[_qp] * (*_prop_dc1db[m][_o])[_qp];

  sum = (*_prop_d2F1dc1db1[0])[_qp] * (*_prop_dc1db[0][_o])[_qp] +
        (*_prop_d2F1dc1db1[1])[_qp] * (*_prop_dc1db[1][_o])[_qp];

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
    // for (unsigned int m = 0; m < _num_c; ++m)
    //   sum += (*_prop_d2F1dc1db1[m])[_qp] * (*_prop_dc1db[m][compvar])[_qp];

    sum = (*_prop_d2F1dc1db1[0])[_qp] * (*_prop_dc1db[0][compvar])[_qp] +
          (*_prop_d2F1dc1db1[1])[_qp] * (*_prop_dc1db[1][compvar])[_qp];

    return sum * _phi[_j][_qp] * _test[_i][_qp];
  }

  // if order parameter is the coupled variables
  if (jvar == _eta_var)
  {
    // for (unsigned int m = 0; m < _num_c; ++m)
    //   sum += (*_prop_d2F1dc1db1[m])[_qp] * (*_prop_dc1deta[m])[_qp];

    sum = (*_prop_d2F1dc1db1[0])[_qp] * (*_prop_dc1deta[0])[_qp] +
          (*_prop_d2F1dc1db1[1])[_qp] * (*_prop_dc1deta[1])[_qp];

    return sum * _phi[_j][_qp] * _test[_i][_qp];
  }

  return 0.0;
}
