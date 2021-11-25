//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KKSMultiACBulkF.h"

registerMooseObject("PhaseFieldApp", KKSMultiACBulkF);

InputParameters
KKSMultiACBulkF::validParams()
{
  InputParameters params = KKSMultiACBulkBase::validParams();
  params.addClassDescription("KKS model kernel (part 1 of 2) for the Bulk Allen-Cahn. This "
                             "includes all terms NOT dependent on chemical potential.");
  params.addRequiredCoupledVar("global_cs", "Global concentrations c, b, etc.");
  params.addRequiredCoupledVar("etas", "Order parameters for all phases.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "c1_names",
      "Phase concentrations in the frist phase of etas. The order must match global_cs, for "
      "example, c1, b1, etc.");
  params.addParam<std::vector<MaterialPropertyName>>(
      "dcidb_names",
      "The phase concentrations taken derivatives wrt global concentrations. i must match the "
      "order of etas. ci and b must match the order of global_cs. First keep the same b and loop "
      "through ci for one species, for example, dc1dc, dc2dc, dc3dc, "
      "dc1db, dc2db, dc3db, db1dc, db2dc, db3dc, db1db, db2db, db3db, etc.");
  params.addParam<std::vector<MaterialPropertyName>>(
      "dcidetaj_names",
      "The phase concentrations taken derivatives wrt kernel variable. ci must match the order in "
      "global_cs and etas, and etaj must match the order in etas, for example, dc1deta1, dc2deta1, "
      "dc3deta1, dc1deta2...dc1deta3...db1deta1...db2deta1...db3deta1..., etc.");
  params.addRequiredParam<Real>("wi", "Double well height parameter.");
  params.addRequiredParam<MaterialPropertyName>(
      "gi_name", "Base name for the double well function g_i(eta_i) for the given phase");
  return params;
}

KKSMultiACBulkF::KKSMultiACBulkF(const InputParameters & parameters)
  : KKSMultiACBulkBase(parameters),
    _c_names(coupledComponents("global_cs")),
    _c_map(getParameterJvarMap("global_cs")),
    _num_c(coupledComponents("global_cs")),
    _eta_names(coupledComponents("etas")),
    _eta_map(getParameterJvarMap("etas")),
    _k(-1),
    _c1_names(getParam<std::vector<MaterialPropertyName>>("c1_names")),

    _dcidb_names(getParam<std::vector<MaterialPropertyName>>("dcidb_names")),
    _prop_dcidb(_num_c),

    _dcidetaj_names(getParam<std::vector<MaterialPropertyName>>("dcidetaj_names")),
    _prop_dcidetaj(_num_c),

    _wi(getParam<Real>("wi")),
    _gi_name(getParam<MaterialPropertyName>("gi_name")),
    _dgi(getMaterialPropertyDerivative<Real>("gi_name", _etai_name)),
    _d2gi(getMaterialPropertyDerivative<Real>("gi_name", _etai_name, _etai_name)),
    _d2hjdetaidetap(_num_j),
    _dF1dc1(_num_c)
{
  for (unsigned int i = 0; i < _num_j; ++i)
  {
    // get order parameter names and variable indices
    _eta_names[i] = getVar("etas", i)->name();

    // Set _k to the position of the nonlinear variable in the list of etaj's
    if (coupled("etas", i) == _var.number())
      _k = i;
  }

  // initialize _prop_dcidb
  for (unsigned int m = 0; m < _num_c; ++m)
  {
    _prop_dcidb[m].resize(_num_j);

    for (unsigned int n = 0; n < _num_j; ++n)
    {
      _prop_dcidb[m][n].resize(_num_c);

      for (unsigned int l = 0; l < _num_c; ++l)
        _prop_dcidb[m][n][l] =
            &getMaterialPropertyByName<Real>(_dcidb_names[m * _num_j * _num_c + n + l * _num_j]);
    }
  }

  // initialize _prop_dcidetaj
  for (unsigned int m = 0; m < _num_c; ++m)
  {
    _prop_dcidetaj[m].resize(_num_j);

    for (unsigned int n = 0; n < _num_j; ++n)
    {
      _prop_dcidetaj[m][n].resize(_num_j);

      for (unsigned int l = 0; l < _num_j; ++l)
        _prop_dcidetaj[m][n][l] =
            &getMaterialPropertyByName<Real>(_dcidetaj_names[m * _num_j * _num_j + n + l * _num_j]);
    }
  }

  // initialize _dF1dc1 as dF1dc1, dF1db1, etc
  for (unsigned int i = 0; i < _num_c; ++i)
    _dF1dc1[i] = &getMaterialPropertyDerivative<Real>(_Fj_names[0], _c1_names[i]);

  // initialize _d2hjdetaidetap
  for (unsigned int m = 0; m < _num_j; ++m)
  {
    _d2hjdetaidetap[m].resize(_num_j);

    // Get the derivative of dhjdetai wrt all order parameters p
    for (unsigned int n = 0; n < _num_j; ++n)
      _d2hjdetaidetap[m][n] =
          &getMaterialPropertyDerivative<Real>(_hj_names[m], _eta_names[_k], _eta_names[n]);
  }
}

Real
KKSMultiACBulkF::computeDFDOP(PFFunctionType type)
{
  Real sum = 0.0;

  switch (type)
  {
    case Residual:
      for (unsigned int m = 0; m < _num_j; ++m)
        sum += (*_prop_dhjdetai[m])[_qp] * (*_prop_Fj[m])[_qp];

      return sum + _wi * _dgi[_qp];

    case Jacobian:
      // For when this kernel is used in the Lagrange multiplier equation
      // In that case the Lagrange multiplier is the nonlinear variable
      if (_etai_var != _var.number())
        return 0.0;

      // For when eta_i is the nonlinear variable
      for (unsigned int m = 0; m < _num_j; ++m)
      {
        Real sum1 = 0.0;

        for (unsigned int n = 0; n < _num_c; ++n)
          sum1 += (*_dF1dc1[n])[_qp] * (*_prop_dcidetaj[n][m][_k])[_qp];

        sum +=
            (*_d2hjdetaidetap[m][_k])[_qp] * (*_prop_Fj[m])[_qp] + (*_prop_dhjdetai[m])[_qp] * sum1;
      }

      return _phi[_j][_qp] * (sum + _wi * _d2gi[_qp]);
  }

  mooseError("Invalid type passed in");
}

Real
KKSMultiACBulkF::computeQpOffDiagJacobian(unsigned int jvar)
{
  // first get dependence of mobility _L on other variables using parent class member function Real
  Real res = ACBulk<Real>::computeQpOffDiagJacobian(jvar);

  Real sum = 0.0;

  // if concentrations are the coupled variables
  auto compvar = mapJvarToCvar(jvar, _c_map);
  if (compvar >= 0)
  {
    for (unsigned int m = 0; m < _num_j; ++m)
    {
      Real sum1 = 0.0;

      for (unsigned int n = 0; n < _num_c; ++n)
        sum1 += (*_dF1dc1[n])[_qp] * (*_prop_dcidb[n][m][compvar])[_qp];

      sum += (*_prop_dhjdetai[m])[_qp] * sum1;
    }

    res += _L[_qp] * sum * _phi[_j][_qp] * _test[_i][_qp];

    return res;
  }

  // if order parameters are the coupled variables
  auto etavar = mapJvarToCvar(jvar, _eta_map);
  if (etavar >= 0)
  {
    for (unsigned int m = 0; m < _num_j; ++m)
    {
      Real sum1 = 0.0;

      for (unsigned int n = 0; n < _num_c; ++n)
        sum1 += (*_dF1dc1[n])[_qp] * (*_prop_dcidetaj[n][m][etavar])[_qp];

      sum += (*_d2hjdetaidetap[m][etavar])[_qp] * (*_prop_Fj[m])[_qp] +
             (*_prop_dhjdetai[m])[_qp] * sum1;
    }

    res += _L[_qp] * sum * _phi[_j][_qp] * _test[_i][_qp];

    return res;
  }

  // // get the coupled variable jvar is referring to
  // const unsigned int cvar = mapJvarToCvar(jvar);
  // // add dependence of KKSMultiACBulkF on other variables
  // for (unsigned int n = 0; n < _num_j; ++n)
  //   sum += (*_prop_d2hjdetaidarg[n][cvar])[_qp] * (*_prop_Fj[n])[_qp] +
  //          (*_prop_dhjdetai[n])[_qp] * (*_prop_dFjdarg[n][cvar])[_qp];

  // Handle the case when this kernel is used in the Lagrange multiplier equation
  // In this case the second derivative of the barrier function contributes
  // to the off-diagonal Jacobian
  if (jvar == _etai_var)
  {
    sum += _wi * _d2gi[_qp];

    res += _L[_qp] * sum * _phi[_j][_qp] * _test[_i][_qp];
  }

  return res;
}
