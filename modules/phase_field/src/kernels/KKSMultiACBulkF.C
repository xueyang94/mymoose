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
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "ci_names",
      "Phase concentrations in the order of s1, s2, s3, l1, l2, l3, etc. The order of the "
      "component species must match that in the global_c.");
  params.addRequiredCoupledVar("etas", "Order parameters for all phases.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "dcidc_names",
      "The names of dci/dc in the order of ds1ds, ds2ds, ds3ds, dl1dl, dl2dl, dl3dl, etc. The "
      "order of compoment species must match that of the global_c.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "dcidetaj_names",
      "The names of dci/detaj in the order of dc1deta1, dc2deta1, dc3deta1, dc1deta2, dc2deta2, "
      "dc3deta2, etc. It must match the order of ci_names and etas.");
  params.addRequiredParam<Real>("wi", "Double well height parameter.");
  params.addRequiredParam<MaterialPropertyName>(
      "gi_name", "Base name for the double well function g_i(eta_i) for the given phase");
  params.addRequiredCoupledVar("global_c", "Global concentrations s, l, etc.");
  return params;
}

KKSMultiACBulkF::KKSMultiACBulkF(const InputParameters & parameters)
  : KKSMultiACBulkBase(parameters),
    _ci_names(getParam<std::vector<MaterialPropertyName>>("ci_names")),
    _eta_names(coupledComponents("etas")),
    _eta_map(getParameterJvarMap("etas")),
    _k(-1),
    _prop_d2hjdetapdetai(_num_j),
    _dcidc_names(getParam<std::vector<MaterialPropertyName>>("dcidc_names")),
    _prop_dcidc(_num_j),

    _dcidetaj_names(getParam<std::vector<MaterialPropertyName>>("dcidetaj_names")),
    _prop_dcidetaj(_num_j),
    _prop_dFidci(_num_j),
    _wi(getParam<Real>("wi")),
    _gi_name(getParam<MaterialPropertyName>("gi_name")),
    _prop_dgi(getMaterialPropertyDerivative<Real>("gi_name", _etai_name)),
    _prop_d2gi(getMaterialPropertyDerivative<Real>("gi_name", _etai_name, _etai_name)),

    _c_names(coupledComponents("global_c")),
    _c_map(getParameterJvarMap("global_c")),
    _num_c(_c_names.size())
{

  for (unsigned int i = 0; i < _num_j; ++i)
  {
    // get order parameter names and variable indices
    _eta_names[i] = getVar("etas", i)->name();

    // Set _k to the position of the nonlinear variable in the list of etaj's
    if (coupled("etas", i) == _var.number())
      _k = i;

    _prop_d2hjdetapdetai[i].resize(_num_j);
    _prop_dcidetaj[i].resize(_num_j);
    _prop_dFidci[i].resize(_num_c);
    _prop_dcidc[i].resize(_num_c);
  }

  for (unsigned int m = 0; m < _num_j; ++m) // loop through phases
  {
    for (unsigned int n = 0; n < _num_c; ++n) // loop through components
    {
      _prop_dFidci[m][n] =
          &getMaterialPropertyDerivative<Real>(_Fj_names[m], _ci_names[n * _num_j + m]);

      _prop_dcidc[m][n] = &getMaterialPropertyByName<Real>(_dcidc_names[n * _num_j + m]);
    }

    // Get the derivative of dhjdetap wrt all order parameters
    for (unsigned int n = 0; n < _num_j; ++n) // loop through phases
      _prop_d2hjdetapdetai[m][n] =
          &getMaterialPropertyDerivative<Real>(_hj_names[m], _eta_names[_k], _eta_names[n]);
  }

  // Get dcidetaj indexes by converting the vector of _dcidetaj_names to the matrix of
  // _prop_dcidetaj, so that _prop_dcidetaj[m][n] is dci[m]/detaj[n]
  for (unsigned int m = 0; m < _num_j; ++m)
  {
    for (unsigned int n = 0; n < _num_j; ++n)
      _prop_dcidetaj[m][n] = &getMaterialPropertyByName<Real>(_dcidetaj_names[m * _num_j + n]);
  }
}

Real
KKSMultiACBulkF::computeDFDOP(PFFunctionType type)
{
  Real sum = 0.0;

  switch (type)
  {
    case Residual:
      for (unsigned int n = 0; n < _num_j; ++n)
        sum += (*_prop_dhjdetai[n])[_qp] * (*_prop_Fj[n])[_qp];

      return sum + _wi * _prop_dgi[_qp];

    case Jacobian:
      // For when this kernel is used in the Lagrange multiplier equation
      // In that case the Lagrange multiplier is the nonlinear variable
      if (_etai_var != _var.number())
        return 0.0;

      // For when eta_i is the nonlinear variable
      for (unsigned int n = 0; n < _num_j; ++n)
        sum +=
            (*_prop_d2hjdetapdetai[n][_k])[_qp] * (*_prop_Fj[n])[_qp] +
            (*_prop_dhjdetai[n])[_qp] * (*_prop_dFidci[n][0])[_qp] * (*_prop_dcidetaj[n][_k])[_qp];

      return _phi[_j][_qp] * (sum + _wi * _prop_d2gi[_qp]);
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
    for (unsigned int n = 0; n < _num_j; ++n)
      sum += (*_prop_dhjdetai[n])[_qp] * (*_prop_dFidci[n][compvar])[_qp] *
             (*_prop_dcidc[n][compvar])[_qp];

    res += _L[_qp] * sum * _phi[_j][_qp] * _test[_i][_qp];

    return res;
  }

  // if order parameters are the coupled variables
  auto etavar = mapJvarToCvar(jvar, _eta_map);
  if (etavar >= 0)
  {
    for (unsigned int n = 0; n < _num_j; ++n)
    {
      sum += (*_prop_d2hjdetapdetai[n][etavar])[_qp] * (*_prop_Fj[n])[_qp] +
             (*_prop_dhjdetai[n])[_qp] * (*_prop_dFidci[n][0])[_qp] *
                 (*_prop_dcidetaj[n][etavar])[_qp];
    }

    res += _L[_qp] * sum * _phi[_j][_qp] * _test[_i][_qp];

    return res;
  }

  // Handle the case when this kernel is used in the Lagrange multiplier equation
  // In this case the second derivative of the barrier function contributes
  // to the off-diagonal Jacobian
  if (jvar == _etai_var)
    sum += _wi * _prop_d2gi[_qp];

  // get the coupled variable jvar is referring to
  const unsigned int cvar = mapJvarToCvar(jvar);

  // add dependence of KKSMultiACBulkF on other variables
  for (unsigned int n = 0; n < _num_j; ++n)
    sum += (*_prop_d2hjdetaidarg[n][cvar])[_qp] * (*_prop_Fj[n])[_qp] +
           (*_prop_dhjdetai[n])[_qp] * (*_prop_dFjdarg[n][cvar])[_qp];

  res += _L[_qp] * sum * _phi[_j][_qp] * _test[_i][_qp];

  return res;
}
