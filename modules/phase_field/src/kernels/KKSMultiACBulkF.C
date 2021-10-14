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
  params.addRequiredParam<std::vector<MaterialPropertyName>>("ci_names", "Phase concentrations");
  params.addCoupledVar("etas", "Order parameters for all phases.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>("dcidc_names", "The names of dci/dc");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "dcidetaj_names",
      "The names of dci/detaj in the order of dc1deta1, dc2deta1, dc3deta1, dc1deta2, dc2deta2, "
      "dc3deta2, etc");
  params.addRequiredParam<Real>("wi", "Double well height parameter");
  params.addRequiredParam<MaterialPropertyName>(
      "gi_name", "Base name for the double well function g_i(eta_i) for the given phase");
  params.addCoupledVar("global_c", "Global concentration.");
  return params;
}

KKSMultiACBulkF::KKSMultiACBulkF(const InputParameters & parameters)
  : KKSMultiACBulkBase(parameters),
    _ci_names(getParam<std::vector<MaterialPropertyName>>("ci_names")),
    _eta_names(coupledComponents("etas")),
    _eta_map(getParameterJvarMap("etas")),
    _k(-1),
    _prop_dhjdetap(_num_j),
    _prop_d2hjdetapdetai(_num_j),
    _dcidc_names(getParam<std::vector<MaterialPropertyName>>("dcidc_names")),
    _prop_dcidc(_num_j),

    _dcidetaj_names(getParam<std::vector<MaterialPropertyName>>("dcidetaj_names")),
    _prop_dcidetaj(_num_j),
    _prop_dFidci(_num_j),
    _wi(getParam<Real>("wi")),
    _gi_name(getParam<MaterialPropertyName>("gi_name")),
    _prop_dgi(getMaterialPropertyDerivative<Real>("gi_name", _etai_name)),
    _prop_d2gpdetapdetai(_num_j),
    _c_var(coupled("global_c"))
{

  for (unsigned int i = 0; i < _num_j; ++i)
  {
    // get order parameter names and variable indices
    _eta_names[i] = getVar("etas", i)->name();

    // Set _k to the position of the nonlinear variable in the list of etaj's
    if (coupled("etas", i) == _var.number())
      _k = i;

    // declare dcidc material properties
    _prop_dcidc[i] = &getMaterialPropertyByName<Real>(_dcidc_names[i]);
  }

  for (unsigned int m = 0; m < _num_j; ++m)
  {
    // Get dFidci
    _prop_dFidci[m] = &getMaterialPropertyDerivative<Real>(_Fj_names[m], _ci_names[m]);

    // Get the derivatives of switching functions wrt phase_eta
    _prop_dhjdetap[m] = &getMaterialPropertyDerivative<Real>(_hj_names[m], _eta_names[_k]);

    // Get the derivative of dgidetap wrt all order parameters
    _prop_d2gpdetapdetai[m] =
        &getMaterialPropertyDerivative<Real>(_gi_name, _eta_names[_k], _eta_names[m]);

    _prop_d2hjdetapdetai[m].resize(_num_j);
    _prop_dcidetaj[m].resize(_num_j);

    for (unsigned int n = 0; n < _num_j; ++n)
    {
      // Get the derivative of dhjdetap wrt all order parameters
      _prop_d2hjdetapdetai[m][n] =
          &getMaterialPropertyDerivative<Real>(_hj_names[m], _eta_names[_k], _eta_names[n]);
    }
  }

  // Get dcidetaj indexes by converting the vector of _dcidetaj_names to the matrix of
  // _prop_dcidetaj, so that _prop_dcidetaj[m][n] is dci[m]/detaj[n]
  for (unsigned int i = 0; i < _num_j * _num_j; ++i)
  {
    if (i >= 0 && i < _num_j)
    {
      _prop_dcidetaj[i][0] = &getMaterialPropertyByName<Real>(_dcidetaj_names[i]);
      continue;
    }
    if (i >= _num_j && i < 2 * _num_j)
    {
      _prop_dcidetaj[i - _num_j][1] = &getMaterialPropertyByName<Real>(_dcidetaj_names[i]);
      continue;
    }
    if (i >= 2 * _num_j && i < _num_j * _num_j)
    {
      _prop_dcidetaj[i - 2 * _num_j][2] = &getMaterialPropertyByName<Real>(_dcidetaj_names[i]);
      continue;
    }
  }
}

Real
KKSMultiACBulkF::computeDFDOP(PFFunctionType type)
{
  switch (type)
  {
    case Residual:
      return (*_prop_dhjdetap[0])[_qp] * (*_prop_Fj[0])[_qp] +
             (*_prop_dhjdetap[1])[_qp] * (*_prop_Fj[1])[_qp] +
             (*_prop_dhjdetap[2])[_qp] * (*_prop_Fj[2])[_qp] + _wi * _prop_dgi[_qp];

    case Jacobian:
      return ((*_prop_d2hjdetapdetai[0][_k])[_qp] * (*_prop_Fj[0])[_qp] +
              (*_prop_dhjdetap[0])[_qp] * (*_prop_dFidci[0])[_qp] * (*_prop_dcidetaj[0][_k])[_qp] +
              (*_prop_d2hjdetapdetai[1][_k])[_qp] * (*_prop_Fj[1])[_qp] +
              (*_prop_dhjdetap[1])[_qp] * (*_prop_dFidci[1])[_qp] * (*_prop_dcidetaj[1][_k])[_qp] +
              (*_prop_d2hjdetapdetai[2][_k])[_qp] * (*_prop_Fj[2])[_qp] +
              (*_prop_dhjdetap[2])[_qp] * (*_prop_dFidci[2])[_qp] * (*_prop_dcidetaj[2][_k])[_qp] +
              _wi * (*_prop_d2gpdetapdetai[_k])[_qp]) *
             _phi[_j][_qp];
  }
  mooseError("Invalid type passed in");
}

Real
KKSMultiACBulkF::computeQpOffDiagJacobian(unsigned int jvar)
{
  // first get dependence of mobility _L on other variables using parent class member function Real
  Real res = ACBulk<Real>::computeQpOffDiagJacobian(jvar);

  Real sum = 0.0;

  // if c is the coupled variable
  if (jvar == _c_var)
  {
    for (unsigned int n = 0; n < _num_j; ++n)
      sum += (*_prop_dhjdetap[n])[_qp] * (*_prop_dFidci[n])[_qp] * (*_prop_dcidc[n])[_qp];

    res += _L[_qp] * sum * _phi[_j][_qp] * _test[_i][_qp];

    return res;
  }

  // if order parameters are the coupled variables
  auto etavar = mapJvarToCvar(jvar, _eta_map);
  if (etavar >= 0)
  {
    for (unsigned int n = 0; n < _num_j; ++n)
    {
      sum +=
          (*_prop_d2hjdetapdetai[n][etavar])[_qp] * (*_prop_Fj[n])[_qp] +
          (*_prop_dhjdetap[n])[_qp] * (*_prop_dFidci[n])[_qp] * (*_prop_dcidetaj[n][etavar])[_qp];
    }

    sum += _wi * (*_prop_d2gpdetapdetai[etavar])[_qp];
  }

  res += _L[_qp] * sum * _phi[_j][_qp] * _test[_i][_qp];

  return res;
}
