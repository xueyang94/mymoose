//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KKSMultiACBulkC.h"

registerMooseObject("PhaseFieldApp", KKSMultiACBulkC);

InputParameters
KKSMultiACBulkC::validParams()
{
  InputParameters params = KKSMultiACBulkBase::validParams();
  params.addClassDescription("Multi-phase KKS model kernel (part 2 of 2) for the Bulk Allen-Cahn. "
                             "This includes all terms dependent on chemical potential.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>("ci_names", "Phase concentrations");
  params.addCoupledVar("etas", "Order parameters for all phases.");
  params.addRequiredParam<MaterialPropertyName>("dc1dc_name", "The name of dc1/dc");
  params.addRequiredParam<MaterialPropertyName>("dc2dc_name", "The name of dc2/dc");
  params.addRequiredParam<MaterialPropertyName>("dc3dc_name", "The name of dc3/dc");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "dcidetaj_names",
      "The name of dci/detaj in the order of dc1deta1, dc2deta1, dc3deta1, dc1deta2, dc2deta2, "
      "dc3deta2, etc");
  params.addRequiredParam<MaterialPropertyName>("F1_name",
                                                "The name of the bulk energy of phase 1");
  params.addCoupledVar("global_c", "Global concentration.");
  return params;
}

KKSMultiACBulkC::KKSMultiACBulkC(const InputParameters & parameters)
  : KKSMultiACBulkBase(parameters),
    _ci_names(getParam<std::vector<MaterialPropertyName>>("ci_names")),
    _prop_ci(_num_j),
    _eta_names(coupledComponents("etas")),
    _eta_map(getParameterJvarMap("etas")),
    _k(-1),
    _prop_dhjdetap(_num_j),
    _prop_d2hjdetapdetai(_num_j),
    _dc1dc(getMaterialProperty<Real>("dc1dc_name")),
    _dc2dc(getMaterialProperty<Real>("dc2dc_name")),
    _dc3dc(getMaterialProperty<Real>("dc3dc_name")),
    _dcidetaj_names(getParam<std::vector<MaterialPropertyName>>("dcidetaj_names")),
    _prop_dcidetaj(_num_j),
    _c1_name("c1"),
    _first_df1(getMaterialPropertyDerivative<Real>("F1_name", _c1_name)),
    _second_df1(getMaterialPropertyDerivative<Real>("F1_name", _c1_name, _c1_name)),
    _c_var(coupled("global_c"))
{
  for (unsigned int i = 0; i < _num_j; ++i)
    _prop_ci[i] = &getMaterialPropertyByName<Real>(_ci_names[i]);

  // get order parameter names and variable indices
  for (unsigned int i = 0; i < _num_j; ++i)
  {
    _eta_names[i] = getVar("etas", i)->name();
    // Set _k to the position of the nonlinear variable in the list of etaj's
    if (coupled("etas", i) == _var.number())
      _k = i;
  }

  for (unsigned int m = 0; m < _num_j; ++m)
  {
    // Get the derivatives of switching functions wrt phase_eta
    _prop_dhjdetap[m] = &getMaterialPropertyDerivative<Real>(_hj_names[m], _eta_names[_k]);

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
    }
    if (i >= _num_j && i < 2 * _num_j)
    {
      _prop_dcidetaj[i - _num_j][1] = &getMaterialPropertyByName<Real>(_dcidetaj_names[i]);
    }
    if (i >= 2 * _num_j && i < _num_j * _num_j)
    {
      _prop_dcidetaj[i - 2 * _num_j][2] = &getMaterialPropertyByName<Real>(_dcidetaj_names[i]);
    }
  }
}

Real
KKSMultiACBulkC::computeDFDOP(PFFunctionType type)
{
  switch (type)
  {
    case Residual:
      return _first_df1[_qp] * ((*_prop_dhjdetap[0])[_qp] * (*_prop_ci[0])[_qp] +
                                (*_prop_dhjdetap[1])[_qp] * (*_prop_ci[1])[_qp] +
                                (*_prop_dhjdetap[2])[_qp] * (*_prop_ci[2])[_qp]);

    case Jacobian:
      return _first_df1[_qp] *
             ((*_prop_d2hjdetapdetai[0][0])[_qp] * (*_prop_ci[0])[_qp] +
              (*_prop_dhjdetap[0])[_qp] * (*_prop_dcidetaj[0][_k])[_qp] +
              (*_prop_d2hjdetapdetai[1][0])[_qp] * (*_prop_ci[1])[_qp] +
              (*_prop_dhjdetap[1])[_qp] * (*_prop_dcidetaj[1][_k])[_qp] +
              (*_prop_d2hjdetapdetai[2][0])[_qp] * (*_prop_ci[2])[_qp] +
              (*_prop_dhjdetap[2])[_qp] * (*_prop_dcidetaj[2][_k])[_qp]) *
             _phi[_j][_qp];
  }
  mooseError("Invalid type passed in");
}

Real
KKSMultiACBulkC::computeQpOffDiagJacobian(unsigned int jvar)
{
  // first get dependence of mobility _L on other variables using parent class member function Real
  Real res = ACBulk<Real>::computeQpOffDiagJacobian(jvar);

  Real sum = 0.0;

  // if c is the coupled variable
  if (jvar == _c_var)
    sum = _second_df1[_qp] * _dc1dc[_qp] *
              ((*_prop_dhjdetap[0])[_qp] * (*_prop_ci[0])[_qp] +
               (*_prop_dhjdetap[1])[_qp] * (*_prop_ci[1])[_qp] +
               (*_prop_dhjdetap[2])[_qp] * (*_prop_ci[2])[_qp]) +
          _first_df1[_qp] *
              ((*_prop_dhjdetap[0])[_qp] * _dc1dc[_qp] + (*_prop_dhjdetap[1])[_qp] * _dc2dc[_qp] +
               (*_prop_dhjdetap[2])[_qp] * _dc3dc[_qp]);

  // if other order parameters are the coupled variables
  auto etavar = mapJvarToCvar(jvar, _eta_map);
  if (etavar >= 0)
  {
    for (unsigned int n = 0; n < _num_j; ++n)
      sum += _first_df1[_qp] * ((*_prop_d2hjdetapdetai[n][etavar])[_qp] * (*_prop_ci[n])[_qp] +
                                (*_prop_dhjdetap[n])[_qp] * (*_prop_dcidetaj[n][etavar])[_qp]);
  }

  res += -_L[_qp] * sum * _phi[_j][_qp] * _test[_i][_qp];

  return res;
}
