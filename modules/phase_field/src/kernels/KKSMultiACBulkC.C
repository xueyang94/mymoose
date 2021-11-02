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
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "ci_names", "Phase concentrations. These must have the same order as Fj_names.");
  params.addRequiredCoupledVar("etas", "Order parameters for all phases.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>("dcidc_names", "The names of dci/dc");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "dcidetaj_names",
      "The names of dci/detaj in the order of dc1deta1, dc2deta1, dc3deta1, dc1deta2, dc2deta2, "
      "dc3deta2, etc");
  params.addParam<std::vector<MaterialPropertyName>>(
      "coupled_dcidb_names",
      "Coupled dcidb in the order of dc1db, dc2db, dc3db, etc. These must have the same order as "
      "Fj_names");
  params.addRequiredCoupledVar(
      "global_c", "The global concentration of the component corresponding to ci_names.");
  params.addRequiredCoupledVar("other_global_c", "The other coupled global concentrations.");
  return params;
}

KKSMultiACBulkC::KKSMultiACBulkC(const InputParameters & parameters)
  : KKSMultiACBulkBase(parameters),
    _ci_names(getParam<std::vector<MaterialPropertyName>>("ci_names")),
    _prop_ci(_num_j),
    _eta_names(coupledComponents("etas")),
    _eta_map(getParameterJvarMap("etas")),
    _k(-1),

    _prop_d2hjdetapdetai(_num_j),
    _dcidc_names(getParam<std::vector<MaterialPropertyName>>("dcidc_names")),
    _prop_dcidc(_num_j),

    _dcidetaj_names(getParam<std::vector<MaterialPropertyName>>("dcidetaj_names")),
    _prop_dcidetaj(_num_j),

    _coupled_dcidb_names(getParam<std::vector<MaterialPropertyName>>("coupled_dcidb_names")),
    _prop_coupled_dcidb(_num_j),

    _first_df1(getMaterialPropertyDerivative<Real>(_Fj_names[0], _ci_names[0])),
    _second_df1(getMaterialPropertyDerivative<Real>(_Fj_names[0], _ci_names[0], _ci_names[0])),
    _c_var(coupled("global_c")),
    _b_var(coupled("other_global_c"))
{
  // initialize coupled dcidb
  for (unsigned int i = 0; i < _num_j; ++i)
    _prop_coupled_dcidb[i] = &getMaterialPropertyByName<Real>(_coupled_dcidb_names[i]);

  // get ci values
  for (unsigned int i = 0; i < _num_j; ++i)
    _prop_ci[i] = &getMaterialPropertyByName<Real>(_ci_names[i]);

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
  // _prop_dcidetaj where _prop_dcidetaj[m][n] is dci[m]/detaj[n]
  for (unsigned int m = 0; m < _num_j; ++m)
  {
    for (unsigned int n = 0; n < _num_j; ++n)
      _prop_dcidetaj[m][n] = &getMaterialPropertyByName<Real>(_dcidetaj_names[m * _num_j + n]);
  }
}

Real
KKSMultiACBulkC::computeDFDOP(PFFunctionType type)
{
  Real sum = 0.0;
  Real sum1 = 0.0;
  Real sum2 = 0.0;

  switch (type)
  {
    case Residual:
      for (unsigned int n = 0; n < _num_j; ++n)
        sum += (*_prop_dhjdetai[n])[_qp] * (*_prop_ci[n])[_qp];

      return -_first_df1[_qp] * sum;

    case Jacobian:
      // For when this kernel is used in the Lagrange multiplier equation
      // In that case the Lagrange multiplier is the nonlinear variable
      if (_etai_var != _var.number())
        return 0.0;

      for (unsigned int n = 0; n < _num_j; ++n)
      {
        sum1 += (*_prop_dhjdetai[n])[_qp] * (*_prop_ci[n])[_qp];
        sum2 += (*_prop_d2hjdetapdetai[n][_k])[_qp] * (*_prop_ci[n])[_qp] +
                (*_prop_dhjdetai[n])[_qp] * (*_prop_dcidetaj[n][_k])[_qp];
      }

      return -(_second_df1[_qp] * (*_prop_dcidetaj[0][_k])[_qp] * sum1 + _first_df1[_qp] * sum2) *
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
  Real sum1 = 0.0;
  Real sum2 = 0.0;

  // if c is the coupled variable
  if (jvar == _c_var)
  {
    for (unsigned int n = 0; n < _num_j; ++n)
    {
      sum1 += (*_prop_dhjdetai[n])[_qp] * (*_prop_ci[n])[_qp];
      sum2 += (*_prop_dhjdetai[n])[_qp] * (*_prop_dcidc[n])[_qp];
    }

    sum = _second_df1[_qp] * (*_prop_dcidc[0])[_qp] * sum1 + _first_df1[_qp] * sum2;

    res += -_L[_qp] * sum * _phi[_j][_qp] * _test[_i][_qp];

    return res;
  }

  // if b is the coupled variable
  if (jvar == _b_var)
  {
    for (unsigned int n = 0; n < _num_j; ++n)
    {
      sum1 += (*_prop_dhjdetai[n])[_qp] * (*_prop_ci[n])[_qp];
      sum2 += (*_prop_dhjdetai[n])[_qp] * (*_prop_coupled_dcidb[n])[_qp];
    }

    res += -_L[_qp] *
           (_second_df1[_qp] * (*_prop_coupled_dcidb[0])[_qp] * sum1 + _first_df1[_qp] * sum2) *
           _phi[_j][_qp] * _test[_i][_qp];

    return res;
  }

  // if order parameters are the coupled variables
  auto etavar = mapJvarToCvar(jvar, _eta_map);
  if (etavar >= 0)
  {
    for (unsigned int n = 0; n < _num_j; ++n)
    {
      sum1 += (*_prop_dhjdetai[n])[_qp] * (*_prop_ci[n])[_qp];
      sum2 += (*_prop_d2hjdetapdetai[n][etavar])[_qp] * (*_prop_ci[n])[_qp] +
              (*_prop_dhjdetai[n])[_qp] * (*_prop_dcidetaj[n][etavar])[_qp];
    }
    sum = _second_df1[_qp] * (*_prop_dcidetaj[0][etavar])[_qp] * sum1 + _first_df1[_qp] * sum2;
  }

  res += -_L[_qp] * sum * _phi[_j][_qp] * _test[_i][_qp];

  return res;
}
