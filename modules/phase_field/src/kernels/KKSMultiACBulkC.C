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
  params.addRequiredCoupledVar(
      "etas", "Order parameters for all phases. Place in the same order as Fj_names.");
  params.addRequiredCoupledVar(
      "global_c", "The global concentration of the component corresponding to ci_names.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "ci_names", "Phase concentrations. These must have the same order as Fj_names.");
  params.addRequiredCoupledVar("other_global_b", "The other coupled global concentrations.");
  params.addParam<std::vector<MaterialPropertyName>>(
      "other_b1_names",
      "The phase concentration of other_global_b in the first phase of Fj_names. The order "
      "must match other_global_b.");

  params.addRequiredParam<std::vector<MaterialPropertyName>>("dcidc_names", "The names of dci/dc");
  params.addParam<std::vector<MaterialPropertyName>>(
      "dcidb_names",
      "Coupled dcidb in the order of dc1db, dc2db, dc3db, etc. These must have the same order as "
      "Fj_names");
  params.addParam<std::vector<MaterialPropertyName>>(
      "db1dc_names",
      "The phase concentration of other_global_b in the frist phase of Fj_names taken derivarive"
      "wrt global_c. The b1 order must match other_global_b.");
  params.addParam<std::vector<MaterialPropertyName>>(
      "db1db_names",
      "The phase concentration of other_global_b in the first phase of Fj_names taken "
      "derivative wrt other_global_b. The order must match other_global_b.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "dcidetaj_names",
      "The names of dci/detaj in the order of dc1deta1, dc2deta1, dc3deta1, dc1deta2, dc2deta2, "
      "dc3deta2, etc");
  params.addParam<std::vector<MaterialPropertyName>>(
      "db1detaj_names",
      "The phase concentration of other_global_b in the first phase of Fj_names taken derivative"
      "wrt eta in the order of db1deta1, db1deta2, db1deta3, de1deta1, de1deta2, de1deta3, etc. b1 "
      "must have the order of other_global_b and eta must match Fj_names.");
  return params;
}

KKSMultiACBulkC::KKSMultiACBulkC(const InputParameters & parameters)
  : KKSMultiACBulkBase(parameters),
    _eta_names(coupledComponents("etas")),
    _eta_map(getParameterJvarMap("etas")),
    _k(-1),
    _c_var(coupled("global_c")),
    _ci_names(getParam<std::vector<MaterialPropertyName>>("ci_names")),
    _prop_ci(_num_j),
    _num_other_b(coupledComponents("other_global_b")),
    _other_b_map(getParameterJvarMap("other_global_b")),
    _other_b1_names(getParam<std::vector<MaterialPropertyName>>("other_b1_names")),

    _dcidc_names(getParam<std::vector<MaterialPropertyName>>("dcidc_names")),
    _prop_dcidc(_num_j),
    _dcidb_names(getParam<std::vector<MaterialPropertyName>>("dcidb_names")),
    _prop_dcidb(_num_j),
    _db1dc_names(getParam<std::vector<MaterialPropertyName>>("db1dc_names")),
    _prop_db1dc(_num_other_b),
    _db1db_names(getParam<std::vector<MaterialPropertyName>>("db1db_names")),
    _prop_db1db(_num_other_b),
    _dcidetaj_names(getParam<std::vector<MaterialPropertyName>>("dcidetaj_names")),
    _prop_dcidetaj(_num_j),
    _db1detaj_names(getParam<std::vector<MaterialPropertyName>>("db1detaj_names")),
    _prop_db1detaj(_num_other_b),

    _prop_d2hjdetaidetap(_num_j),
    _first_df1(getMaterialPropertyDerivative<Real>(_Fj_names[0], _ci_names[0])),
    _second_df1(getMaterialPropertyDerivative<Real>(_Fj_names[0], _ci_names[0], _ci_names[0])),
    _d2F1dc1db1(_num_other_b)
{
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
  }

  // initialize _prop_dcidc
  for (unsigned int i = 0; i < _num_j; ++i)
    _prop_dcidc[i] = &getMaterialPropertyByName<Real>(_dcidc_names[i]);

  // initialize _prop_dcidb[m][n] where m is phase i and n is b index
  for (unsigned int m = 0; m < _num_j; ++m)
  {
    _prop_dcidb[m].resize(_num_other_b);
    for (unsigned int n = 0; n < _num_other_b; ++n)
      _prop_dcidb[m][n] = &getMaterialPropertyByName<Real>(_dcidb_names[m * _num_other_b + n]);
  }

  // initialize _prop_db1dc
  for (unsigned int i = 0; i < _num_other_b; ++i)
    _prop_db1dc[i] = &getMaterialPropertyByName<Real>(_db1dc_names[i]);

  // initialize _prop_db1db
  for (unsigned int i = 0; i < _num_other_b; ++i)
    _prop_db1db[i] = &getMaterialPropertyByName<Real>(_db1db_names[i]);

  // initialize _prop_dcidetaj
  for (unsigned int m = 0; m < _num_j; ++m)
  {
    _prop_dcidetaj[m].resize(_num_j);
    for (unsigned int n = 0; n < _num_j; ++n)
      _prop_dcidetaj[m][n] = &getMaterialPropertyByName<Real>(_dcidetaj_names[m * _num_j + n]);
  }

  // initialize _prop_db1detaj
  for (unsigned int m = 0; m < _num_other_b; ++m)
  {
    _prop_db1detaj[m].resize(_num_j);
    for (unsigned int n = 0; n < _num_j; ++n)
      _prop_db1detaj[m][n] = &getMaterialPropertyByName<Real>(_db1detaj_names[m * _num_j + n]);
  }

  // initialize _prop_d2hjdetaidetap
  for (unsigned int m = 0; m < _num_j; ++m)
  {
    _prop_d2hjdetaidetap[m].resize(_num_j);
    for (unsigned int n = 0; n < _num_j; ++n)
    {
      // Get the derivative of dhjdetai wrt all order parameters
      _prop_d2hjdetaidetap[m][n] =
          &getMaterialPropertyDerivative<Real>(_hj_names[m], _eta_names[_k], _eta_names[n]);
    }
  }

  // initialize _d2F1dc1db1
  for (unsigned int i = 0; i < _num_other_b; ++i)
    _d2F1dc1db1[i] =
        &getMaterialPropertyDerivative<Real>(_Fj_names[0], _ci_names[0], _other_b1_names[i]);
}

Real
KKSMultiACBulkC::computeDFDOP(PFFunctionType type)
{
  Real sum = 0.0;
  Real sum1 = 0.0;
  Real sum2 = 0.0;
  Real factor = 0.0; // the first multiplier term

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

      factor = _second_df1[_qp] * (*_prop_dcidetaj[0][_k])[_qp];

      for (unsigned int n = 0; n < _num_other_b; ++n)
        factor += (*_d2F1dc1db1[n])[_qp] * (*_prop_db1detaj[n][_k])[_qp];

      for (unsigned int n = 0; n < _num_j; ++n)
      {
        sum1 += (*_prop_dhjdetai[n])[_qp] * (*_prop_ci[n])[_qp];
        sum2 += (*_prop_d2hjdetaidetap[n][_k])[_qp] * (*_prop_ci[n])[_qp] +
                (*_prop_dhjdetai[n])[_qp] * (*_prop_dcidetaj[n][_k])[_qp];
      }

      return -(factor * sum1 + _first_df1[_qp] * sum2) * _phi[_j][_qp];
  }
  mooseError("Invalid type passed in");
}

Real
KKSMultiACBulkC::computeQpOffDiagJacobian(unsigned int jvar)
{
  // first get dependence of mobility _L on other variables using parent class member function Real
  Real res = ACBulk<Real>::computeQpOffDiagJacobian(jvar);

  Real sum1 = 0.0;
  Real sum2 = 0.0;
  Real factor = 0.0; // the first multiplier term

  // if c is the coupled variable
  if (jvar == _c_var)
  {
    factor = _second_df1[_qp] * (*_prop_dcidc[0])[_qp];

    for (unsigned int n = 0; n < _num_other_b; ++n)
      factor += (*_d2F1dc1db1[n])[_qp] * (*_prop_db1dc[n])[_qp];

    for (unsigned int n = 0; n < _num_j; ++n)
    {
      sum1 += (*_prop_dhjdetai[n])[_qp] * (*_prop_ci[n])[_qp];
      sum2 += (*_prop_dhjdetai[n])[_qp] * (*_prop_dcidc[n])[_qp];
    }

    res += -_L[_qp] * (factor * sum1 + _first_df1[_qp] * sum2) * _phi[_j][_qp] * _test[_i][_qp];

    return res;
  }

  // if other cs are the coupled variables
  auto other_b_var = mapJvarToCvar(jvar, _other_b_map);
  if (other_b_var >= 0)
  {
    factor = _second_df1[_qp] * (*_prop_dcidb[0][other_b_var])[_qp];

    for (unsigned int m = 0; m < _num_other_b; ++m)
      factor += (*_d2F1dc1db1[other_b_var])[_qp] * (*_prop_db1db[m])[_qp];

    for (unsigned int n = 0; n < _num_j; ++n)
    {
      sum1 += (*_prop_dhjdetai[n])[_qp] * (*_prop_ci[n])[_qp];
      sum2 += (*_prop_dhjdetai[n])[_qp] * (*_prop_dcidb[n][other_b_var])[_qp];
    }

    res += -_L[_qp] * (factor * sum1 + _first_df1[_qp] * sum2) * _phi[_j][_qp] * _test[_i][_qp];

    return res;
  }

  // if order parameters are the coupled variables
  auto etavar = mapJvarToCvar(jvar, _eta_map);
  if (etavar >= 0)
  {
    factor = _second_df1[_qp] * (*_prop_dcidetaj[0][etavar])[_qp];

    for (unsigned int n = 0; n < _num_other_b; ++n)
      factor += (*_d2F1dc1db1[n])[_qp] * (*_prop_db1detaj[n][etavar])[_qp];

    for (unsigned int n = 0; n < _num_j; ++n)
    {
      sum1 += (*_prop_dhjdetai[n])[_qp] * (*_prop_ci[n])[_qp];
      sum2 += (*_prop_d2hjdetaidetap[n][etavar])[_qp] * (*_prop_ci[n])[_qp] +
              (*_prop_dhjdetai[n])[_qp] * (*_prop_dcidetaj[n][etavar])[_qp];
    }
  }

  res += -_L[_qp] * (factor * sum1 + _first_df1[_qp] * sum2) * _phi[_j][_qp] * _test[_i][_qp];

  return res;
}
