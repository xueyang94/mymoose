//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KKSACBulkF.h"

registerMooseObject("PhaseFieldApp", KKSACBulkF);

InputParameters
KKSACBulkF::validParams()
{
  InputParameters params = KKSACBulkBase::validParams();
  params.addClassDescription("KKS model kernel (part 1 of 2) for the Bulk Allen-Cahn. This "
                             "includes all terms NOT dependent on chemical potential.");
  params.addRequiredCoupledVar("global_cs", "Global concentrations c, b, etc.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "ci_names",
      "Phase concentrations. These must have the same order as Fj_names and global_cs, for "
      "example, c1, c2, b1, b2.");
  params.addParam<std::vector<MaterialPropertyName>>(
      "dcidb_names",
      "Coupled dcidb in the order of dc1dc, dc2dc, dc1db, dc2db, db1dc, db2dc, db1db, db2db. These "
      "must have the same order as Fj_names and ci_names");
  params.addParam<std::vector<MaterialPropertyName>>(
      "dcideta_names",
      "The phase concentrations taken derivatives wrt kernel variable. ci must match the order in "
      "ci_names, for example, dc1deta, dc2deta, db1deta, db2deta, etc");
  params.addRequiredParam<MaterialPropertyName>(
      "fb_name",
      "Base name of the free energy function F (f_base in the corresponding KKSBaseMaterial)");
  params.addParam<MaterialPropertyName>(
      "g_name", "g", "Base name for the double well function g(eta)");
  params.addRequiredParam<MaterialPropertyName>("L_name", "The name of the Allen-Cahn mobility");
  params.addRequiredParam<Real>("w", "Double well height parameter");
  return params;
}

KKSACBulkF::KKSACBulkF(const InputParameters & parameters)
  : KKSACBulkBase(parameters),
    _c_map(getParameterJvarMap("global_cs")),
    _num_c(coupledComponents("global_cs")),
    _num_j(2),
    _ci_names(getParam<std::vector<MaterialPropertyName>>("ci_names")),
    _ci_name_matrix(_num_c),
    _dcidb_names(getParam<std::vector<MaterialPropertyName>>("dcidb_names")),
    _prop_dcidb(_num_c),
    _dcideta_names(getParam<std::vector<MaterialPropertyName>>("dcideta_names")),
    _prop_dcideta(_num_c),
    _Fa_name(getParam<MaterialPropertyName>("fa_name")),
    _first_dFa(_num_c),
    _Fb_name(getParam<MaterialPropertyName>("fb_name")),
    _prop_Fb(getMaterialProperty<Real>("fb_name")),
    _first_dFb(_num_c),
    _prop_dg(getMaterialPropertyDerivative<Real>("g_name", _eta_name)),
    _prop_d2g(getMaterialPropertyDerivative<Real>("g_name", _eta_name, _eta_name)),
    _L(getMaterialProperty<Real>("L_name")),
    _w(getParam<Real>("w"))
{
  // initialize _ci_name_matrix
  for (unsigned int m = 0; m < _num_c; ++m)
  {
    _ci_name_matrix[m].resize(_num_j);

    for (unsigned int n = 0; n < _num_j; ++n)
    {
      _ci_name_matrix[m][n] = _ci_names[m * _num_j + n];
    }
  }

  // declare _prop_dcidb. m is the numerator species (ci or bi), n is the phase of the numerator i,
  // l is the denominator species (c or b)
  for (unsigned int m = 0; m < _num_c; ++m)
  {
    _prop_dcidb[m].resize(_num_j);

    for (unsigned int n = 0; n < _num_j; ++n)
    {
      _prop_dcidb[m][n].resize(_num_c);

      for (unsigned int l = 0; l < _num_c; ++l)
      {
        _prop_dcidb[m][n][l] =
            &getMaterialPropertyByName<Real>(_dcidb_names[m * _num_j * _num_c + n + l * _num_j]);
      }
    }
  }

  // declare _prop_dcideta. m is te numerator species (ci or bi), n is the phase of the numerator
  // i
  for (unsigned int m = 0; m < _num_c; ++m)
  {
    _prop_dcideta[m].resize(_num_j);

    for (unsigned int n = 0; n < _num_j; ++n)
    {
      _prop_dcideta[m][n] = &getMaterialPropertyByName<Real>(_dcideta_names[m * _num_j + n]);
    }
  }

  // initialize _first_dFa
  for (unsigned int m = 0; m < _num_c; ++m)
    _first_dFa[m] = &getMaterialPropertyDerivative<Real>(_Fa_name, _ci_name_matrix[m][0]);

  // initialize _first_dFb
  for (unsigned int m = 0; m < _num_c; ++m)
    _first_dFb[m] = &getMaterialPropertyDerivative<Real>(_Fb_name, _ci_name_matrix[m][1]);
}

Real
KKSACBulkF::computeDFDOP(PFFunctionType type)
{
  switch (type)
  {
    case Residual:

      return -_prop_dh[_qp] * (_prop_Fa[_qp] - _prop_Fb[_qp]) + _w * _prop_dg[_qp];

    case Jacobian:

      Real sum1 = 0.0;
      for (unsigned int m = 0; m < _num_c; ++m)
        sum1 += (*_first_dFa[m])[_qp] * (*_prop_dcideta[m][0])[_qp] -
                (*_first_dFb[m])[_qp] * (*_prop_dcideta[m][1])[_qp];

      return (-(_prop_d2h[_qp] * (_prop_Fa[_qp] - _prop_Fb[_qp]) + _prop_dh[_qp] * sum1) +
              _w * _prop_d2g[_qp]) *
             _phi[_j][_qp];
  }

  mooseError("Invalid type passed in");
}

Real
KKSACBulkF::computeQpOffDiagJacobian(unsigned int jvar)
{
  // c is the coupled variable
  auto compvar = mapJvarToCvar(jvar, _c_map);
  if (compvar >= 0)
  {
    Real sum1 = 0.0;
    for (unsigned int m = 0; m < _num_c; ++m)
      sum1 += (*_first_dFa[m])[_qp] * (*_prop_dcidb[m][0][compvar])[_qp] -
              (*_first_dFb[m])[_qp] * (*_prop_dcidb[m][1][compvar])[_qp];

    return -_L[_qp] * _prop_dh[_qp] * sum1 * _phi[_j][_qp] * _test[_i][_qp];
  }

  return 0.0;
}
