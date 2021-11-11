//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PhaseConcentrationDerivatives.h"
#include "MatrixTools.h"

registerMooseObject("PhaseFieldApp", PhaseConcentrationDerivatives);

InputParameters
PhaseConcentrationDerivatives::validParams()
{
  InputParameters params = DerivativeMaterialInterface<Material>::validParams();
  params.addClassDescription(
      "Computes the KKS phase concentration derivatives wrt global concentrations and order "
      "parameters, which are used in the chain rules in the KKS kernels.");
  params.addRequiredCoupledVar("all_cs", "The interpolated concentration c, b, etc");
  params.addRequiredCoupledVar("all_etas", "Vector of all order parameters for all phases");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "Fj_names", "Phase energies in the same order as the all_etas.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "hj_names", "Names of the switching functions in the same order of the all_etas");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "ci_names",
      "Phase concentrations. The phase order must match Fi_names and global_c, for example, c1, "
      "c2, b1, b2, etc");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "dcidb_names",
      "The derivative of phase concentrations wrt global concentrations. The order must match "
      "Fj_names and global_c, for example, dc1dc, dc2dc, dc1db, dc2db, db1dc, db2dc, db1db, "
      "db2db, etc");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "dcidetaj_names",
      "The derivative of phase concentrations wrt order parameters. The order must match "
      "ci_names and Fj_names, for example, dc1deta1, dc1deta2, dc2deta1, dc2deta2, db1deta1, "
      "db1deta2, db2deta1, db2deta2, etc.");
  return params;
}

PhaseConcentrationDerivatives::PhaseConcentrationDerivatives(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _num_c(coupledComponents("all_cs")),
    _c_names(coupledComponents("all_cs")),
    _prop_ci(_num_c),
    _eta_names(coupledNames("all_etas")),
    _num_eta(coupledComponents("all_etas")),
    _Fj_names(getParam<std::vector<MaterialPropertyName>>("Fj_names")),
    _prop_d2Fjdcjdbj(_num_eta),
    _hj_names(getParam<std::vector<MaterialPropertyName>>("hj_names")),
    _prop_hj(_num_eta),
    _prop_dhjdetai(_num_eta),
    _ci_names(getParam<std::vector<MaterialPropertyName>>("ci_names")),
    _ci_name_matrix(_num_c),
    _dcidb_names(getParam<std::vector<MaterialPropertyName>>("dcidb_names")),
    _prop_dcidb(_num_c),
    _dcidetaj_names(getParam<std::vector<MaterialPropertyName>>("dcidetaj_names")),
    _prop_dcidetaj(_num_c)

{
  // declare _ci_name_matrix
  for (unsigned int m = 0; m < _num_c; ++m)
  {
    _ci_name_matrix[m].resize(_num_eta);

    for (unsigned int n = 0; n < _num_eta; ++n)
    {
      _ci_name_matrix[m][n] = _ci_names[m * _num_eta + n];
    }
  }

  // declare _prop_ci
  for (unsigned int m = 0; m < _num_c; ++m)
  {
    _prop_ci[m].resize(_num_eta);

    for (unsigned int n = 0; n < _num_eta; ++n)
    {
      _prop_ci[m][n] = &getMaterialPropertyByName<Real>(_ci_names[m * _num_eta + n]);
    }
  }

  for (unsigned int m = 0; m < _num_eta; ++m)
  {
    // declare _prop_hj
    _prop_hj[m] = &getMaterialPropertyByName<Real>(_hj_names[m]);

    _prop_dhjdetai[m].resize(_num_eta);

    // declare _prop_dhjetai
    for (unsigned int n = 0; n < _num_eta; ++n)
    {
      _prop_dhjdetai[m][n] = &getMaterialPropertyDerivative<Real>(_hj_names[m], _eta_names[n]);
    }
  }

  // declare _prop_d2Fjdcjdbj, m is phase, n is the first phase concentration species in phase m
  // (cm), l is the second phase concentration species in phase m (bm)
  for (unsigned int m = 0; m < _num_eta; ++m)
  {
    _prop_d2Fjdcjdbj[m].resize(_num_c);

    for (unsigned int n = 0; n < _num_c; ++n)
    {
      _prop_d2Fjdcjdbj[m][n].resize(_num_c);

      for (unsigned int l = 0; l < _num_c; ++l)
      {
        _prop_d2Fjdcjdbj[m][n][l] = &getMaterialPropertyDerivative<Real>(
            _Fj_names[m], _ci_name_matrix[n][m], _ci_name_matrix[l][m]);
      }
    }
  }

  // declare _prop_dcidb. m is the numerator species (ci or bi), n is the denominator species (c or
  // b), l is the phase of the numerator
  for (unsigned int m = 0; m < _num_c; ++m)
  {
    _prop_dcidb[m].resize(_num_c);

    for (unsigned int n = 0; n < _num_c; ++n)
    {
      _prop_dcidb[m][n].resize(_num_eta);

      for (unsigned int l = 0; l < _num_eta; ++l)
      {
        _prop_dcidb[m][n][l] =
            &declareProperty<Real>(_dcidb_names[m * (_num_c * _num_eta) + n * _num_eta + l]);
      }
    }
  }

  // declare _prop_dcidetaj
  for (unsigned int m = 0; m < _num_c; ++m)
  {
    _prop_dcidetaj[m].resize(_num_eta);

    for (unsigned int n = 0; n < _num_eta; ++n)
    {
      _prop_dcidetaj[m][n].resize(_num_eta);

      for (unsigned int l = 0; l < _num_eta; ++l)
      {
        _prop_dcidetaj[m][n][l] =
            &declareProperty<Real>(_dcidetaj_names[m * (_num_eta * _num_eta) + n + l * _num_eta]);
      }
    }
  }
}

void
PhaseConcentrationDerivatives::computeQpProperties()
{
  //////////////////////////////////////////////////// solve linear system of constraint derivatives wrt c for computing dcidc and dbidc
  std::vector<std::vector<Real>> A_c(6);
  for (auto & row : A_c)
    row.resize(6);

  A_c[0][0] = (*_prop_d2Fjdcjdbj[0][0][0])[_qp];
  A_c[0][1] = -(*_prop_d2Fjdcjdbj[1][1][1])[_qp];
  A_c[0][2] = 0;
  A_c[0][3] = (*_prop_d2Fjdcjdbj[0][0][1])[_qp];
  A_c[0][4] = -(*_prop_d2Fjdcjdbj[1][0][1])[_qp];
  A_c[0][5] = 0;
  A_c[1][0] = 0;
  A_c[1][1] = (*_prop_d2Fjdcjdbj[1][0][0])[_qp];
  A_c[1][2] = -(*_prop_d2Fjdcjdbj[2][0][0])[_qp];
  A_c[1][3] = 0;
  A_c[1][4] = (*_prop_d2Fjdcjdbj[1][0][1])[_qp];
  A_c[1][5] = -(*_prop_d2Fjdcjdbj[2][0][1])[_qp];
  A_c[2][0] = (*_prop_hj[0])[_qp];
  A_c[2][1] = (*_prop_hj[1])[_qp];
  A_c[2][2] = (*_prop_hj[2])[_qp];
  A_c[2][3] = 0;
  A_c[2][4] = 0;
  A_c[2][5] = 0;

  A_c[3][0] = (*_prop_d2Fjdcjdbj[0][1][0])[_qp];
  A_c[3][1] = -(*_prop_d2Fjdcjdbj[1][1][0])[_qp];
  A_c[3][2] = 0;
  A_c[3][3] = (*_prop_d2Fjdcjdbj[0][1][1])[_qp];
  A_c[3][4] = -(*_prop_d2Fjdcjdbj[1][1][1])[_qp];
  A_c[3][5] = 0;
  A_c[4][0] = 0;
  A_c[4][1] = (*_prop_d2Fjdcjdbj[1][1][0])[_qp];
  A_c[4][2] = -(*_prop_d2Fjdcjdbj[2][1][0])[_qp];
  A_c[4][3] = 0;
  A_c[4][4] = (*_prop_d2Fjdcjdbj[1][1][1])[_qp];
  A_c[4][5] = -(*_prop_d2Fjdcjdbj[2][1][1])[_qp];
  A_c[5][0] = 0;
  A_c[5][1] = 0;
  A_c[5][2] = 0;
  A_c[5][3] = (*_prop_hj[0])[_qp];
  A_c[5][4] = (*_prop_hj[1])[_qp];
  A_c[5][5] = (*_prop_hj[2])[_qp];

  MatrixTools::inverse(A_c, A_c);

  // compute dcidc and dbidc
  std::vector<Real> x_dcondc(6);
  std::vector<Real> k_dcondc{0, 0, 1, 0, 0, 0};

  for (unsigned int i = 0; i < 6; ++i)
  {
    x_dcondc[i] = A_c[i][0] * k_dcondc[0] + A_c[i][1] * k_dcondc[1] + A_c[i][2] * k_dcondc[2] +
                  A_c[i][3] * k_dcondc[3] + A_c[i][4] * k_dcondc[4] + A_c[i][5] * k_dcondc[5];
  }

  (*_prop_dcidb[0][0][0])[_qp] = x_dcondc[0];
  (*_prop_dcidb[0][0][1])[_qp] = x_dcondc[1];
  (*_prop_dcidb[0][0][2])[_qp] = x_dcondc[2];
  (*_prop_dcidb[1][0][0])[_qp] = x_dcondc[3];
  (*_prop_dcidb[1][0][1])[_qp] = x_dcondc[4];
  (*_prop_dcidb[1][0][2])[_qp] = x_dcondc[5];

  ////////////////////////////////////////////////// solve linear system of constraint derivatives wrt b for computing dcidb and dbidb
  std::vector<std::vector<Real>> A_b(6);
  for (auto & row : A_b)
    row.resize(6);

  A_b[0][0] = (*_prop_d2Fjdcjdbj[0][1][0])[_qp];
  A_b[0][1] = -(*_prop_d2Fjdcjdbj[1][1][0])[_qp];
  A_b[0][2] = 0;
  A_b[0][3] = (*_prop_d2Fjdcjdbj[0][1][1])[_qp];
  A_b[0][4] = -(*_prop_d2Fjdcjdbj[1][1][1])[_qp];
  A_b[0][5] = 0;
  A_b[1][0] = 0;
  A_b[1][1] = (*_prop_d2Fjdcjdbj[1][1][0])[_qp];
  A_b[1][2] = -(*_prop_d2Fjdcjdbj[2][1][0])[_qp];
  A_b[1][3] = 0;
  A_b[1][4] = (*_prop_d2Fjdcjdbj[1][1][1])[_qp];
  A_b[1][5] = -(*_prop_d2Fjdcjdbj[2][1][1])[_qp];
  A_b[2][0] = 0;
  A_b[2][1] = 0;
  A_b[2][2] = 0;
  A_b[2][3] = (*_prop_hj[0])[_qp];
  A_b[2][4] = (*_prop_hj[1])[_qp];
  A_b[2][5] = (*_prop_hj[2])[_qp];

  A_b[3][0] = (*_prop_d2Fjdcjdbj[0][0][0])[_qp];
  A_b[3][1] = -(*_prop_d2Fjdcjdbj[1][0][0])[_qp];
  A_b[3][2] = 0;
  A_b[3][3] = (*_prop_d2Fjdcjdbj[0][0][1])[_qp];
  A_b[3][4] = -(*_prop_d2Fjdcjdbj[1][0][1])[_qp];
  A_b[3][5] = 0;
  A_b[4][0] = 0;
  A_b[4][1] = (*_prop_d2Fjdcjdbj[1][0][0])[_qp];
  A_b[4][2] = -(*_prop_d2Fjdcjdbj[2][0][0])[_qp];
  A_b[4][3] = 0;
  A_b[4][4] = (*_prop_d2Fjdcjdbj[1][0][1])[_qp];
  A_b[4][5] = -(*_prop_d2Fjdcjdbj[2][0][1])[_qp];
  A_b[5][0] = (*_prop_hj[0])[_qp];
  A_b[5][1] = (*_prop_hj[1])[_qp];
  A_b[5][2] = (*_prop_hj[2])[_qp];
  A_b[5][3] = 0;
  A_b[5][4] = 0;
  A_b[5][5] = 0;

  MatrixTools::inverse(A_b, A_b);

  // compute dcidb and dbidb
  std::vector<Real> x_dcondb(6);
  std::vector<Real> k_dcondb{0, 0, 1, 0, 0, 0};

  for (unsigned int i = 0; i < 6; ++i)
  {
    x_dcondb[i] = A_b[i][0] * k_dcondb[0] + A_b[i][1] * k_dcondb[1] + A_b[i][2] * k_dcondb[2] +
                  A_b[i][3] * k_dcondb[3] + A_b[i][4] * k_dcondb[4] + A_b[i][5] * k_dcondb[5];
  }

  (*_prop_dcidb[0][1][0])[_qp] = x_dcondb[0];
  (*_prop_dcidb[0][1][1])[_qp] = x_dcondb[1];
  (*_prop_dcidb[0][1][2])[_qp] = x_dcondb[2];
  (*_prop_dcidb[1][1][0])[_qp] = x_dcondb[3];
  (*_prop_dcidb[1][1][1])[_qp] = x_dcondb[4];
  (*_prop_dcidb[1][1][2])[_qp] = x_dcondb[5];

  ////////////////////////////////////////////////// solve linear system of constraint derivatives wrt etaj for computing dcidetaj
  std::vector<std::vector<Real>> A_eta(6);
  for (auto & row : A_eta)
    row.resize(6);

  A_eta[0][0] = (*_prop_d2Fjdcjdbj[0][0][0])[_qp];
  A_eta[0][1] = -(*_prop_d2Fjdcjdbj[1][0][0])[_qp];
  A_eta[0][2] = 0;
  A_eta[0][3] = (*_prop_d2Fjdcjdbj[0][0][1])[_qp];
  A_eta[0][4] = -(*_prop_d2Fjdcjdbj[1][0][1])[_qp];
  A_eta[0][5] = 0;
  A_eta[1][0] = 0;
  A_eta[1][1] = (*_prop_d2Fjdcjdbj[1][0][0])[_qp];
  A_eta[1][2] = -(*_prop_d2Fjdcjdbj[2][0][0])[_qp];
  A_eta[1][3] = 0;
  A_eta[1][4] = (*_prop_d2Fjdcjdbj[1][0][1])[_qp];
  A_eta[1][5] = -(*_prop_d2Fjdcjdbj[2][0][1])[_qp];
  A_eta[2][0] = (*_prop_hj[0])[_qp];
  A_eta[2][1] = (*_prop_hj[1])[_qp];
  A_eta[2][2] = (*_prop_hj[2])[_qp];
  A_eta[2][3] = 0;
  A_eta[2][4] = 0;
  A_eta[2][5] = 0;

  A_eta[3][0] = (*_prop_d2Fjdcjdbj[0][1][0])[_qp];
  A_eta[3][1] = -(*_prop_d2Fjdcjdbj[1][1][0])[_qp];
  A_eta[3][2] = 0;
  A_eta[3][3] = (*_prop_d2Fjdcjdbj[0][1][1])[_qp];
  A_eta[3][4] = -(*_prop_d2Fjdcjdbj[1][1][1])[_qp];
  A_eta[3][5] = 0;
  A_eta[4][0] = 0;
  A_eta[4][1] = (*_prop_d2Fjdcjdbj[1][1][0])[_qp];
  A_eta[4][2] = -(*_prop_d2Fjdcjdbj[2][1][0])[_qp];
  A_eta[4][3] = 0;
  A_eta[4][4] = (*_prop_d2Fjdcjdbj[1][1][1])[_qp];
  A_eta[4][5] = -(*_prop_d2Fjdcjdbj[2][1][1])[_qp];
  A_eta[5][0] = 0;
  A_eta[5][1] = 0;
  A_eta[5][2] = 0;
  A_eta[5][3] = (*_prop_hj[0])[_qp];
  A_eta[5][4] = (*_prop_hj[1])[_qp];
  A_eta[5][5] = (*_prop_hj[2])[_qp];

  MatrixTools::inverse(A_eta, A_eta);

  for (unsigned int m = 0; m < _num_eta; ++m) // loop through all the order parameters that the
                                              // phase concentrations take derivative to
  {
    std::vector<Real> x_dcidetaj(6);
    std::vector<Real> k_dcidetaj{0,
                                 0,
                                 -((*_prop_dhjdetai[0][m])[_qp] * (*_prop_ci[0][m])[_qp] +
                                   (*_prop_dhjdetai[1][m])[_qp] * (*_prop_ci[0][m])[_qp] +
                                   (*_prop_dhjdetai[2][m])[_qp] * (*_prop_ci[0][m])[_qp]),
                                 0,
                                 0,
                                 -((*_prop_dhjdetai[0][m])[_qp] * (*_prop_ci[1][m])[_qp] +
                                   (*_prop_dhjdetai[1][m])[_qp] * (*_prop_ci[1][m])[_qp] +
                                   (*_prop_dhjdetai[2][m])[_qp] * (*_prop_ci[1][m])[_qp])};

    for (unsigned int i = 0; i < 6; ++i)
    {
      x_dcidetaj[i] = A_eta[i][0] * k_dcidetaj[0] + A_eta[i][1] * k_dcidetaj[1] +
                      A_eta[i][2] * k_dcidetaj[2] + A_eta[i][3] * k_dcidetaj[3] +
                      A_eta[i][4] * k_dcidetaj[4] + A_eta[i][5] * k_dcidetaj[5];
    }

    (*_prop_dcidetaj[0][0][m])[_qp] = x_dcidetaj[0];
    (*_prop_dcidetaj[0][1][m])[_qp] = x_dcidetaj[1];
    (*_prop_dcidetaj[0][2][m])[_qp] = x_dcidetaj[2];
    (*_prop_dcidetaj[1][0][m])[_qp] = x_dcidetaj[3];
    (*_prop_dcidetaj[1][1][m])[_qp] = x_dcidetaj[4];
    (*_prop_dcidetaj[1][2][m])[_qp] = x_dcidetaj[5];
  }
}
