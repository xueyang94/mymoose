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
    _d2Fjdcjdbj(_num_eta),
    _hj_names(getParam<std::vector<MaterialPropertyName>>("hj_names")),
    _prop_hj(_num_eta),
    _dhjdetap(_num_eta),
    _ci_names(getParam<std::vector<MaterialPropertyName>>("ci_names")),
    _ci_name_matrix(_num_c),
    _dcidb_names(getParam<std::vector<MaterialPropertyName>>("dcidb_names")),
    _prop_dcidb(_num_c),
    _dcidetaj_names(getParam<std::vector<MaterialPropertyName>>("dcidetaj_names")),
    _prop_dcidetaj(_num_c)
{
  // initialize _ci_name_matrix
  for (unsigned int m = 0; m < _num_c; ++m)
  {
    _ci_name_matrix[m].resize(_num_eta);

    for (unsigned int n = 0; n < _num_eta; ++n)
    {
      _ci_name_matrix[m][n] = _ci_names[m * _num_eta + n];
    }
  }

  // initialize _prop_ci
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
    // initialize _prop_hj
    _prop_hj[m] = &getMaterialPropertyByName<Real>(_hj_names[m]);

    _dhjdetap[m].resize(_num_eta);

    // initialize _prop_dhjetap
    for (unsigned int n = 0; n < _num_eta; ++n)
    {
      _dhjdetap[m][n] = &getMaterialPropertyDerivative<Real>(_hj_names[m], _eta_names[n]);
    }
  }

  // initialize _d2Fjdcjdbj[m][n][l], m is phase, n is the first phase concentration species in
  // phase m, l is the second phase concentration species in phase m
  for (unsigned int m = 0; m < _num_eta; ++m)
  {
    _d2Fjdcjdbj[m].resize(_num_c);

    for (unsigned int n = 0; n < _num_c; ++n)
    {
      _d2Fjdcjdbj[m][n].resize(_num_c);

      for (unsigned int l = 0; l < _num_c; ++l)
      {
        _d2Fjdcjdbj[m][n][l] = &getMaterialPropertyDerivative<Real>(
            _Fj_names[m], _ci_name_matrix[n][m], _ci_name_matrix[l][m]);
      }
    }
  }

  // declare _prop_dcidb. m is the numerator species (ci or bi), n is the phase of the numerator i,
  // l is the denominator species (c or b)
  for (unsigned int m = 0; m < _num_c; ++m)
  {
    _prop_dcidb[m].resize(_num_eta);

    for (unsigned int n = 0; n < _num_eta; ++n)
    {
      _prop_dcidb[m][n].resize(_num_c);

      for (unsigned int l = 0; l < _num_c; ++l)
      {
        _prop_dcidb[m][n][l] =
            &declareProperty<Real>(_dcidb_names[m * _num_eta * _num_c + n + l * _num_eta]);
      }
    }
  }

  // declare _prop_dcidetaj. m is the numerator species (ci or bi), n is the phase of the numerator
  // i, l is the phase of denominator j
  for (unsigned int m = 0; m < _num_c; ++m)
  {
    _prop_dcidetaj[m].resize(_num_eta);

    for (unsigned int n = 0; n < _num_eta; ++n)
    {
      _prop_dcidetaj[m][n].resize(_num_eta);

      for (unsigned int l = 0; l < _num_eta; ++l)
      {
        _prop_dcidetaj[m][n][l] =
            &declareProperty<Real>(_dcidetaj_names[m * _num_eta * _num_eta + n + l * _num_eta]);
      }
    }
  }
}

void
PhaseConcentrationDerivatives::computeQpProperties()
{
  /////////////////////////////////////////////////////////////////////////////////////////////////// solve linear system of constraint derivatives wrt c for computing dcidc and dbidc
  // loop through taking derivative wrt the ith component, each i constructs a A_c
  for (unsigned int i = 0; i < _num_c; ++i)
  {
    // declare A_c
    std::vector<std::vector<Real>> A_c(_num_eta * _num_c);
    for (auto & row : A_c)
      row.resize(_num_eta * _num_c);

    // initialize all elements in A_c to be zero
    for (unsigned int m = 0; m < _num_eta * _num_c; ++m)
    {
      for (unsigned int n = 0; n < _num_eta * _num_c; ++n)
        A_c[m][n] = 0;
    }

    // now fill in the non-zero elements in A_c
    // first assign the elements in A that come from the mu equality derivative equations
    // loop through the constraint equation sets of the mth component
    for (unsigned int m = 0; m < _num_c; ++m)
    {
      // loop through the nth constraint equation in the constraint set of one component
      for (unsigned int n = 0; n < (_num_eta - 1); ++n)
      {
        // loop through the lth chain rule terms in one constrain equation
        for (unsigned int l = 0; l < _num_c; ++l)
        {
          A_c[m * _num_eta + n][n + l * _num_eta] = (*_d2Fjdcjdbj[n][m][l])[_qp];
          A_c[m * _num_eta + n][n + l * _num_eta + 1] = -(*_d2Fjdcjdbj[n + 1][m][l])[_qp];
        }
      }
    }

    // then assign the elements in A that come from the concentration conservation equations
    for (unsigned int m = 0; m < _num_c; ++m)
    {
      for (unsigned int n = 0; n < _num_eta; ++n)
        A_c[m * _num_eta + _num_c][m * _num_eta + n] = (*_prop_hj[n])[_qp];
    }

    MatrixTools::inverse(A_c, A_c);

    std::vector<Real> k_c(_num_eta * _num_c); // of component i
    std::vector<Real> x_c(_num_eta * _num_c); // of component i

    // initialize all elements in k_c to be zero
    for (unsigned int m = 0; m < (_num_eta * _num_c); ++m)
      k_c[m] = 0;

    // assign non-zero elements in k_c
    k_c[i * _num_eta + _num_c] = 1;

    // compute x_c
    for (unsigned int m = 0; m < (_num_eta * _num_c); ++m)
    {
      for (unsigned int n = 0; n < (_num_eta * _num_c); ++n)
        x_c[m] += A_c[m][n] * k_c[n];
    }

    // assign the values in x_c to _prop_dcidb
    for (unsigned int m = 0; m < _num_c; ++m)
    {
      for (unsigned int n = 0; n < _num_eta; ++n)
      {
        (*_prop_dcidb[m][n][i])[_qp] = x_c[m * _num_eta + n];
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////// solve linear system of constraint derivatives wrt etaj for computing dcidetaj
  std::vector<std::vector<Real>> A_eta(_num_eta * _num_c);
  for (auto & row : A_eta)
    row.resize(_num_eta * _num_c);

  // initialize all elements in A to be zero
  for (unsigned int m = 0; m < _num_eta * _num_c; ++m)
  {
    for (unsigned int n = 0; n < _num_eta * _num_c; ++n)
      A_eta[m][n] = 0;
  }

  // now fill in the non-zero elements in A_c
  // first assign the elements in A that come from the mu equality derivative equations
  for (unsigned int m = 0; m < _num_c; ++m)
  {
    for (unsigned int n = 0; n < (_num_eta - 1); ++n)
    {
      for (unsigned int l = 0; l < _num_c; ++l)
      {
        A_eta[m * _num_eta + n][n + l * _num_eta] = (*_d2Fjdcjdbj[n][m][l])[_qp];
        A_eta[m * _num_eta + n][n + l * _num_eta + 1] = -(*_d2Fjdcjdbj[n + 1][m][l])[_qp];
      }
    }
  }

  // then assign the elements in A that come from the concentration conservation equations
  for (unsigned int m = 0; m < _num_c; ++m)
  {
    for (unsigned int n = 0; n < _num_eta; ++n)
      A_eta[m * _num_eta + _num_c][m * _num_eta + n] = (*_prop_hj[n])[_qp];
  }

  MatrixTools::inverse(A_eta, A_eta);

  for (unsigned int i = 0; i < _num_eta; ++i) // loop through all the order parameters that the
                                              // phase concentrations take derivative to
  {
    std::vector<Real> k_eta(_num_eta * _num_c); // of phase i
    std::vector<Real> x_eta(_num_eta * _num_c); // of phase i

    // initialize all elements in k_eta to be zero
    for (unsigned int m = 0; m < (_num_eta * _num_c); ++m)
      k_eta[m] = 0;

    // assign non-zero elements in k_eta
    for (unsigned int m = 0; m < _num_c; ++m)
    {
      Real sum = 0.0;

      for (unsigned int n = 0; n < _num_eta; ++n)
        sum += (*_dhjdetap[n][i])[_qp] * (*_prop_ci[m][n])[_qp];

      k_eta[m * _num_eta + _num_c] = -sum;
    }

    // compute x_eta
    for (unsigned int m = 0; m < (_num_eta * _num_c); ++m)
    {
      for (unsigned int n = 0; n < (_num_eta * _num_c); ++n)
        x_eta[m] += A_eta[m][n] * k_eta[n];
    }

    // assign the values in x_eta to _prop_dcidetaj
    for (unsigned int m = 0; m < _num_c; ++m)
    {
      for (unsigned int n = 0; n < _num_eta; ++n)
      {
        (*_prop_dcidetaj[m][n][i])[_qp] = x_eta[m * _num_eta + n];
      }
    }
  }
}
