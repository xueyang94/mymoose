//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KKSPhaseConcentrationMultiPhaseDerivatives.h"
#include "MatrixTools.h"

registerMooseObject("PhaseFieldApp", KKSPhaseConcentrationMultiPhaseDerivatives);

InputParameters
KKSPhaseConcentrationMultiPhaseDerivatives::validParams()
{
  InputParameters params = DerivativeMaterialInterface<Material>::validParams();
  params.addClassDescription(
      "Computes the KKS phase concentration derivatives wrt global concentrations and order "
      "parameters, which are used for the chain rule in the KKS kernels.");
  params.addRequiredCoupledVar("global_cs", "Global concentrations, for example, c, b.");
  params.addRequiredCoupledVar(
      "all_etas", "Order parameters for all phases. Place in the same order as Fj_names.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "ci_names",
      "Phase concentrations. These must have the same order as Fj_names and global_cs, for "
      "example, c1, c2, b1, b2.");
  params.addParam<std::vector<MaterialPropertyName>>(
      "dcidb_names",
      "The derivative of phase concentration wrt global concentration. They must have the same "
      "order as Fj_names and ci_names, for example, dc1dc, dc2dc, dc1db, dc2db, db1dc, db2dc, "
      "db1db, db2db.");
  params.addParam<std::vector<MaterialPropertyName>>(
      "dcidetaj_names",
      "THe derivative of phase concentration wrt order parameter. ci must have the order as "
      "global_c and all_etas, and etaj must match the order in all_etas, for example, dc1deta1, "
      "dc2deta1, dc1deta2, dc2deta2, db1deta1, db2deta1, db1deta2, db2deta2.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "hj_names", "Names of the switching functions in the same order of the all_etas");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "Fj_names", "Phase energies in the same order as the all_etas.");
  return params;
}

KKSPhaseConcentrationMultiPhaseDerivatives::KKSPhaseConcentrationMultiPhaseDerivatives(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _num_c(coupledComponents("global_cs")),
    _eta_names(coupledNames("all_etas")),
    _num_j(coupledComponents("all_etas")),
    _prop_ci(_num_c * _num_j),
    _ci_names(getParam<std::vector<MaterialPropertyName>>("ci_names")),
    _dcidb_names(getParam<std::vector<MaterialPropertyName>>("dcidb_names")),
    _prop_dcidb(_num_c),
    _dcidetaj_names(getParam<std::vector<MaterialPropertyName>>("dcidetaj_names")),
    _prop_dcidetaj(_num_c),
    _hj_names(getParam<std::vector<MaterialPropertyName>>("hj_names")),
    _prop_hj(_num_j),
    _dhjdetap(_num_j),
    _Fj_names(getParam<std::vector<MaterialPropertyName>>("Fj_names")),
    _d2Fjdcjdbj(_num_j)
{
  for (unsigned int m = 0; m < _num_c * _num_j; ++m)
    _prop_ci[m] = &getMaterialPropertyByName<Real>(_ci_names[m]);

  for (unsigned int m = 0; m < _num_c; ++m)
  {
    _prop_dcidb[m].resize(_num_j);
    _prop_dcidetaj[m].resize(_num_j);

    for (unsigned int n = 0; n < _num_j; ++n)
    {
      _prop_dcidb[m][n].resize(_num_c);
      _prop_dcidetaj[m][n].resize(_num_j);

      // declare _prop_dcidb. m is the numerator species, n is the phase of the numerator i, l is
      // the denominator species
      for (unsigned int l = 0; l < _num_c; ++l)
        _prop_dcidb[m][n][l] =
            &declareProperty<Real>(_dcidb_names[m * _num_j * _num_c + n + l * _num_j]);

      // declare _prop_dcidetaj. m is the numerator species, n is the phase of the numerator i, l is
      // the phase of denominator j
      for (unsigned int l = 0; l < _num_j; ++l)
        _prop_dcidetaj[m][n][l] =
            &declareProperty<Real>(_dcidetaj_names[m * _num_j * _num_j + n + l * _num_j]);
    }
  }

  for (unsigned int m = 0; m < _num_j; ++m)
  {
    _prop_hj[m] = &getMaterialPropertyByName<Real>(_hj_names[m]);

    _dhjdetap[m].resize(_num_j);

    for (unsigned int n = 0; n < _num_j; ++n)
      _dhjdetap[m][n] = &getMaterialPropertyDerivative<Real>(_hj_names[m], _eta_names[n]);
  }

  // initialize _d2Fjdcjdbj[m][n][l], m is phase, n is the first phase concentration species in
  // phase m, l is the second phase concentration species in phase m
  for (unsigned int m = 0; m < _num_j; ++m)
  {
    _d2Fjdcjdbj[m].resize(_num_c);

    for (unsigned int n = 0; n < _num_c; ++n)
    {
      _d2Fjdcjdbj[m][n].resize(_num_c);

      for (unsigned int l = 0; l < _num_c; ++l)
        _d2Fjdcjdbj[m][n][l] = &getMaterialPropertyDerivative<Real>(
            _Fj_names[m], _ci_names[m + n * _num_j], _ci_names[m + l * _num_j]);
    }
  }
}

void
KKSPhaseConcentrationMultiPhaseDerivatives::computeQpProperties()
{
  // declare Jacobian matrix A
  std::vector<std::vector<Real>> A(_num_j * _num_c);
  for (auto & row : A)
    row.resize(_num_j * _num_c);

  // initialize all elements in A_c to be zero
  for (unsigned int m = 0; m < _num_j * _num_c; ++m)
  {
    for (unsigned int n = 0; n < _num_j * _num_c; ++n)
      A[m][n] = 0;
  }

  // now fill in the non-zero elements in A
  // first assign the elements in A that come from the mu equality derivative equations
  // loop through the constraint equation sets of the mth component
  for (unsigned int m = 0; m < _num_c; ++m)
  {
    // loop through the nth constraint equation in the constraint set of one component
    for (unsigned int n = 0; n < (_num_j - 1); ++n)
    {
      // loop through the lth chain rule terms in one constrain equation
      for (unsigned int l = 0; l < _num_c; ++l)
      {
        A[m * _num_j + n][n + l * _num_j] = (*_d2Fjdcjdbj[n][m][l])[_qp];
        A[m * _num_j + n][n + l * _num_j + 1] = -(*_d2Fjdcjdbj[n + 1][m][l])[_qp];
      }
    }
  }

  // then assign the elements in A that come from the concentration conservation equations
  for (unsigned int m = 0; m < _num_c; ++m)
  {
    for (unsigned int n = 0; n < _num_j; ++n)
      A[(m + 1) * _num_j - 1][m * _num_j + n] = (*_prop_hj[n])[_qp];
  }

  MatrixTools::inverse(A, A);

  ///////////////////////////////////////////////////////////// solve linear system of constraint derivatives wrt c for computing dcidc and dbidc
  // loop through taking derivative wrt the ith component
  for (unsigned int i = 0; i < _num_c; ++i)
  {
    std::vector<Real> k_c(_num_j * _num_c);

    // assign non-zero elements in k_c
    k_c[i * _num_j + _num_j - 1] = 1;

    std::vector<Real> x_c(_num_j * _num_c);

    // compute x_c
    for (unsigned int m = 0; m < (_num_j * _num_c); ++m)
    {
      for (unsigned int n = 0; n < (_num_j * _num_c); ++n)
        x_c[m] += A[m][n] * k_c[n];
    }

    // assign the values in x_c to _prop_dcidb
    for (unsigned int m = 0; m < _num_c; ++m)
    {
      for (unsigned int n = 0; n < _num_j; ++n)
        (*_prop_dcidb[m][n][i])[_qp] = x_c[m * _num_j + n];
    }
  }

  ///////////////////////////////////////////////////////////// solve linear system of constraint derivatives wrt etaj for computing dcidetaj
  for (unsigned int i = 0; i < _num_j; ++i) // loop through all the order parameters that the
                                            // phase concentrations take derivative to
  {
    std::vector<Real> k_eta(_num_j * _num_c);

    // assign non-zero elements in k_eta
    for (unsigned int m = 0; m < _num_c; ++m)
    {
      Real sum = 0.0;

      for (unsigned int n = 0; n < _num_j; ++n)
        sum += (*_dhjdetap[n][i])[_qp] * (*_prop_ci[m * _num_j + n])[_qp];

      k_eta[m * _num_j + _num_j - 1] = -sum;
    }

    std::vector<Real> x_eta(_num_j * _num_c);

    // compute x_eta
    for (unsigned int m = 0; m < (_num_j * _num_c); ++m)
    {
      for (unsigned int n = 0; n < (_num_j * _num_c); ++n)
        x_eta[m] += A[m][n] * k_eta[n];
    }

    // assign the values in x_eta to _prop_dcidetaj
    for (unsigned int m = 0; m < _num_c; ++m)
    {
      for (unsigned int n = 0; n < _num_j; ++n)
        (*_prop_dcidetaj[m][n][i])[_qp] = x_eta[m * _num_j + n];
    }
  }

  // std::cout << "dc1dc " << (*_prop_dcidb[0][0][0])[_qp] << std::endl;
  // std::cout << "dc2dc " << (*_prop_dcidb[0][1][0])[_qp] << std::endl;
  // std::cout << "dc3dc " << (*_prop_dcidb[0][2][0])[_qp] << std::endl;
  // std::cout << "db1dc " << (*_prop_dcidb[1][0][0])[_qp] << std::endl;
  // std::cout << "db2dc " << (*_prop_dcidb[1][1][0])[_qp] << std::endl;
  // std::cout << "db3dc " << (*_prop_dcidb[1][2][0])[_qp] << std::endl;
  // std::cout << "dc1db " << (*_prop_dcidb[0][0][1])[_qp] << std::endl;
  // std::cout << "dc2db " << (*_prop_dcidb[0][1][1])[_qp] << std::endl;
  // std::cout << "dc3db " << (*_prop_dcidb[0][2][1])[_qp] << std::endl;
  // std::cout << "db1db " << (*_prop_dcidb[1][0][1])[_qp] << std::endl;
  // std::cout << "db2db " << (*_prop_dcidb[1][1][1])[_qp] << std::endl;
  // std::cout << "db3db " << (*_prop_dcidb[1][2][1])[_qp] << std::endl;
  // std::cout << '\n' << std::endl;
}
