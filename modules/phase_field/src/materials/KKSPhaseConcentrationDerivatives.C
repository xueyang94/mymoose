//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KKSPhaseConcentrationDerivatives.h"
#include "MatrixTools.h"

registerMooseObject("PhaseFieldApp", KKSPhaseConcentrationDerivatives);

InputParameters
KKSPhaseConcentrationDerivatives::validParams()
{
  InputParameters params = DerivativeMaterialInterface<Material>::validParams();
  params.addClassDescription(
      "Computes the KKS phase concentration derivatives wrt global concentrations and order "
      "parameters, which are used in the chain rules in the KKS kernels.");
  params.addRequiredCoupledVar("global_cs", "The interpolated concentration c, b, etc");
  params.addRequiredCoupledVar("eta", "Order parameter");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "ci_names",
      "Phase concentrations. The order must match Fj_names and global_c, for example, c1, "
      "c2, b1, b2, etc");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "dcidb_names",
      "The derivative of phase concentrations wrt global concentrations. The order must match "
      "Fj_names and ci_names, for example, dc1dc, dc2dc, dc1db, dc2db, db1dc, db2dc, db1db, "
      "db2db, etc");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "dcideta_names",
      "The derivative of phase concentrations wrt the order parameter. The order must match "
      "ci_names, for example, dc1deta, dc2deta, db1deta, db2deta, etc.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>("Fj_names", "F1 and F2");
  params.addParam<MaterialPropertyName>(
      "h_name", "h", "Base name for the switching function h(eta)");
  return params;
}

KKSPhaseConcentrationDerivatives::KKSPhaseConcentrationDerivatives(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _num_c(coupledComponents("global_cs")),
    _c_names(coupledComponents("global_cs")),
    _prop_ci(_num_c * 2),
    _eta_name(getVar("eta", 0)->name()),
    _ci_names(getParam<std::vector<MaterialPropertyName>>("ci_names")),
    _dcidb_names(getParam<std::vector<MaterialPropertyName>>("dcidb_names")),
    _prop_dcidb(_num_c),
    _dcideta_names(getParam<std::vector<MaterialPropertyName>>("dcideta_names")),
    _prop_dcideta(_num_c),
    _Fj_names(getParam<std::vector<MaterialPropertyName>>("Fj_names")),
    _d2Fjdcjdbj(2),
    _prop_h(getMaterialProperty<Real>("h_name")),
    _prop_dh(getMaterialPropertyDerivative<Real>("h_name", _eta_name))
{
  // initialize _prop_ci
  for (unsigned int m = 0; m < _num_c * 2; ++m)
    _prop_ci[m] = &getMaterialPropertyByName<Real>(_ci_names[m]);

  // declare _prop_dcidb. m is the numerator species (ci or bi), n is the phase of the numerator i,
  // l is the denominator species (c or b)
  for (unsigned int m = 0; m < _num_c; ++m)
  {
    _prop_dcidb[m].resize(2);
    _prop_dcideta[m].resize(2);

    for (unsigned int n = 0; n < 2; ++n)
    {
      _prop_dcideta[m][n] = &declareProperty<Real>(_dcideta_names[m * 2 + n]);

      _prop_dcidb[m][n].resize(_num_c);

      for (unsigned int l = 0; l < _num_c; ++l)
      {
        _prop_dcidb[m][n][l] = &declareProperty<Real>(_dcidb_names[m * 2 * _num_c + n + l * 2]);
      }
    }
  }

  // initialize _d2Fjdcjdbj[m][n][l], m is phase, n is the first phase concentration species in
  // phase m, l is the second phase concentration species in phase m
  for (unsigned int m = 0; m < 2; ++m)
  {
    _d2Fjdcjdbj[m].resize(_num_c);

    for (unsigned int n = 0; n < _num_c; ++n)
    {
      _d2Fjdcjdbj[m][n].resize(_num_c);

      for (unsigned int l = 0; l < _num_c; ++l)
      {
        _d2Fjdcjdbj[m][n][l] = &getMaterialPropertyDerivative<Real>(
            _Fj_names[m], _ci_names[m + n * 2], _ci_names[m + l * 2]);
      }
    }
  }
}

void
KKSPhaseConcentrationDerivatives::computeQpProperties()
{
  /////////////////////////////////////////////////////////////////////////////////////////////////// solve linear system of constraint derivatives wrt c for computing dcidb.
  // declare A
  std::vector<std::vector<Real>> A(_num_c * 2);
  for (auto & row : A)
    row.resize(_num_c * 2);

  // initialize all elements in A to be zero
  for (unsigned int m = 0; m < _num_c * 2; ++m)
  {
    for (unsigned int n = 0; n < _num_c * 2; ++n)
      A[m][n] = 0;
  }

  // fill in the non-zero elements in A
  for (unsigned int m = 0; m < _num_c; ++m)
  {
    for (unsigned int n = 0; n < _num_c; ++n)
    {
      // equal chemical potential derivative equations
      A[m * 2][n * 2] = (*_d2Fjdcjdbj[0][m][n])[_qp];
      A[m * 2][n * 2 + 1] = -(*_d2Fjdcjdbj[1][m][n])[_qp];
    }

    // concentration conservation derivative equations
    A[m * 2 + 1][m * 2] = 1 - _prop_h[_qp];
    A[m * 2 + 1][m * 2 + 1] = _prop_h[_qp];
  }

  MatrixTools::inverse(A, A);

  // loop through taking derivative wrt the ith component, they have the same A, but have
  // different k_c
  for (unsigned int i = 0; i < _num_c; ++i)
  {
    std::vector<Real> k_c(_num_c * 2);
    std::vector<Real> x_c(_num_c * 2);

    // assign non-zero elements in k_c
    k_c[i * 2 + 1] = 1;

    // compute x_c
    for (unsigned int m = 0; m < _num_c * 2; ++m)
    {
      for (unsigned int n = 0; n < _num_c * 2; ++n)
        x_c[m] += A[m][n] * k_c[n];
    }

    // assign the values in x_c to _prop_dcidb
    for (unsigned int m = 0; m < _num_c; ++m)
    {
      for (unsigned int n = 0; n < 2; ++n)
        (*_prop_dcidb[m][n][i])[_qp] = x_c[m * 2 + n];
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////// solve linear system of constraint derivatives wrt eta for computing dcideta
  // using the same linear matrix as computring dcidb
  std::vector<Real> k_eta(_num_c * 2);
  std::vector<Real> x_eta(_num_c * 2);

  // fill in k_eta
  for (unsigned int m = 0; m < _num_c; ++m)
  {
    k_eta[m * 2] = 0;
    k_eta[m * 2 + 1] = _prop_dh[_qp] * ((*_prop_ci[m * 2])[_qp] - (*_prop_ci[m * 2 + 1])[_qp]);
  }

  // compute x_eta
  for (unsigned int m = 0; m < _num_c * 2; ++m)
  {
    for (unsigned int n = 0; n < _num_c * 2; ++n)
      x_eta[m] += A[m][n] * k_eta[n];
  }

  // assign the values in x_eta to _prop_dcidetaj
  for (unsigned int m = 0; m < _num_c; ++m)
  {
    for (unsigned int n = 0; n < 2; ++n)
      (*_prop_dcideta[m][n])[_qp] = x_eta[m * 2 + n];
  }
}
