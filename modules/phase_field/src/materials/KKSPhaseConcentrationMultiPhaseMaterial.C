//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KKSPhaseConcentrationMultiPhaseMaterial.h"
#include "MatrixTools.h"
#include <Eigen/Dense>
#include <Eigen/Core>

registerMooseObject("PhaseFieldApp", KKSPhaseConcentrationMultiPhaseMaterial);

InputParameters
KKSPhaseConcentrationMultiPhaseMaterial::validParams()
{
  InputParameters params = DerivativeMaterialInterface<Material>::validParams();
  params.addClassDescription(
      "Computes the KKS phase concentrations by using a nested Newton iteration "
      "to solve the equal chemical potential and concentration conservation equations");
  params.addRequiredCoupledVar("global_cs", "Global concentrations, for example, c, b.");
  params.addRequiredCoupledVar("all_etas", "Order parameters for all phases.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "hj_names", "Names of the switching functions in the same order of the all_etas");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "ci_names",
      "Phase concentrations. These must have the same order as Fj_names and global_cs, for "
      "example, c1, c2, b1, b2.");
  params.addRequiredParam<std::vector<Real>>("ci_IC",
                                             "Initial values of ci in the same order of ci_names");

  params.addRequiredParam<MaterialName>("F1_material", "F1");
  params.addRequiredParam<MaterialName>("F2_material", "F2");
  params.addRequiredParam<MaterialName>("F3_material", "F3");
  // params.addRequiredParam<MaterialName>("F4_material", "F4");
  // params.addRequiredParam<MaterialName>("F5_material", "F5");

  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "Fi_names", "Phase energies in the same order as all_etas");
  params.addParam<MaterialPropertyName>(
      "nested_iterations", "The number of nested Newton iterations at each quadrature point");
  params.set<unsigned int>("min_iterations") = 10;
  params.set<unsigned int>("max_iterations") = 1000;
  params.set<Real>("absolute_tolerance") = 1e-13;
  params.set<Real>("relative_tolerance") = 1e-8;
  return params;
}

KKSPhaseConcentrationMultiPhaseMaterial::KKSPhaseConcentrationMultiPhaseMaterial(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _prop_c(coupledValues("global_cs")),
    _num_c(coupledComponents("global_cs")),
    _num_j(coupledComponents("all_etas")),
    _eta_names(coupledNames("all_etas")),
    _hj_names(getParam<std::vector<MaterialPropertyName>>("hj_names")),
    _prop_hj(_num_j),
    _ci_names(getParam<std::vector<MaterialPropertyName>>("ci_names")),
    _prop_ci(_num_c * _num_j),
    _ci_old(_num_c * _num_j),
    _ci_IC(getParam<std::vector<Real>>("ci_IC")),

    _Fi_names(getParam<std::vector<MaterialPropertyName>>("Fi_names")),
    _first_dFi(_num_j),
    _second_dFi(_num_j),

    _iter(declareProperty<Real>("nested_iterations")),
    _abs_tol(getParam<Real>("absolute_tolerance")),
    _rel_tol(getParam<Real>("relative_tolerance")),
    _nested_solve(NestedSolve(parameters))

{
  for (unsigned int m = 0; m < _num_c * _num_j; ++m)
  {
    _ci_old[m] = &getMaterialPropertyOld<Real>(_ci_names[m]);
    _prop_ci[m] = &declareProperty<Real>(_ci_names[m]);
  }

  for (unsigned int m = 0; m < _num_j; ++m)
  {
    _first_dFi[m].resize(_num_c);
    _second_dFi[m].resize(_num_c);

    _prop_hj[m] = &getMaterialPropertyByName<Real>(_hj_names[m]);

    for (unsigned int n = 0; n < _num_c; ++n)
    {
      _first_dFi[m][n] =
          &getMaterialPropertyDerivative<Real>(_Fi_names[m], _ci_names[m + n * _num_j]);

      _second_dFi[m][n].resize(_num_c);

      for (unsigned int l = 0; l < _num_c; ++l)
        _second_dFi[m][n][l] = &getMaterialPropertyDerivative<Real>(
            _Fi_names[m], _ci_names[m + n * _num_j], _ci_names[m + l * _num_j]);
    }
  }
}

void
KKSPhaseConcentrationMultiPhaseMaterial::initQpStatefulProperties()
{
  for (unsigned int m = 0; m < _num_c * _num_j; ++m)
    (*_prop_ci[m])[_qp] = _ci_IC[m];
}

void
KKSPhaseConcentrationMultiPhaseMaterial::initialSetup()
{
  _F1 = &getMaterial("F1_material");
  _F2 = &getMaterial("F2_material");
  _F3 = &getMaterial("F3_material");
}

void
KKSPhaseConcentrationMultiPhaseMaterial::computeQpProperties()
{
  // dynamicaly sized vector class from the Eigen library
  NestedSolve::Value<> solution(_num_c * _num_j);

  ////////////////////////////////////// need change start
  // solution << (*_ci_old[0])[_qp], (*_ci_old[1])[_qp], (*_ci_old[2])[_qp];

  solution << (*_ci_old[0])[_qp], (*_ci_old[1])[_qp], (*_ci_old[2])[_qp], (*_ci_old[3])[_qp],
      (*_ci_old[4])[_qp], (*_ci_old[5])[_qp];

  // solution << (*_ci_old[0])[_qp], (*_ci_old[1])[_qp], (*_ci_old[2])[_qp], (*_ci_old[3])[_qp],
  //     (*_ci_old[4])[_qp], (*_ci_old[5])[_qp], (*_ci_old[6])[_qp], (*_ci_old[7])[_qp],
  //     (*_ci_old[8])[_qp], (*_ci_old[9])[_qp], (*_ci_old[10])[_qp], (*_ci_old[11])[_qp];

  // solution << (*_ci_old[0])[_qp], (*_ci_old[1])[_qp], (*_ci_old[2])[_qp], (*_ci_old[3])[_qp],
  //     (*_ci_old[4])[_qp], (*_ci_old[5])[_qp], (*_ci_old[6])[_qp], (*_ci_old[7])[_qp],
  //     (*_ci_old[8])[_qp], (*_ci_old[9])[_qp], (*_ci_old[10])[_qp], (*_ci_old[11])[_qp],
  //     (*_ci_old[12])[_qp], (*_ci_old[13])[_qp], (*_ci_old[14])[_qp], (*_ci_old[15])[_qp],
  //     (*_ci_old[16])[_qp], (*_ci_old[17])[_qp], (*_ci_old[18])[_qp], (*_ci_old[19])[_qp];
  ////////////////////////////////////// need change end

  _nested_solve.setAbsoluteTolerance(_abs_tol);
  _nested_solve.setRelativeTolerance(_rel_tol);

  auto compute = [&](const NestedSolve::Value<> & guess,
                     NestedSolve::Value<> & residual,
                     NestedSolve::Jacobian<> & jacobian) {
    for (unsigned int m = 0; m < _num_c * _num_j; ++m)
      (*_prop_ci[m])[_qp] = guess(m);

    ////////////////////////////////////// need change start
    _F1->computePropertiesAtQp(_qp);
    _F2->computePropertiesAtQp(_qp);
    _F3->computePropertiesAtQp(_qp);
    ////////////////////////////////////// need change end

    // assign residual functions
    for (unsigned int m = 0; m < _num_c; ++m)
    {
      for (unsigned int n = 0; n < _num_j - 1; ++n)
        residual(m * _num_j + n) = (*_first_dFi[n][m])[_qp] - (*_first_dFi[n + 1][m])[_qp];

      residual((m + 1) * _num_j - 1) = -(*_prop_c[m])[_qp];

      for (unsigned int l = 0; l < _num_j; ++l)
        residual((m + 1) * _num_j - 1) += (*_prop_hj[l])[_qp] * (*_prop_ci[m * _num_j + l])[_qp];
    }

    // initialize all terms in jacobian to be zero
    for (unsigned int m = 0; m < _num_j * _num_c; ++m)
    {
      for (unsigned int n = 0; n < _num_j * _num_c; ++n)
        jacobian(m, n) = 0;
    }

    // fill in the non-zero terms in jacobian
    // first assign the terms in jacobian that come from the mu equality derivative equations
    // loop through the constraint equation sets of the mth component
    for (unsigned int m = 0; m < _num_c; ++m)
    {
      // loop through the nth constraint equation in the constraint set of one component
      for (unsigned int n = 0; n < (_num_j - 1); ++n)
      {
        // loop through the lth chain rule terms in one constrain equation
        for (unsigned int l = 0; l < _num_c; ++l)
        {
          jacobian(m * _num_j + n, n + l * _num_j) = (*_second_dFi[n][m][l])[_qp];
          jacobian(m * _num_j + n, n + l * _num_j + 1) = -(*_second_dFi[n + 1][m][l])[_qp];
        }
      }
    }

    // then assign the terms in jacobian that come from the concentration conservation equations
    for (unsigned int m = 0; m < _num_c; ++m)
    {
      for (unsigned int n = 0; n < _num_j; ++n)
        jacobian((m + 1) * _num_j - 1, m * _num_j + n) = (*_prop_hj[n])[_qp];
    }
  };

  _nested_solve.nonlinear(solution, compute);
  _iter[_qp] = _nested_solve.getIterations();

  if (_nested_solve.getState() == NestedSolve::State::NOT_CONVERGED)
  {
    std::cout << "Newton iteration did not converge." << std::endl;
  }

  // assign solution to ci
  for (unsigned int m = 0; m < _num_c * _num_j; ++m)
    (*_prop_ci[m])[_qp] = solution[m];
}
