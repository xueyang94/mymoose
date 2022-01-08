//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KKSPhaseConcentrationMaterial.h"
#include "MatrixTools.h"

registerMooseObject("PhaseFieldApp", KKSPhaseConcentrationMaterial);

InputParameters
KKSPhaseConcentrationMaterial::validParams()
{
  InputParameters params = DerivativeMaterialInterface<Material>::validParams();
  params.addClassDescription(
      "Computes the KKS phase concentrations by using Newton iteration to solve the equal chemical "
      "potential and concentration conservation equations");
  params.addRequiredCoupledVar("global_cs", "The interpolated concentrations");
  params.addRequiredParam<MaterialPropertyName>("h_name", "Name of the switching function");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "ci_names", "Phase concentrations. The phase order must match Fi_names");
  params.addRequiredParam<std::vector<Real>>("ci_IC",
                                             "Initial values of ci in the same order of ci_names");
  params.addRequiredParam<MaterialName>("F1_material", "F1");
  params.addRequiredParam<MaterialName>("F2_material", "F2");
  params.addRequiredParam<std::vector<MaterialPropertyName>>("Fi_names", "Fi");
  params.addParam<MaterialPropertyName>(
      "nested_iterations", "The number of nested Newton iterations at each quadrature point");
  params.set<unsigned int>("min_iterations") = 10;
  params.set<unsigned int>("max_iterations") = 1000;
  params.set<Real>("absolute_tolerance") = 1e-13;
  params.set<Real>("relative_tolerance") = 1e-8;
  return params;
}

KKSPhaseConcentrationMaterial::KKSPhaseConcentrationMaterial(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _prop_c(coupledValues("global_cs")),
    _num_c(coupledComponents("global_cs")),
    _prop_h(getMaterialProperty<Real>("h_name")),
    _ci_names(getParam<std::vector<MaterialPropertyName>>("ci_names")),
    _prop_ci(_num_c * 2),
    _ci_old(_num_c * 2),
    _ci_IC(getParam<std::vector<Real>>("ci_IC")),
    _f1(getMaterial("F1_material")),
    _f2(getMaterial("F2_material")),
    _Fi_names(getParam<std::vector<MaterialPropertyName>>("Fi_names")),
    _first_dFi(2),
    _second_dFi(2),
    _iter(declareProperty<Real>("nested_iterations")),
    _abs_tol(getParam<Real>("absolute_tolerance")),
    _rel_tol(getParam<Real>("relative_tolerance")),
    _nested_solve(NestedSolve(parameters))

{
  for (unsigned int m = 0; m < _num_c * 2; ++m)
  {
    _prop_ci[m] = &declareProperty<Real>(_ci_names[m]);

    _ci_old[m] = &getMaterialPropertyOld<Real>(_ci_names[m]);
  }

  for (unsigned int m = 0; m < 2; ++m)
  {
    _first_dFi[m].resize(_num_c);
    _second_dFi[m].resize(_num_c);

    for (unsigned int n = 0; n < _num_c; ++n)
    {
      _first_dFi[m][n] = &getMaterialPropertyDerivative<Real>(_Fi_names[m], _ci_names[m + n * 2]);

      _second_dFi[m][n].resize(_num_c);

      for (unsigned int l = 0; l < _num_c; ++l)
        _second_dFi[m][n][l] = &getMaterialPropertyDerivative<Real>(
            _Fi_names[m], _ci_names[m + n * 2], _ci_names[m + l * 2]);
    }
  }
}

void
KKSPhaseConcentrationMaterial::initQpStatefulProperties()
{
  for (unsigned int m = 0; m < _num_c * 2; ++m)
    (*_prop_ci[m])[_qp] = _ci_IC[m];
}

void
KKSPhaseConcentrationMaterial::computeQpProperties()
{
  NestedSolve::Value<> solution(_num_c * 2); // dynamicaly sized vector class from the Eigen library

  solution << (*_ci_old[0])[_qp], (*_ci_old[1])[_qp], (*_ci_old[2])[_qp], (*_ci_old[3])[_qp];

  _nested_solve.setAbsoluteTolerance(_abs_tol);
  _nested_solve.setRelativeTolerance(_rel_tol);

  auto compute = [&](const NestedSolve::Value<> & guess,
                     NestedSolve::Value<> & residual,
                     NestedSolve::Jacobian<> & jacobian) {
    for (unsigned int m = 0; m < _num_c * 2; ++m)
      (*_prop_ci[m])[_qp] = guess(m);

    _f1.computePropertiesAtQp(_qp);
    _f2.computePropertiesAtQp(_qp);

    // assign residual functions
    for (unsigned int m = 0; m < _num_c; ++m)
    {
      residual(m * 2) = (*_first_dFi[0][m])[_qp] - (*_first_dFi[1][m])[_qp];
      residual(m * 2 + 1) = (1 - _prop_h[_qp]) * (*_prop_ci[m * 2])[_qp] +
                            _prop_h[_qp] * (*_prop_ci[m * 2 + 1])[_qp] - (*_prop_c[m])[_qp];
    }

    // initialize all elements in jacobian to be zero
    for (unsigned int m = 0; m < _num_c * 2; ++m)
    {
      for (unsigned int n = 0; n < _num_c * 2; ++n)
        jacobian(m, n) = 0;
    }

    // fill in the non-zero elements in jacobian
    for (unsigned int m = 0; m < _num_c; ++m)
    {
      for (unsigned int n = 0; n < _num_c; ++n)
      {
        // equal chemical potential derivative equations
        jacobian(m * 2, n * 2) = (*_second_dFi[0][m][n])[_qp];
        jacobian(m * 2, n * 2 + 1) = -(*_second_dFi[1][m][n])[_qp];
      }

      // concentration conservation derivative equations
      jacobian(m * 2 + 1, m * 2) = 1 - _prop_h[_qp];
      jacobian(m * 2 + 1, m * 2 + 1) = _prop_h[_qp];
    }
  };

  _nested_solve.nonlinear(solution, compute);
  _iter[_qp] = _nested_solve.getIterations();

  if (_nested_solve.getState() == NestedSolve::State::NOT_CONVERGED)
    std::cout << "Newton iteration did not converge." << std::endl;

  // assign solution to ci
  for (unsigned int m = 0; m < _num_c * 2; ++m)
    (*_prop_ci[m])[_qp] = solution[m];
}
