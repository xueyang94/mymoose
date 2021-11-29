//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SubConcentration.h"
#include "MatrixTools.h"

registerMooseObject("PhaseFieldApp", SubConcentration);

InputParameters
SubConcentration::validParams()
{
  InputParameters params = DerivativeMaterialInterface<Material>::validParams();
  params.addClassDescription(
      "Computes the KKS phase concentrations by using Newton iteration to solve the equal chemical "
      "potential and concentration conservation equations");
  params.addRequiredCoupledVar("global_c", "The interpolated concentration");
  params.addRequiredParam<MaterialPropertyName>("h_name", "Name of the switching function");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "ci_names", "Phase concentrations. The phase order must match Fi_names");

  params.addRequiredParam<MaterialPropertyName>("c1_name", "c1 name");
  params.addRequiredParam<MaterialPropertyName>("c2_name", "c2 name");
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

SubConcentration::SubConcentration(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _c(coupledValue("global_c")),
    _prop_h(getMaterialProperty<Real>("h_name")),
    _ci_names(getParam<std::vector<MaterialPropertyName>>("ci_names")),
    _prop_ci(2),
    _c1_old(getMaterialPropertyOld<Real>("c1_name")), // old
    _c2_old(getMaterialPropertyOld<Real>("c2_name")), // old
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
  // declare the first and second derivative of phase energy wrt phase concentrations
  for (unsigned int n = 0; n < 2; ++n)
  {
    _first_dFi[n] = &getMaterialPropertyDerivative<Real>(_Fi_names[n], _ci_names[n]);
    _second_dFi[n] = &getMaterialPropertyDerivative<Real>(_Fi_names[n], _ci_names[n], _ci_names[n]);
  }

  // declare ci material properties
  for (unsigned int i = 0; i < 2; ++i)
    _prop_ci[i] = &declareProperty<Real>(_ci_names[i]);
}

void
SubConcentration::initQpStatefulProperties()
{
  for (unsigned int i = 0; i < 2; ++i)
    (*_prop_ci[i])[_qp] = _ci_IC[i];
}

void
SubConcentration::computeQpProperties()
{
  NestedSolve::Value<> solution(2); // dynamicly sized vector class from the Eigen library
  solution << _c1_old[_qp], _c2_old[_qp];

  _nested_solve.setAbsoluteTolerance(_abs_tol);
  _nested_solve.setRelativeTolerance(_rel_tol);

  auto compute = [&](const NestedSolve::Value<> & guess,
                     NestedSolve::Value<> & residual,
                     NestedSolve::Jacobian<> & jacobian) {
    for (unsigned int i = 0; i < 2; ++i)
      (*_prop_ci[i])[_qp] = guess(i);

    _f1.computePropertiesAtQp(_qp);
    _f2.computePropertiesAtQp(_qp);

    residual(0) = (*_first_dFi[0])[_qp] - (*_first_dFi[1])[_qp];
    residual(1) =
        (1 - _prop_h[_qp]) * (*_prop_ci[0])[_qp] + _prop_h[_qp] * (*_prop_ci[1])[_qp] - _c[_qp];

    jacobian(0, 0) = (*_second_dFi[0])[_qp];
    jacobian(0, 1) = -(*_second_dFi[1])[_qp];
    jacobian(1, 0) = 1 - _prop_h[_qp];
    jacobian(1, 1) = _prop_h[_qp];
  };

  _nested_solve.nonlinear(solution, compute);
  _iter[_qp] = _nested_solve.getIterations();

  if (_nested_solve.getState() == NestedSolve::State::NOT_CONVERGED)
    std::cout << "Newton iteration did not converge." << std::endl;

  for (unsigned int i = 0; i < 2; ++i)
    (*_prop_ci[i])[_qp] = solution[i];
}
