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

  params.addRequiredCoupledVar("all_etas", "Vector of all order parameters for all phases");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "hj_names", "Names of the switching functions in the same order of the all_etas");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "ci_names", "Phase concentrations. The phase order must match Fi_names");

  params.addRequiredParam<MaterialPropertyName>("c1_name", "c1 name");
  params.addRequiredParam<MaterialPropertyName>("c2_name", "c2 name");
  params.addRequiredParam<MaterialPropertyName>("c3_name", "c3 name");
  params.addRequiredParam<std::vector<Real>>("ci_IC",
                                             "Initial values of ci in the same order of ci_names");

  params.addRequiredParam<MaterialName>("F1_material", "F1");
  params.addRequiredParam<MaterialName>("F2_material", "F2");
  params.addRequiredParam<MaterialName>("F3_material", "F3");

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
    _num_eta(coupledComponents("all_etas")),
    _eta_names(coupledNames("all_etas")),
    _hj_names(getParam<std::vector<MaterialPropertyName>>("hj_names")),
    _prop_hj(_num_eta),
    _prop_dhjdetai(_num_eta),

    _ci_names(getParam<std::vector<MaterialPropertyName>>("ci_names")),
    _ci_prop(_num_eta),

    _c1_old(getMaterialPropertyOld<Real>("c1_name")), // old
    _c2_old(getMaterialPropertyOld<Real>("c2_name")), // old
    _c3_old(getMaterialPropertyOld<Real>("c3_name")), // old
    _ci_IC(getParam<std::vector<Real>>("ci_IC")),

    // _Fi_material(getParam<std::vector<MaterialName>>("Fi_material_names")),
    // _fi_material(_num_eta),
    // _f1(getMaterial(_Fi_material[0])),
    // _f2(getMaterial(_Fi_material[1])),
    // _f3(getMaterial(_Fi_material[2])),

    _f1(getMaterial("F1_material")),
    _f2(getMaterial("F2_material")),
    _f3(getMaterial("F3_material")),

    _Fi_names(getParam<std::vector<MaterialPropertyName>>("Fi_names")),
    _first_dFi(_num_eta),
    _second_dFi(_num_eta),

    _iter(declareProperty<Real>("nested_iterations")),
    _abs_tol(getParam<Real>("absolute_tolerance")),
    _rel_tol(getParam<Real>("relative_tolerance")),
    _nested_solve(NestedSolve(parameters))

{
  // declare the first and second derivative of phase energy wrt phase concentrations
  for (unsigned int n = 0; n < _num_eta; ++n)
  {
    _first_dFi[n] = &getMaterialPropertyDerivative<Real>(_Fi_names[n], _ci_names[n]);
    _second_dFi[n] = &getMaterialPropertyDerivative<Real>(_Fi_names[n], _ci_names[n], _ci_names[n]);
  }

  // error check parameters
  if (_hj_names.size() != 0 && _hj_names.size() != _num_eta)
    paramError("hj_names",
               "Specify either as many entries are eta values or none at all for auto-naming the "
               "hj material properties.");

  if (_ci_names.size() != 0 && _ci_names.size() != _num_eta)
    paramError("ci_names",
               "Specify either as many entries are eta values or none at all for auto-naming the "
               "ci material properties.");

  // declare ci material properties
  for (unsigned int i = 0; i < _num_eta; ++i)
    _ci_prop[i] = &declareProperty<Real>(_ci_names[i]);

  for (unsigned int m = 0; m < _num_eta; ++m)
  {
    // _ci_old_prop[m] = &declareProperty<Real>(_ci_old[i]);

    // declare h and dh material properties
    _prop_hj[m] = &getMaterialPropertyByName<Real>(_hj_names[m]);

    // _prop_dhjdetai[m][n] is the derivative of h_j w.r.t. eta_i
    _prop_dhjdetai[m].resize(_num_eta);
    for (unsigned int n = 0; n < _num_eta; ++n)
      _prop_dhjdetai[m][n] = &getMaterialPropertyDerivative<Real>(_hj_names[m], _eta_names[n]);
  }
}

void
SubConcentration::initQpStatefulProperties()
{
  for (unsigned int i = 0; i < _num_eta; ++i)
    (*_ci_prop[i])[_qp] = _ci_IC[i];
}

// // This function is also defined in NestedSolve.C (protected)
// // And the inverse() method does not give the inverse of A, strange
// void
// linear(RankTwoTensor A, RealVectorValue & x, RealVectorValue b)
// {
//   x = A.inverse() * b;
// }

void
SubConcentration::computeQpProperties()
{
  NestedSolve::Value<> solution(3); // dynamicly sized vector class from the Eigen library
  // solution << _ci_IC[0], _ci_IC[1], _ci_IC[2];
  solution << _c1_old[_qp], _c2_old[_qp], _c3_old[_qp];
  // solution << (*_ci_old_prop[0])[_qp], (*_ci_old_prop[1])[_qp], (*_ci_old_prop[2])[_qp];
  // solution << _ci_old[_qp][0], _ci_old[_qp][1], _ci_old[_qp][2];

  _nested_solve.setAbsoluteTolerance(_abs_tol);
  _nested_solve.setRelativeTolerance(_rel_tol);

  auto compute = [&](const NestedSolve::Value<> & guess,
                     NestedSolve::Value<> & residual,
                     NestedSolve::Jacobian<> & jacobian) {
    for (unsigned int i = 0; i < _num_eta; ++i)
      (*_ci_prop[i])[_qp] = guess(i);

    _f1.computePropertiesAtQp(_qp);
    _f2.computePropertiesAtQp(_qp);
    _f3.computePropertiesAtQp(_qp);
    // _fi_material[0]->computePropertiesAtQp(_qp);
    // _fi_material[1]->computePropertiesAtQp(_qp);
    // _fi_material[2]->computePropertiesAtQp(_qp);

    residual(0) = (*_first_dFi[0])[_qp] - (*_first_dFi[1])[_qp];
    residual(1) = (*_first_dFi[1])[_qp] - (*_first_dFi[2])[_qp];
    residual(2) = (*_prop_hj[0])[_qp] * guess(0) + (*_prop_hj[1])[_qp] * guess(1) +
                  (*_prop_hj[2])[_qp] * guess(2) - _c[_qp];

    jacobian(0, 0) = (*_second_dFi[0])[_qp];
    jacobian(0, 1) = -(*_second_dFi[1])[_qp];
    jacobian(0, 2) = 0;
    jacobian(1, 0) = 0;
    jacobian(1, 1) = (*_second_dFi[1])[_qp];
    jacobian(1, 2) = -(*_second_dFi[2])[_qp];
    jacobian(2, 0) = (*_prop_hj[0])[_qp];
    jacobian(2, 1) = (*_prop_hj[1])[_qp];
    jacobian(2, 2) = (*_prop_hj[2])[_qp];
  };

  _nested_solve.nonlinear(solution, compute);
  _iter[_qp] = _nested_solve.getIterations();

  if (_nested_solve.getState() == NestedSolve::State::NOT_CONVERGED)
  {
    std::cout << "Newton iteration did not converge." << std::endl;
  }

  for (unsigned int i = 0; i < _num_eta; ++i)
    (*_ci_prop[i])[_qp] = solution[i];
}
