//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SubConcentration.h"
// #include "libmesh/utility.h"
// #include <cmath>
// #include "libmesh/vector_value.h"
// #include "RankTwoTensor.h"
// #include "gtest/gtest.h"
// #include "Conversion.h"
// #include "IndirectSort.h"

registerMooseObject("PhaseFieldApp", SubConcentration);

InputParameters
SubConcentration::validParams()
{
  InputParameters params = DerivativeMaterialInterface<Material>::validParams();
  params.addClassDescription(
      "Computes the KKS sub-concentrations by using Newton iteration to solve the equal chemical "
      "potential and concentration conservation equations");
  params.addRequiredCoupledVar("global_c", "The interpolated concentration");
  params.addRequiredCoupledVar("all_etas", "Vector of all order parameters for all phases");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "hj_names", "Names of the switching functions in the same order of the all_etas");
  params.addParam<std::vector<MaterialPropertyName>>("ci_names", "Phase concentrations");
  params.addRequiredParam<std::vector<Real>>("ci_IC",
                                             "Initial values of ci in the same order of ci_names");

  params.addRequiredParam<MaterialPropertyName>("dc1dc_name",
                                                "The first derivative of c1 w.r.t. c");
  params.addRequiredParam<MaterialPropertyName>("dc2dc_name",
                                                "The first derivative of c2 w.r.t. c");
  params.addRequiredParam<MaterialPropertyName>("dc3dc_name",
                                                "The first derivative of c3 w.r.t. c");

  params.addRequiredParam<MaterialPropertyName>("dc1deta1_name",
                                                "The first derivative of c1 w.r.t. eta1");
  params.addRequiredParam<MaterialPropertyName>("dc1deta2_name",
                                                "The first derivative of c1 w.r.t. eta2");
  params.addRequiredParam<MaterialPropertyName>("dc1deta3_name",
                                                "The first derivative of c1 w.r.t. eta3");
  params.addRequiredParam<MaterialPropertyName>("dc2deta1_name",
                                                "The first derivative of c2 w.r.t. eta1");
  params.addRequiredParam<MaterialPropertyName>("dc2deta2_name",
                                                "The first derivative of c2 w.r.t. eta2");
  params.addRequiredParam<MaterialPropertyName>("dc2deta3_name",
                                                "The first derivative of c2 w.r.t. eta3");
  params.addRequiredParam<MaterialPropertyName>("dc3deta1_name",
                                                "The first derivative of c3 w.r.t. eta1");
  params.addRequiredParam<MaterialPropertyName>("dc3deta2_name",
                                                "The first derivative of c3 w.r.t. eta2");
  params.addRequiredParam<MaterialPropertyName>("dc3deta3_name",
                                                "The first derivative of c3 w.r.t. eta3");

  params.addRequiredParam<MaterialName>("F1_material", "F1");
  params.addRequiredParam<MaterialName>("F2_material", "F2");
  params.addRequiredParam<MaterialName>("F3_material", "F3");
  params.addRequiredParam<MaterialPropertyName>("F1_name", "F1");
  params.addRequiredParam<MaterialPropertyName>("F2_name", "F2");
  params.addRequiredParam<MaterialPropertyName>("F3_name", "F3");
  // params.addParam<Real>("absolute_tol_value", 1e-9, "Absolute tolerance of the Newton
  // iteration"); params.addParam<Real>("relative_tol_value", 1e-9, "Relative tolerance of the
  // Newton iteration"); params.addParam<Real>("max_iteration", 100, "The maximum number of Newton
  // iterations");

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

    _ci_IC(getParam<std::vector<Real>>("ci_IC")),

    _dc1dc(declareProperty<Real>("dc1dc_name")),
    _dc2dc(declareProperty<Real>("dc2dc_name")),
    _dc3dc(declareProperty<Real>("dc3dc_name")),

    _dc1deta1(declareProperty<Real>("dc1deta1_name")),
    _dc1deta2(declareProperty<Real>("dc1deta2_name")),
    _dc1deta3(declareProperty<Real>("dc1deta3_name")),
    _dc2deta1(declareProperty<Real>("dc2deta1_name")),
    _dc2deta2(declareProperty<Real>("dc2deta2_name")),
    _dc2deta3(declareProperty<Real>("dc2deta3_name")),
    _dc3deta1(declareProperty<Real>("dc3deta1_name")),
    _dc3deta2(declareProperty<Real>("dc3deta2_name")),
    _dc3deta3(declareProperty<Real>("dc3deta3_name")),

    _f1(getMaterial("F1_material")),
    _f2(getMaterial("F2_material")),
    _f3(getMaterial("F3_material")),

    _c1_name("c1"),
    _c2_name("c2"),
    _c3_name("c3"),
    _first_df1(getMaterialPropertyDerivative<Real>("F1_name", _c1_name)),
    _first_df2(getMaterialPropertyDerivative<Real>("F2_name", _c2_name)),
    _first_df3(getMaterialPropertyDerivative<Real>("F3_name", _c3_name)),
    _second_df1(getMaterialPropertyDerivative<Real>("F1_name", _c1_name, _c1_name)),
    _second_df2(getMaterialPropertyDerivative<Real>("F2_name", _c2_name, _c2_name)),
    _second_df3(getMaterialPropertyDerivative<Real>("F3_name", _c3_name, _c3_name)),

    _iter(declareProperty<Real>("nested_iterations")),
    _abs_tol(getParam<Real>("absolute_tolerance")),
    _rel_tol(getParam<Real>("relative_tolerance")),
    _nested_solve(NestedSolve(parameters))

{
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
    // declare h and dh material properties
    _prop_hj[m] = &getMaterialPropertyByName<Real>(_hj_names[m]);

    // _prop_dhjdetai[m][n] is the derivative of h_j w.r.t. eta_i
    _prop_dhjdetai[m].resize(_num_eta);

    for (unsigned int n = 0; n < _num_eta; ++n)
      _prop_dhjdetai[m][n] = &getMaterialPropertyDerivative<Real>(_hj_names[m], _eta_names[n]);
  }
}

// This function is also defined in NestedSolve.C (protected)
void
linear(RankTwoTensor A, RealVectorValue & x, RealVectorValue b)
{
  x = A.inverse() * b;
}

void
SubConcentration::computeQpProperties()
{
  // NestedSolve solver;
  NestedSolve::Value<> solution(3); // dynamicly sized vector class from the Eigen library
  solution << _ci_IC[0], _ci_IC[1], _ci_IC[2];
  // solver.setAbsoluteTolerance(_abs_tol);
  // solver.setRelativeTolerance(_rel_tol);
  _nested_solve.setAbsoluteTolerance(_abs_tol);
  _nested_solve.setRelativeTolerance(_rel_tol);

  auto compute = [&](const NestedSolve::Value<> & guess,
                     NestedSolve::Value<> & residual,
                     NestedSolve::Jacobian<> & jacobian) {
    (*_ci_prop[0])[_qp] = guess(0);
    (*_ci_prop[1])[_qp] = guess(1);
    (*_ci_prop[2])[_qp] = guess(2);
    _f1.computePropertiesAtQp(_qp);
    _f2.computePropertiesAtQp(_qp);
    _f3.computePropertiesAtQp(_qp);

    residual(0) = _first_df1[_qp] - _first_df2[_qp];
    residual(1) = _first_df2[_qp] - _first_df3[_qp];
    residual(2) = (*_prop_hj[0])[_qp] * guess(0) + (*_prop_hj[1])[_qp] * guess(1) +
                  (*_prop_hj[2])[_qp] * guess(2) - _c[_qp];

    jacobian(0, 0) = _second_df1[_qp];
    jacobian(0, 1) = -_second_df2[_qp];
    jacobian(0, 2) = 0;
    jacobian(1, 0) = 0;
    jacobian(1, 1) = _second_df2[_qp];
    jacobian(1, 2) = -_second_df3[_qp];
    jacobian(2, 0) = (*_prop_hj[0])[_qp];
    jacobian(2, 1) = (*_prop_hj[1])[_qp];
    jacobian(2, 2) = (*_prop_hj[2])[_qp];
  };

  // solver.nonlinear(solution, compute);
  _nested_solve.nonlinear(solution, compute);
  _iter[_qp] = _nested_solve.getIterations();

  // if (solver.getState() == NestedSolve::State::NOT_CONVERGED)
  if (_nested_solve.getState() == NestedSolve::State::NOT_CONVERGED)
  {
    std::cout << "Newton iteration did not converge." << std::endl;
  }

  (*_ci_prop[0])[_qp] = solution[0];
  (*_ci_prop[1])[_qp] = solution[1];
  (*_ci_prop[2])[_qp] = solution[2];
  _f1.computePropertiesAtQp(_qp);
  _f2.computePropertiesAtQp(_qp);
  _f3.computePropertiesAtQp(_qp);

  ////////////////////////////////////////////////////////////////////////////////////////// compute dc1dc, dc2dc, and dc3dc
  // The matrix A used to compute dcidc and dcidetai are the same as the jacobian matrix
  RankTwoTensor A;
  RealVectorValue x_dcidc;
  RealVectorValue b_dcidc{0, 0, 1};

  A(0, 0) = _second_df1[_qp];
  A(0, 1) = -_second_df2[_qp];
  A(0, 2) = 0;
  A(1, 0) = 0;
  A(1, 1) = _second_df2[_qp];
  A(1, 2) = -_second_df3[_qp];
  A(2, 0) = (*_prop_hj[0])[_qp];
  A(2, 1) = (*_prop_hj[1])[_qp];
  A(2, 2) = (*_prop_hj[2])[_qp];

  linear(A, x_dcidc, b_dcidc);
  _dc1dc[_qp] = x_dcidc(0);
  _dc2dc[_qp] = x_dcidc(1);
  _dc3dc[_qp] = x_dcidc(2);

  //////////////////////////////////////////////////////////////////////////////////////////// compute dc1deta1, dc2deta1, and dc3deta1
  RealVectorValue x_cideta1;
  RealVectorValue b_cideta1{0,
                            0,
                            -(*_prop_dhjdetai[0][0])[_qp] * (*_ci_prop[0])[_qp] -
                                (*_prop_dhjdetai[1][0])[_qp] * (*_ci_prop[1])[_qp] -
                                (*_prop_dhjdetai[2][0])[_qp] * (*_ci_prop[2])[_qp]};

  linear(A, x_cideta1, b_cideta1);
  _dc1deta1[_qp] = x_cideta1(0);
  _dc2deta1[_qp] = x_cideta1(1);
  _dc3deta1[_qp] = x_cideta1(2);

  //////////////////////////////////////////////////////////////////////////////////////////// compute dc1deta2, dc2deta2, and dc3deta2
  RealVectorValue x_cideta2;
  RealVectorValue b_cideta2{0,
                            0,
                            -(*_prop_dhjdetai[0][1])[_qp] * (*_ci_prop[0])[_qp] -
                                (*_prop_dhjdetai[1][1])[_qp] * (*_ci_prop[1])[_qp] -
                                (*_prop_dhjdetai[2][1])[_qp] * (*_ci_prop[2])[_qp]};

  linear(A, x_cideta2, b_cideta2);
  _dc1deta2[_qp] = x_cideta2(0);
  _dc2deta2[_qp] = x_cideta2(1);
  _dc3deta2[_qp] = x_cideta2(2);

  //////////////////////////////////////////////////////////////////////////////////////////// compute dc1deta3, dc2deta3, and dc3deta3
  RealVectorValue x_cideta3;
  RealVectorValue b_cideta3{0,
                            0,
                            -(*_prop_dhjdetai[0][2])[_qp] * (*_ci_prop[0])[_qp] -
                                (*_prop_dhjdetai[1][2])[_qp] * (*_ci_prop[1])[_qp] -
                                (*_prop_dhjdetai[2][2])[_qp] * (*_ci_prop[2])[_qp]};

  linear(A, x_cideta3, b_cideta3);
  _dc1deta3[_qp] = x_cideta3(0);
  _dc2deta3[_qp] = x_cideta3(1);
  _dc3deta3[_qp] = x_cideta3(2);
}
