//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SubConcentration.h"
#include "libmesh/fparser_ad.hh"
#include "libmesh/utility.h"
#include <cmath>
#include "NestedSolve.h"
#include "libmesh/vector_value.h"
#include "RankTwoTensor.h"
#include "gtest/gtest.h"

registerMooseObject("PhaseFieldApp", SubConcentration);

InputParameters
SubConcentration::validParams()
{
  InputParameters params = DerivativeMaterialInterface<Material>::validParams();
  params.addClassDescription(
      "Computes the KKS sub-concentrations by using Newton iteration to solve the equal chemical "
      "potential and concentration conservation equations");
  params.addRequiredCoupledVar("global_c", "The interpolated concentration");
  params.addRequiredCoupledVar("eta", "The order parameter");
  params.addRequiredParam<MaterialPropertyName>("h_name", "Name of the switching function");
  params.addRequiredParam<MaterialPropertyName>("c1_name", "Name of the first phase concentration");
  params.addRequiredParam<MaterialPropertyName>("c2_name",
                                                "Name of the second phase concentration");
  params.addRequiredParam<Real>("c1_IC", "The initial value of c1");
  params.addRequiredParam<Real>("c2_IC", "The initial value of c2");
  params.addParam<Real>("absolute_tol_value", 1e-9, "Absolute tolerance of the Newton iteration");
  params.addParam<Real>("relative_tol_value", 1e-9, "Relative tolerance of the Newton iteration");
  params.addParam<Real>("max_iteration", 100, "The maximum number of Newton iterations");
  params.addRequiredParam<MaterialPropertyName>("dc1dc_name",
                                                "The first derivative of c1 w.r.t. c");
  params.addRequiredParam<MaterialPropertyName>("dc2dc_name",
                                                "The first derivative of c2 w.r.t. c");
  params.addRequiredParam<MaterialPropertyName>("dc1deta_name",
                                                "The first derivative of c1 w.r.t. eta");
  params.addRequiredParam<MaterialPropertyName>("dc2deta_name",
                                                "The first derivative of c2 w.r.t. eta");
  params.addRequiredParam<MaterialName>("F1_material", "F1");
  params.addRequiredParam<MaterialName>("F2_material", "F2");
  params.addRequiredParam<MaterialPropertyName>("F1_name", "F1");
  params.addRequiredParam<MaterialPropertyName>("F2_name", "F2");
  // params.addParam<Real>("num_iter", "The number of nested Newton iteration");
  return params;
}

SubConcentration::SubConcentration(const InputParameters & parameters)
  // : Material(parameters),
  : DerivativeMaterialInterface<Material>(parameters),
    _c(coupledValue("global_c")),
    _eta(coupledValue("eta")),
    _h(getMaterialProperty<Real>("h_name")),
    _c1(declareProperty<Real>("c1_name")),
    _c2(declareProperty<Real>("c2_name")),
    // _c1_old(getMaterialPropertyOld<Real>("c1_name")), // old
    // _c2_old(getMaterialPropertyOld<Real>("c2_name")), // old
    _c1_initial(getParam<Real>("c1_IC")),
    _c2_initial(getParam<Real>("c2_IC")),
    _abs_tol(getParam<Real>("absolute_tol_value")),
    _rel_tol(getParam<Real>("relative_tol_value")),
    _maxiter(getParam<Real>("max_iteration")),
    _dc1dc(declareProperty<Real>("dc1dc_name")),
    _dc2dc(declareProperty<Real>("dc2dc_name")),
    _dc1deta(declareProperty<Real>("dc1deta_name")),
    _dc2deta(declareProperty<Real>("dc2deta_name")),
    _c1_name("c1"),
    _c2_name("c2"),
    _f1(getMaterial("F1_material")),
    _f2(getMaterial("F2_material")),
    _first_df1(getMaterialPropertyDerivative<Real>("F1_name", _c1_name)),
    _first_df2(getMaterialPropertyDerivative<Real>("F2_name", _c2_name)),
    _second_df1(getMaterialPropertyDerivative<Real>("F1_name", _c1_name, _c1_name)),
    _second_df2(getMaterialPropertyDerivative<Real>("F2_name", _c2_name, _c2_name))
// _num_iter(declareProperty<Real>("num_iter"))

{
}

// void
// SubConcentration::initQpStatefulProperties()
// {
//   _c1[_qp] = _c1_initial;
//   _c2[_qp] = _c2_initial;
// }

void
SubConcentration::computeQpProperties()
{
  Real n = _eta[_qp];

  NestedSolve solver;
  // NestedSolve::Value<2> solution{_c1_initial, _c2_initial};
  NestedSolve::Value<> solution(2); // dynamicly sized vector class from the Eigen library
  solution << _c1_initial, _c2_initial;
  solver.setRelativeTolerance(1e-9);

  auto compute = [&](const NestedSolve::Value<> & guess,
                     NestedSolve::Value<> & residual,
                     NestedSolve::Jacobian<> & jacobian) {
    _c1[_qp] = guess(0);
    _c2[_qp] = guess(1);
    _f1.computePropertiesAtQp(_qp);
    _f2.computePropertiesAtQp(_qp);

    residual(0) = _first_df1[_qp] - _first_df2[_qp];
    residual(1) = _c[_qp] - guess(1) * _h[_qp] + guess(0) * (_h[_qp] - 1);

    jacobian(0, 0) = _second_df1[_qp];
    jacobian(0, 1) = -_second_df2[_qp];
    jacobian(1, 0) = _h[_qp] - 1;
    jacobian(1, 1) = -_h[_qp];
  };

  solver.nonlinear(solution, compute);

  if (solver.getState() == NestedSolve::State::NOT_CONVERGED)
  {
    std::cout << "Newton iteration did not converge." << std::endl;
  }

  // if (solver.getState() == NestedSolve::State::CONVERGED_REL)
  // {
  //   std::cout << "Newton iteration converged." << std::endl;
  // }

  _c1[_qp] = solution[0];
  _c2[_qp] = solution[1];
  _f1.computePropertiesAtQp(_qp);
  _f2.computePropertiesAtQp(_qp);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////// Powell's dogleg begin
  // auto computeResidual = [&](const NestedSolve::Value<> & guess, NestedSolve::Value<> & residual)
  // {
  //   _c1[_qp] = guess(0);
  //   _c2[_qp] = guess(1);
  //   _f1.computePropertiesAtQp(_qp);
  //   _f2.computePropertiesAtQp(_qp);
  //
  //   residual(0) = _first_df1[_qp] - _first_df2[_qp];
  //   residual(1) = _c[_qp] - guess(1) * _h[_qp] + guess(0) * (_h[_qp] - 1);
  // };
  //
  // auto computeJacobian = [&](const NestedSolve::Value<> & guess,
  //                            NestedSolve::Jacobian<> & jacobian) {
  //   _c1[_qp] = guess(0);
  //   _c2[_qp] = guess(1);
  //   _f1.computePropertiesAtQp(_qp);
  //   _f2.computePropertiesAtQp(_qp);
  //
  //   jacobian(0, 0) = _second_df1[_qp];
  //   jacobian(0, 1) = -_second_df2[_qp];
  //   jacobian(1, 0) = _h[_qp] - 1;
  //   jacobian(1, 1) = -_h[_qp];
  // };
  //
  // solver.nonlinear(solution, computeResidual, computeJacobian);
  //
  // if (solver.getState() == NestedSolve::State::NOT_CONVERGED)
  // {
  //   std::cout << "Newton iteration did not converge." << std::endl;
  // }
  //
  // _c1[_qp] = solution[0];
  // _c2[_qp] = solution[1];
  // _f1.computePropertiesAtQp(_qp);
  // _f2.computePropertiesAtQp(_qp);
  //////////////////////////////////////////////////////////////////////////////////////////////////////////// Powell's dogleg end

  // compute dc1dc, dc2dc, dc1deta, and dc2deta
  _dc1dc[_qp] =
      _second_df2[_qp] / (_second_df2[_qp] + _h[_qp] * (_second_df1[_qp] - _second_df2[_qp]));
  _dc2dc[_qp] =
      _second_df1[_qp] / (_second_df2[_qp] + _h[_qp] * (_second_df1[_qp] - _second_df2[_qp]));

  _dc1deta[_qp] = _second_df2[_qp] * (_c1[_qp] - _c2[_qp]) *
                  (30.0 * n * n * (n * n - 2.0 * n + 1.0)) /
                  (_second_df2[_qp] + _h[_qp] * (_second_df1[_qp] - _second_df2[_qp]));
  _dc2deta[_qp] = _second_df1[_qp] * (_c1[_qp] - _c2[_qp]) *
                  (30.0 * n * n * (n * n - 2.0 * n + 1.0)) /
                  (_second_df2[_qp] + _h[_qp] * (_second_df1[_qp] - _second_df2[_qp]));
}
