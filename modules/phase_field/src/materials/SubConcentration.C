//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SubConcentration.h"
#include "libmesh/utility.h"
// #include "MathUtils.h"
#include <cmath>

registerMooseObject("PhaseFieldApp", SubConcentration);

InputParameters
SubConcentration::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Computes the KKS sub-concentrations by using Newton iteration to solve the equal chemical "
      "potential and concentration conservation equations");
  params.addRequiredCoupledVar("global_c", "The interpolated concentration");
  params.addRequiredCoupledVar("eta_name", "The order parameter");
  params.addRequiredParam<MaterialPropertyName>("h_name", "Name of the switching function");
  params.addRequiredParam<MaterialPropertyName>("c1_name", "Name of the first phase concentration");
  params.addRequiredParam<MaterialPropertyName>("c2_name",
                                                "Name of the second phase concentration");
  params.addParam<Real>("absolute_tol_value", 1e-9, "Absolute tolerance of the Newton iteration");
  params.addParam<Real>("relative_tol_value", 1e-9, "Relative tolerance of the Newton iteration");
  params.addParam<Real>("max_iteration", 100, "The maximum number of Newton iterations");
  params.addRequiredParam<MaterialPropertyName>("dc1dc_name",
                                                "The approximated first derivative of c1 w.r.t. c");
  params.addRequiredParam<MaterialPropertyName>(
      "dc1deta_name", "The approximated first derivative of c1 w.r.t. eta");
  params.addRequiredParam<MaterialPropertyName>("dc2dc_name",
                                                "The approximated first derivative of c2 w.r.t. c");
  params.addRequiredParam<MaterialPropertyName>(
      "dc2deta_name", "The approximated first derivative of c2 w.r.t. eta");
  return params;
}

SubConcentration::SubConcentration(const InputParameters & parameters)
  : Material(parameters),
    _c(coupledValue("global_c")),
    _eta(coupledValue("eta_name")),
    _h(getMaterialProperty<Real>("h_name")),
    _c1(declareProperty<Real>("c1_name")),
    _c2(declareProperty<Real>("c2_name")),
    // _c1_old(getMaterialPropertyOld<Real>("c1_name")), // old
    // _c2_old(getMaterialPropertyOld<Real>("c2_name")), // old
    _abs_tol(getParam<Real>("absolute_tol_value")),
    _rel_tol(getParam<Real>("relative_tol_value")),
    _maxiter(getParam<Real>("max_iteration")),
    _dc1dc(declareProperty<Real>("dc1dc_name")),
    _dc1deta(declareProperty<Real>("dc1deta_name")),
    _dc2dc(declareProperty<Real>("dc2dc_name")),
    _dc2deta(declareProperty<Real>("dc2deta_name"))

{
}

// void
// SubConcentration::initQpStatefulProperties()
// {
//   _c1[_qp] = 0.6;
//   _c2[_qp] = 0.4;
//
//   // _c_mat[_qp] = _c[_qp];
//   // _eta_mat[_qp] = _eta[_qp];
//   _c_mat[_qp] = 0.5;
//   _eta_mat[_qp] = 0.5;
// }

void
SubConcentration::computeQpProperties()
{

  // old ci inisde Newton iteration
  std::vector<Real> old_ci_Newton(2);
  old_ci_Newton[0] = 0.6;
  old_ci_Newton[1] = 0.4;

  std::vector<Real> init_err(2);
  std::vector<Real> abs_err(2);

  Real init_err_norm;
  Real abs_err_norm;
  Real rel_err_norm;

  init_err[0] = 400 * log(1 - old_ci_Newton[1]) - 400 * log(1 - old_ci_Newton[0]) +
                400 * log(old_ci_Newton[0]) - 400 * log(old_ci_Newton[1]) - 67;
  init_err[1] = _c[_qp] - old_ci_Newton[1] * _h[_qp] + old_ci_Newton[0] * (_h[_qp] - 1);
  init_err_norm = std::sqrt(Utility::pow<2>(init_err[0]) + Utility::pow<2>(init_err[1]));

  // Real count;
  // count = 0;

  for (unsigned int nloop = 0; nloop < _maxiter; ++nloop)
  {
    _c1[_qp] =
        old_ci_Newton[0] -
        ((old_ci_Newton[0] * (old_ci_Newton[0] - 1) *
          (_c[_qp] - old_ci_Newton[1] * _h[_qp] + old_ci_Newton[0] * (_h[_qp] - 1))) /
             (old_ci_Newton[0] - old_ci_Newton[0] * _h[_qp] + old_ci_Newton[1] * _h[_qp] +
              Utility::pow<2>(old_ci_Newton[0]) * _h[_qp] -
              Utility::pow<2>(old_ci_Newton[1]) * _h[_qp] - Utility::pow<2>(old_ci_Newton[0])) -
         (old_ci_Newton[0] * old_ci_Newton[1] * _h[_qp] * (old_ci_Newton[0] - 1) *
          (old_ci_Newton[1] - 1) *
          (400 * log(1 - old_ci_Newton[0]) - 400 * log(1 - old_ci_Newton[1]) -
           400 * log(old_ci_Newton[0]) + 400 * log(old_ci_Newton[1]) + 67)) /
             (400 *
              (old_ci_Newton[0] - old_ci_Newton[0] * _h[_qp] + old_ci_Newton[1] * _h[_qp] +
               Utility::pow<2>(old_ci_Newton[0]) * _h[_qp] -
               Utility::pow<2>(old_ci_Newton[1]) * _h[_qp] - Utility::pow<2>(old_ci_Newton[0]))));

    _c2[_qp] =
        old_ci_Newton[1] -
        ((old_ci_Newton[1] * (old_ci_Newton[1] - 1) *
          (_c[_qp] - old_ci_Newton[1] * _h[_qp] + old_ci_Newton[0] * (_h[_qp] - 1))) /
             (old_ci_Newton[0] - old_ci_Newton[0] * _h[_qp] + old_ci_Newton[1] * _h[_qp] +
              Utility::pow<2>(old_ci_Newton[0]) * _h[_qp] -
              Utility::pow<2>(old_ci_Newton[1]) * _h[_qp] - Utility::pow<2>(old_ci_Newton[0])) -
         (old_ci_Newton[0] * old_ci_Newton[1] * (old_ci_Newton[0] - 1) * (old_ci_Newton[1] - 1) *
          (_h[_qp] - 1) *
          (400 * log(1 - old_ci_Newton[0]) - 400 * log(1 - old_ci_Newton[1]) -
           400 * log(old_ci_Newton[0]) + 400 * log(old_ci_Newton[1]) + 67)) /
             (400 *
              (old_ci_Newton[0] - old_ci_Newton[0] * _h[_qp] + old_ci_Newton[1] * _h[_qp] +
               Utility::pow<2>(old_ci_Newton[0]) * _h[_qp] -
               Utility::pow<2>(old_ci_Newton[1]) * _h[_qp] - Utility::pow<2>(old_ci_Newton[0]))));

    // std::cout << "eta is " << _eta[_qp] << std::endl;
    // std::cout << "c1 is " << _c1[_qp] << std::endl;
    // std::cout << "c2 is " << _c2[_qp] << std::endl;
    // std::cout << "dc1dc is " << _dc1dc[_qp] << std::endl;
    // std::cout << "dc2dc is " << _dc2dc[_qp] << std::endl;
    // std::cout << "dc1deta is " << _dc1deta[_qp] << std::endl;
    // std::cout << "dc2deta is " << _dc2deta[_qp] << std::endl;

    // compute the absolute Newton error
    abs_err[0] = 400 * log(1 - _c2[_qp]) - 400 * log(1 - _c1[_qp]) + 400 * log(_c1[_qp]) -
                 400 * log(_c2[_qp]) - 67;
    abs_err[1] = _c[_qp] - _c2[_qp] * _h[_qp] + _c1[_qp] * (_h[_qp] - 1);
    abs_err_norm = std::sqrt(Utility::pow<2>(abs_err[0]) + Utility::pow<2>(abs_err[1]));

    // compute the relative Newton error
    rel_err_norm = std::abs(abs_err_norm / init_err_norm);

    // std::cout << "current count is " << count << std::endl;
    // count += 1;

    // std::cout << "The initial error is " << init_err_norm << std::endl;
    // std::cout << "The absolute error is " << abs_err_norm << std::endl;
    // std::cout << "The relative error is " << rel_err_norm << std::endl;

    // update ci
    old_ci_Newton[0] = _c1[_qp];
    old_ci_Newton[1] = _c2[_qp];

    // Newton iteration convergence criterion
    if (abs_err_norm < _abs_tol)
      break;
    else if (rel_err_norm < _rel_tol)
      break;
  }

  // compute dc1/dc, dc1/deta, dc2/dc, and dc2/deta
  _dc1dc[_qp] =
      -(_c1[_qp] * (_c1[_qp] - 1)) /
      (_c1[_qp] - _c1[_qp] * _h[_qp] + _c2[_qp] * _h[_qp] + Utility::pow<2>(_c1[_qp]) * _h[_qp] -
       Utility::pow<2>(_c2[_qp]) * _h[_qp] - Utility::pow<2>(_c1[_qp]));

  _dc1deta[_qp] =
      -(_c1[_qp] * (_c1[_qp] - _c2[_qp]) * 30.0 * _eta[_qp] * _eta[_qp] *
        (_eta[_qp] * _eta[_qp] - 2.0 * _eta[_qp] + 1.0) * (_c1[_qp] - 1)) /
      (_c1[_qp] - _c1[_qp] * _h[_qp] + _c2[_qp] * _h[_qp] + Utility::pow<2>(_c1[_qp]) * _h[_qp] -
       Utility::pow<2>(_c2[_qp]) * _h[_qp] - Utility::pow<2>(_c1[_qp]));

  _dc2dc[_qp] =
      -(_c2[_qp] * (_c2[_qp] - 1)) /
      (_c1[_qp] - _c1[_qp] * _h[_qp] + _c2[_qp] * _h[_qp] + Utility::pow<2>(_c1[_qp]) * _h[_qp] -
       Utility::pow<2>(_c2[_qp]) * _h[_qp] - Utility::pow<2>(_c1[_qp]));

  _dc2deta[_qp] =
      -(_c2[_qp] * (_c1[_qp] - _c2[_qp]) * 30.0 * _eta[_qp] * _eta[_qp] *
        (_eta[_qp] * _eta[_qp] - 2.0 * _eta[_qp] + 1.0) * (_c2[_qp] - 1)) /
      (_c1[_qp] - _c1[_qp] * _h[_qp] + _c2[_qp] * _h[_qp] + Utility::pow<2>(_c1[_qp]) * _h[_qp] -
       Utility::pow<2>(_c2[_qp]) * _h[_qp] - Utility::pow<2>(_c1[_qp]));
}
