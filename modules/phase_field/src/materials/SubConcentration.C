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
  params.addRequiredParam<MaterialPropertyName>("h_name", "Name of the switching function");
  params.addRequiredParam<MaterialPropertyName>("c1_name", "Name of the first phase concentration");
  params.addRequiredParam<MaterialPropertyName>("c2_name",
                                                "Name of the second phase concentration");
  params.addParam<Real>("absolute_tol_value", 1e-9, "Absolute tolerance of the Newton iteration");
  params.addParam<Real>("relative_tol_value", 1e-9, "Relative tolerance of the Newton iteration");
  params.addParam<Real>("max_iteration", 100, "The maximum number of Newton iterations");
  return params;
}

SubConcentration::SubConcentration(const InputParameters & parameters)
  : Material(parameters),
    _c(coupledValue("global_c")),
    _h(getMaterialProperty<Real>("h_name")),
    _c1(declareProperty<Real>("c1_name")),
    _c2(declareProperty<Real>("c2_name")),
    // _c1_old(getMaterialPropertyOld<Real>("c1_name")),
    // _c2_old(getMaterialPropertyOld<Real>("c2_name")),
    _abs_tol(getParam<Real>("absolute_tol_value")),
    _rel_tol(getParam<Real>("relative_tol_value")),
    _maxiter(getParam<Real>("max_iteration"))

{
}

// void
// SubConcentration::initQpStatefulProperties()
// {
//   _c1[_qp] = 0.4;
//   _c2[_qp] = 0.6;
// }

void
SubConcentration::computeQpProperties()
{
  std::vector<Real> old_ci_Newton(2);

  // old_ci_Newton[0] = _c1_old[_qp];
  // old_ci_Newton[1] = _c2_old[_qp];

  old_ci_Newton[0] = 0.4;
  old_ci_Newton[1] = 0.6;

  std::vector<Real> init_err(2);
  std::vector<Real> abs_err(2);

  Real init_err_norm;
  Real abs_err_norm;
  Real rel_err_norm;

  // initial error norm
  init_err[0] = 200 * old_ci_Newton[0] - 200 * old_ci_Newton[1] + 80;
  init_err[1] = _c[_qp] - old_ci_Newton[1] * _h[_qp] + old_ci_Newton[0] * (_h[_qp] - 1);
  init_err_norm = std::sqrt(Utility::pow<2>(init_err[0]) + Utility::pow<2>(init_err[1]));

  for (unsigned int nloop = 0; nloop < _maxiter; ++nloop)
  {
    _c1[_qp] = old_ci_Newton[0] - (old_ci_Newton[0] - _c[_qp] + (2 * _h[_qp]) / 5);

    _c2[_qp] = old_ci_Newton[1] - (old_ci_Newton[1] - _c[_qp] + (2 * _h[_qp]) / 5 - 0.4);

    // absolute error norm
    abs_err[0] = 200 * _c1[_qp] - 200 * _c2[_qp] + 80;
    abs_err[1] = _c[_qp] - _c2[_qp] * _h[_qp] + _c1[_qp] * (_h[_qp] - 1);
    abs_err_norm = std::sqrt(Utility::pow<2>(abs_err[0]) + Utility::pow<2>(abs_err[1]));

    // relative error norm
    rel_err_norm = std::abs(abs_err_norm / init_err_norm);

    // update ci
    old_ci_Newton[0] = _c1[_qp];
    old_ci_Newton[1] = _c2[_qp];

    // std::cout << "Newton iteration loop " << nloop << ", the initial error norm is "
    //           << init_err_norm << ", the absolute error norm is " << abs_err_norm
    //           << ", and the relative error norm is " << rel_err_norm << std::endl;
    //
    // std::cout << "c1 is " << _c1[_qp] << ", and c2 is " << _c2[_qp] << std::endl;

    if (abs_err_norm < _abs_tol)
      break;
    else if (rel_err_norm < _rel_tol)
      break;
  }
}
