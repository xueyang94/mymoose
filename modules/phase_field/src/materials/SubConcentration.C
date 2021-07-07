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
  // InputParameters params = emptyInputParameters();
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Computes the KKS sub-concentrations by using Newton iteration to solve the equal chemical "
      "potential and concentration conservation equations");
  params.addRequiredCoupledVar("global_c", "The interpolated concentration");
  params.addRequiredParam<MaterialPropertyName>("h_name", "Name of the switching function");
  params.addRequiredParam<std::vector<MaterialPropertyName>>("ci_names",
                                                             "Names of the sub-concentrations");
  params.addParam<Real>("absolute_tol_value", 1e-9, "Absolute tolerance of the Newton iteration");
  params.addParam<Real>("relative_tol_value", 1e-9, "Relative tolerance of the Newton iteration");
  params.addParam<Real>("max_iteration", 100, "The maximum number of Newton iterations");
  return params;
}

SubConcentration::SubConcentration(const InputParameters & parameters)
  : Material(parameters),
    // _c(coupledValues("global_c")),
    _c(coupledValue("global_c")),
    _h(getMaterialProperty<Real>("h_name")),
    _ci_name(getParam<std::vector<MaterialPropertyName>>("ci_names")),
    _ci_prop(2),
    _abs_tol(getParam<Real>("absolute_tol_value")),
    _rel_tol(getParam<Real>("relative_tol_value")),
    _maxiter(getParam<Real>("max_iteration"))

{
  // declare ci material properties
  for (unsigned int i = 0; i < 2; i++)
    _ci_prop[i] = &declareProperty<Real>(_ci_name[i]);
}

// unsigned int n_qp{0};

void
SubConcentration::computeQpProperties()
{
  // n_qp += 1;
  // std::cout << "qp point " << n_qp << std::endl;

  std::vector<Real> old_ci{0.4, 0.6};

  std::vector<Real> new_ci(2);

  std::vector<Real> init_err(2);
  std::vector<Real> abs_err(2);
  // std::vector<Real> rel_err(2);
  float init_err_norm;
  float abs_err_norm;
  float rel_err_norm;

  init_err[0] = 200 * old_ci[0] - 200 * old_ci[1] + 80;
  init_err[1] = _c[_qp] - old_ci[1] * _h[_qp] + old_ci[0] * (_h[_qp] - 1);
  init_err_norm = std::sqrt(Utility::pow<2>(init_err[0]) + Utility::pow<2>(init_err[1]));

  for (unsigned int nloop = 0; nloop < _maxiter; ++nloop)
  {
    new_ci[0] = old_ci[0] - (old_ci[0] - _c[_qp] + (2 * _h[_qp]) / 5);

    new_ci[1] = old_ci[1] - (old_ci[1] - _c[_qp] + (2 * _h[_qp]) / 5 - 0.4);

    abs_err[0] = 200 * new_ci[0] - 200 * new_ci[1] + 80;
    abs_err[1] = _c[_qp] - new_ci[1] * _h[_qp] + new_ci[0] * (_h[_qp] - 1);
    abs_err_norm = std::sqrt(Utility::pow<2>(abs_err[0]) + Utility::pow<2>(abs_err[1]));

    rel_err_norm = std::abs(abs_err_norm / init_err_norm);

    // std::cout << "Newton iteration loop " << nloop << ", the absolute error norm is "
    //           << abs_err_norm << ", and the relative error norm is " << rel_err_norm <<
    //           std::endl;

    old_ci = new_ci; // update ci

    if (abs_err_norm < _abs_tol)
      break;
    else if (rel_err_norm < _rel_tol)
      break;
  }

  // std::cout << "c1 is " << new_ci[0] << ", and c2 is " << new_ci[1] << ", c is " << _c[_qp]
  //           << ", and h is " << _h[_qp] << std::endl;

  for (unsigned int i = 0; i < 2; ++i)
    (*_ci_prop[i])[_qp] = new_ci[i];
}
