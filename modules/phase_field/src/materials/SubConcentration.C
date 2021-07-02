//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SubConcentration.h"

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
  params.addParam<Real>("tol_value", 1e-9, "Tolerance of the Newton iteration");
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
    _tol(getParam<Real>("tol_value")),
    _maxiter(getParam<Real>("max_iteration"))

{
  // declare ci material properties
  for (unsigned int i = 0; i < 2; i++)
    _ci_prop[i] = &declareProperty<Real>(_ci_name[i]);
}

void
SubConcentration::computeQpProperties()
{
  // std::vector<Real> old_ci{0.668, 0.000045399}; // IC of ci
  std::vector<Real> old_ci{0.4, 0.8}; // IC of ci

  std::vector<Real> new_ci(2);

  float err;

  for (unsigned int nloop = 0; nloop < _maxiter; ++nloop)
  {

    // new_ci[0] = old_ci[0] - ((old_ci[0] * (old_ci[0] - 1) *
    //                           (_c[_qp] - old_ci[1] * _h[_qp] + old_ci[0] * (_h[_qp] - 1))) /
    //                              (old_ci[0] - old_ci[0] * _h[_qp] + old_ci[1] * _h[_qp] +
    //                               old_ci[0] * old_ci[0] * _h[_qp] -
    //                               old_ci[1] * old_ci[1] * _h[_qp] - old_ci[0] * old_ci[0]) -
    //                          (old_ci[0] * old_ci[1] * _h[_qp] * (old_ci[0] - 1) * (old_ci[1] - 1)
    //                          *
    //                           (400 * log(1 - old_ci[0]) - 400 * log(1 - old_ci[1]) -
    //                            400 * log(old_ci[0]) + 400 * log(old_ci[1]) + 427999 / 100)) /
    //                              (400 * (old_ci[0] - old_ci[0] * _h[_qp] + old_ci[1] * _h[_qp] +
    //                                      old_ci[0] * old_ci[0] * _h[_qp] -
    //                                      old_ci[1] * old_ci[1] * _h[_qp] - old_ci[0] *
    //                                      old_ci[0])));
    // new_ci[1] =
    //     old_ci[1] - ((old_ci[1] * (old_ci[1] - 1) *
    //                   (_c[_qp] - old_ci[1] * _h[_qp] + old_ci[0] * (_h[_qp] - 1))) /
    //                      (old_ci[0] - old_ci[0] * _h[_qp] + old_ci[1] * _h[_qp] +
    //                       old_ci[0] * old_ci[0] * _h[_qp] - old_ci[1] * old_ci[1] * _h[_qp] -
    //                       old_ci[0] * old_ci[0]) -
    //                  (old_ci[0] * old_ci[1] * (old_ci[0] - 1) * (old_ci[1] - 1) * (_h[_qp] - 1) *
    //                   (400 * log(1 - old_ci[0]) - 400 * log(1 - old_ci[1]) - 400 * log(old_ci[0])
    //                   +
    //                    400 * log(old_ci[1]) + 427999 / 100)) /
    //                      (400 * (old_ci[0] - old_ci[0] * _h[_qp] + old_ci[1] * _h[_qp] +
    //                              old_ci[0] * old_ci[0] * _h[_qp] - old_ci[1] * old_ci[1] *
    //                              _h[_qp] - old_ci[0] * old_ci[0])));

    new_ci[0] = old_ci[0] - (old_ci[1] * _h[_qp] - _c[_qp] - old_ci[0] * (_h[_qp] - 1) +
                             (_h[_qp] * (200 * old_ci[0] - 200 * old_ci[1] + 80)) / 200);

    new_ci[1] = old_ci[1] - (old_ci[1] * _h[_qp] - _c[_qp] +
                             (_h[_qp] / 200 - 1 / 200) * (200 * old_ci[0] - 200 * old_ci[1] + 80) -
                             old_ci[0] * (_h[_qp] - 1));

    err = sqrt((new_ci[0] - old_ci[0]) * (new_ci[0] - old_ci[0]) +
               (new_ci[1] - old_ci[1]) * (new_ci[1] - old_ci[1]));

    old_ci = new_ci; // update ci

    // std::cout << "Newton iteration loop " << nloop << ", error norm is " << err << std::endl;

    if (err < _tol)
      break;
  }

  // std::cout << "c1 is " << new_ci[0] << ", and c2 is " << new_ci[1] << std::endl;
  //
  // std::cout << "Next qp point" << std::endl;

  for (unsigned int i = 0; i < 2; ++i)
    (*_ci_prop[i])[_qp] = new_ci[i];

  // (*_ci_prop[0])[_qp] = 0.668;
  // (*_ci_prop[1])[_qp] = 0.000045399;
}
