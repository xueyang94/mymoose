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

// define a log Taylor expansion function
// Real
// tlog(Real x)
// {
//   Real approx;
//   approx =
//       (x - 1) - Utility::pow<2>(x - 1) / 2 + Utility::pow<3>(x - 1) / 3 -
//       Utility::pow<4>(x - 1) / 4 + Utility::pow<5>(x - 1) / 5 - Utility::pow<6>(x - 1) / 6 +
//       Utility::pow<7>(x - 1) / 7 - Utility::pow<8>(x - 1) / 8 + Utility::pow<9>(x - 1) / 9 -
//       Utility::pow<10>(x - 1) / 10 + Utility::pow<11>(x - 1) / 11 - Utility::pow<12>(x - 1) / 12
//       + Utility::pow<13>(x - 1) / 13 - Utility::pow<14>(x - 1) / 14 + Utility::pow<15>(x - 1) /
//       15 - Utility::pow<16>(x - 1) / 16 + Utility::pow<17>(x - 1) / 17 - Utility::pow<18>(x - 1)
//       / 18 + Utility::pow<19>(x - 1) / 19 - Utility::pow<20>(x - 1) / 20 + Utility::pow<21>(x -
//       1) / 21 - Utility::pow<22>(x - 1) / 22 + Utility::pow<23>(x - 1) / 23 - Utility::pow<24>(x
//       - 1) / 24 + Utility::pow<25>(x - 1) / 25 - Utility::pow<26>(x - 1) / 26 +
//       Utility::pow<27>(x - 1) / 27 - Utility::pow<28>(x - 1) / 28 + Utility::pow<29>(x - 1) / 29
//       - Utility::pow<30>(x - 1) / 30;
//   return approx;
// }

Real
tlog(Real x)
{
  Real approx;
  approx = (x - 1) - Utility::pow<2>(x - 1) / 2;
  return approx;
}

unsigned int n_qp{0};

void
SubConcentration::computeQpProperties()
{

  // std::vector<Real> old_ci{0.668, 0.000045399}; // IC of ci
  std::vector<Real> old_ci{0.6, 0.1}; // IC of ci
  // std::vector<Real> old_ci{0.6, 0.4};

  std::vector<Real> new_ci(2);

  std::vector<Real> init_err(2);
  std::vector<Real> abs_err(2);
  // std::vector<Real> rel_err(2);
  float init_err_norm;
  float abs_err_norm;
  float rel_err_norm;

  // init_err[0] = 400 * log(1 - old_ci[1]) - 400 * log(1 - old_ci[0]) + 400 * log(old_ci[0]) -
  //               400 * log(old_ci[1]) - 4279.99;
  // init_err[1] = _c[_qp] - old_ci[1] * _h[_qp] + old_ci[0] * (_h[_qp] - 1);
  init_err[0] = 800 * old_ci[0] - 800 * old_ci[1] +
                400 * old_ci[0] * (Utility::pow<2>(old_ci[0] - 1) - old_ci[0] + 2) -
                400 * old_ci[1] * (Utility::pow<2>(old_ci[1] - 1) - old_ci[1] + 2) -
                200 * Utility::pow<2>(old_ci[0] - 1) + (400 * Utility::pow<3>(old_ci[0] - 1)) / 3 +
                200 * Utility::pow<2>(old_ci[1] - 1) - (400 * Utility::pow<3>(old_ci[1] - 1)) / 3 +
                400 * (old_ci[0] - 1) * (Utility::pow<2>(old_ci[0]) + old_ci[0] + 1) -
                400 * (old_ci[1] - 1) * (Utility::pow<2>(old_ci[1]) + old_ci[1] + 1) +
                200 * Utility::pow<2>(old_ci[0]) + (400 * Utility::pow<3>(old_ci[0])) / 3 -
                200 * Utility::pow<2>(old_ci[1]) - (400 * Utility::pow<3>(old_ci[1])) / 3 - 4279.99;
  init_err[1] = _c[_qp] - old_ci[1] * _h[_qp] + old_ci[0] * (_h[_qp] - 1);
  init_err_norm = std::sqrt(Utility::pow<2>(init_err[0]) + Utility::pow<2>(init_err[1]));

  // init_err[0] = 400 * log(1 - old_ci[1]) - 400 * log(1 - old_ci[0]) + 400 * log(old_ci[0]) -
  //               400 * log(old_ci[1]) - 67;
  // init_err[1] = _c[_qp] - old_ci[1] * _h[_qp] + old_ci[0] * (_h[_qp] - 1);
  // init_err_norm = std::sqrt(Utility::pow<2>(init_err[0]) + Utility::pow<2>(init_err[1]));

  for (unsigned int nloop = 0; nloop < _maxiter; ++nloop)
  {
    // new_ci[0] = old_ci[0] -
    //             ((old_ci[0] * (old_ci[0] - 1) *
    //               (_c[_qp] - old_ci[1] * _h[_qp] + old_ci[0] * (_h[_qp] - 1))) /
    //                  (old_ci[0] - old_ci[0] * _h[_qp] + old_ci[1] * _h[_qp] +
    //                   Utility::pow<2>(old_ci[0]) * _h[_qp] - Utility::pow<2>(old_ci[1]) * _h[_qp]
    //                   - Utility::pow<2>(old_ci[0])) -
    //              (old_ci[0] * old_ci[1] * _h[_qp] * (old_ci[0] - 1) * (old_ci[1] - 1) *
    //               (400 * log(1 - old_ci[0]) - 400 * log(1 - old_ci[1]) - 400 * log(old_ci[0]) +
    //                400 * log(old_ci[1]) + 4279.99)) /
    //                  (400 * (old_ci[0] - old_ci[0] * _h[_qp] + old_ci[1] * _h[_qp] +
    //                          Utility::pow<2>(old_ci[0]) * _h[_qp] -
    //                          Utility::pow<2>(old_ci[1]) * _h[_qp] -
    //                          Utility::pow<2>(old_ci[0]))));

    new_ci[0] = old_ci[0] - (-(840000 * _c[_qp] - 840000 * old_ci[0] + 1283997 * _h[_qp] -
                               960000 * _c[_qp] * old_ci[1] + 960000 * old_ci[0] * old_ci[1] +
                               960000 * _c[_qp] * Utility::pow<2>(old_ci[1]) -
                               960000 * old_ci[0] * Utility::pow<2>(old_ci[1]) +
                               480000 * Utility::pow<2>(old_ci[0]) * _h[_qp] -
                               320000 * Utility::pow<3>(old_ci[0]) * _h[_qp] +
                               480000 * Utility::pow<2>(old_ci[1]) * _h[_qp] -
                               640000 * Utility::pow<3>(old_ci[1]) * _h[_qp] -
                               960000 * old_ci[0] * old_ci[1] * _h[_qp] +
                               960000 * old_ci[0] * Utility::pow<2>(old_ci[1]) * _h[_qp]) /
                             (120000 * (8 * old_ci[1] * _h[_qp] - 8 * old_ci[0] * _h[_qp] -
                                        8 * old_ci[1] + 8 * Utility::pow<2>(old_ci[0]) * _h[_qp] -
                                        8 * Utility::pow<2>(old_ci[1]) * _h[_qp] +
                                        8 * Utility::pow<2>(old_ci[1]) + 7)));

    // new_ci[1] = old_ci[1] -
    //             ((old_ci[1] * (old_ci[1] - 1) *
    //               (_c[_qp] - old_ci[1] * _h[_qp] + old_ci[0] * (_h[_qp] - 1))) /
    //                  (old_ci[0] - old_ci[0] * _h[_qp] + old_ci[1] * _h[_qp] +
    //                   Utility::pow<2>(old_ci[0]) * _h[_qp] - Utility::pow<2>(old_ci[1]) * _h[_qp]
    //                   - Utility::pow<2>(old_ci[0])) -
    //              (old_ci[0] * old_ci[1] * (old_ci[0] - 1) * (old_ci[1] - 1) * (_h[_qp] - 1) *
    //               (400 * log(1 - old_ci[0]) - 400 * log(1 - old_ci[1]) - 400 * log(old_ci[0]) +
    //                400 * log(old_ci[1]) + 4279.99)) /
    //                  (400 * (old_ci[0] - old_ci[0] * _h[_qp] + old_ci[1] * _h[_qp] +
    //                          Utility::pow<2>(old_ci[0]) * _h[_qp] -
    //                          Utility::pow<2>(old_ci[1]) * _h[_qp] -
    //                          Utility::pow<2>(old_ci[0]))));

    new_ci[1] =
        old_ci[1] -
        (-(840000 * _c[_qp] - 840000 * old_ci[1] + 1283997 * _h[_qp] -
           960000 * _c[_qp] * old_ci[0] + 960000 * _c[_qp] * Utility::pow<2>(old_ci[0]) -
           480000 * Utility::pow<2>(old_ci[0]) * _h[_qp] +
           640000 * Utility::pow<3>(old_ci[0]) * _h[_qp] -
           480000 * Utility::pow<2>(old_ci[1]) * _h[_qp] +
           320000 * Utility::pow<3>(old_ci[1]) * _h[_qp] + 480000 * Utility::pow<2>(old_ci[0]) -
           640000 * Utility::pow<3>(old_ci[0]) + 480000 * Utility::pow<2>(old_ci[1]) -
           320000 * Utility::pow<3>(old_ci[1]) + 960000 * old_ci[0] * old_ci[1] * _h[_qp] -
           960000 * Utility::pow<2>(old_ci[0]) * old_ci[1] * _h[_qp] - 1283997) /
         (120000 *
          (8 * old_ci[1] * _h[_qp] - 8 * old_ci[0] * _h[_qp] - 8 * old_ci[1] +
           8 * Utility::pow<2>(old_ci[0]) * _h[_qp] - 8 * Utility::pow<2>(old_ci[1]) * _h[_qp] +
           8 * Utility::pow<2>(old_ci[1]) + 7)));

    // new_ci[0] = old_ci[0] -
    //             ((old_ci[0] * (old_ci[0] - 1) *
    //               (_c[_qp] - old_ci[1] * _h[_qp] + old_ci[0] * (_h[_qp] - 1))) /
    //                  (old_ci[0] - old_ci[0] * _h[_qp] + old_ci[1] * _h[_qp] +
    //                   Utility::pow<2>(old_ci[0]) * _h[_qp] - Utility::pow<2>(old_ci[1]) * _h[_qp]
    //                   - Utility::pow<2>(old_ci[0])) -
    //              (old_ci[0] * old_ci[1] * _h[_qp] * (old_ci[0] - 1) * (old_ci[1] - 1) *
    //               (400 * log(1 - old_ci[0]) - 400 * log(1 - old_ci[1]) - 400 * log(old_ci[0]) +
    //                400 * log(old_ci[1]) + 67)) /
    //                  (400 * (old_ci[0] - old_ci[0] * _h[_qp] + old_ci[1] * _h[_qp] +
    //                          Utility::pow<2>(old_ci[0]) * _h[_qp] -
    //                          Utility::pow<2>(old_ci[1]) * _h[_qp] -
    //                          Utility::pow<2>(old_ci[0]))));
    //
    // new_ci[1] = old_ci[1] -
    //             ((old_ci[1] * (old_ci[1] - 1) *
    //               (_c[_qp] - old_ci[1] * _h[_qp] + old_ci[0] * (_h[_qp] - 1))) /
    //                  (old_ci[0] - old_ci[0] * _h[_qp] + old_ci[1] * _h[_qp] +
    //                   Utility::pow<2>(old_ci[0]) * _h[_qp] - Utility::pow<2>(old_ci[1]) * _h[_qp]
    //                   - Utility::pow<2>(old_ci[0])) -
    //              (old_ci[0] * old_ci[1] * (old_ci[0] - 1) * (old_ci[1] - 1) * (_h[_qp] - 1) *
    //               (400 * log(1 - old_ci[0]) - 400 * log(1 - old_ci[1]) - 400 * log(old_ci[0]) +
    //                400 * log(old_ci[1]) + 67)) /
    //                  (400 * (old_ci[0] - old_ci[0] * _h[_qp] + old_ci[1] * _h[_qp] +
    //                          Utility::pow<2>(old_ci[0]) * _h[_qp] -
    //                          Utility::pow<2>(old_ci[1]) * _h[_qp] -
    //                          Utility::pow<2>(old_ci[0]))));

    // abs_err[0] = 400 * log(1 - new_ci[1]) - 400 * log(1 - new_ci[0]) + 400 * log(new_ci[0]) -
    //              400 * log(new_ci[1]) - 4279.99;
    // abs_err[1] = _c[_qp] - new_ci[1] * _h[_qp] + new_ci[0] * (_h[_qp] - 1);
    // two terms of Taylor
    // abs_err[0] = 800 * old_ci[0] - 800 * old_ci[1] +
    //               400 * old_ci[0] * (Utility::pow<2>(old_ci[0] - 1) - old_ci[0] + 2) -
    //               400 * old_ci[1] * (Utility::pow<2>(old_ci[1] - 1) - old_ci[1] + 2) -
    //               200 * Utility::pow<2>(old_ci[0] - 1) + (400 * Utility::pow<3>(old_ci[0] - 1)) /
    //               3 + 200 * Utility::pow<2>(old_ci[1] - 1) - (400 * Utility::pow<3>(old_ci[1] -
    //               1)) / 3 + 400 * (old_ci[0] - 1) * (Utility::pow<2>(old_ci[0]) + old_ci[0] + 1)
    //               - 400 * (old_ci[1] - 1) * (Utility::pow<2>(old_ci[1]) + old_ci[1] + 1) + 200 *
    //               Utility::pow<2>(old_ci[0]) + (400 * Utility::pow<3>(old_ci[0])) / 3 - 200 *
    //               Utility::pow<2>(old_ci[1]) - (400 * Utility::pow<3>(old_ci[1])) / 3 - 4279.99;
    // three terms of Taylor
    abs_err[0] = 800 * new_ci[0] - 800 * new_ci[1] +
                 400 * new_ci[0] * (Utility::pow<2>(new_ci[0] - 1) - new_ci[0] + 2) -
                 400 * new_ci[1] * (Utility::pow<2>(new_ci[1] - 1) - new_ci[1] + 2) -
                 200 * Utility::pow<2>(new_ci[0] - 1) + (400 * Utility::pow<3>(new_ci[0] - 1)) / 3 +
                 200 * Utility::pow<2>(new_ci[1] - 1) - (400 * Utility::pow<3>(new_ci[1] - 1)) / 3 +
                 400 * (new_ci[0] - 1) * (Utility::pow<2>(new_ci[0]) + new_ci[0] + 1) -
                 400 * (new_ci[1] - 1) * (Utility::pow<2>(new_ci[1]) + new_ci[1] + 1) +
                 200 * Utility::pow<2>(new_ci[0]) + (400 * Utility::pow<3>(new_ci[0])) / 3 -
                 200 * Utility::pow<2>(new_ci[1]) - (400 * Utility::pow<3>(new_ci[1])) / 3 -
                 4279.99;
    abs_err[1] = _c[_qp] - new_ci[1] * _h[_qp] + new_ci[0] * (_h[_qp] - 1);
    abs_err_norm = std::sqrt(Utility::pow<2>(abs_err[0]) + Utility::pow<2>(abs_err[1]));

    // abs_err[0] = 400 * log(1 - new_ci[1]) - 400 * log(1 - new_ci[0]) + 400 * log(new_ci[0]) -
    //              400 * log(new_ci[1]) - 67;
    // abs_err[1] = _c[_qp] - new_ci[1] * _h[_qp] + new_ci[0] * (_h[_qp] - 1);
    // abs_err_norm = std::sqrt(Utility::pow<2>(abs_err[0]) + Utility::pow<2>(abs_err[1]));

    rel_err_norm = std::abs(abs_err_norm / init_err_norm);

    std::cout << "Newton iteration loop " << nloop << ", the absolute error norm is "
              << abs_err_norm << ", and the relative error norm is " << rel_err_norm << std::endl;

    old_ci = new_ci; // update ci

    if (abs_err_norm < _abs_tol)
      break;
    else if (rel_err_norm < _rel_tol)
      break;
  }

  std::cout << "c1 is " << new_ci[0] << ", and c2 is " << new_ci[1] << std::endl;

  n_qp += 1;

  std::cout << "qp point " << n_qp << std::endl;

  for (unsigned int i = 0; i < 2; ++i)
    (*_ci_prop[i])[_qp] = new_ci[i];
}
