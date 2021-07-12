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
  params.addRequiredParam<MaterialPropertyName>("c_material",
                                                "The global concentration as a material property");
  params.addRequiredParam<MaterialPropertyName>("eta_material",
                                                "The order parameter as a material property");
  params.addRequiredParam<MaterialPropertyName>("h_name", "Name of the switching function");
  params.addRequiredParam<MaterialPropertyName>("c1_name", "Name of the first phase concentration");
  params.addRequiredParam<MaterialPropertyName>("c2_name",
                                                "Name of the second phase concentration");
  params.addParam<Real>("absolute_tol_value", 1e-9, "Absolute tolerance of the Newton iteration");
  params.addParam<Real>("relative_tol_value", 1e-9, "Relative tolerance of the Newton iteration");
  params.addParam<Real>("max_iteration", 100, "The maximum number of Newton iterations");
  params.addRequiredParam<MaterialPropertyName>("dc1dc",
                                                "The approximated first derivative of c1 w.r.t. c");
  params.addRequiredParam<MaterialPropertyName>(
      "dc1deta", "The approximated first derivative of c1 w.r.t. eta");
  params.addRequiredParam<MaterialPropertyName>("dc2dc",
                                                "The approximated first derivative of c2 w.r.t. c");
  params.addRequiredParam<MaterialPropertyName>(
      "dc2deta", "The approximated first derivative of c2 w.r.t. eta");
  return params;
}

SubConcentration::SubConcentration(const InputParameters & parameters)
  : Material(parameters),
    _c(coupledValue("global_c")),
    _c_mat(getMaterialProperty<Real>("c_material")),
    _c_mat_old(getMaterialPropertyOld<Real>("c_material")),
    _eta_mat(getMaterialProperty<Real>("eta_material")),
    _eta_mat_old(getMaterialPropertyOld<Real>("eta_material")),
    _h(getMaterialProperty<Real>("h_name")),
    _c1(declareProperty<Real>("c1_name")),
    _c2(declareProperty<Real>("c2_name")),
    _c1_old(getMaterialPropertyOld<Real>("c1_name")),
    _c2_old(getMaterialPropertyOld<Real>("c2_name")),
    _abs_tol(getParam<Real>("absolute_tol_value")),
    _rel_tol(getParam<Real>("relative_tol_value")),
    _maxiter(getParam<Real>("max_iteration")),
    _dc1dc(declareProperty<Real>("dc1dc")),
    _dc1deta(declareProperty<Real>("dc1deta")),
    _dc2dc(declareProperty<Real>("dc2dc")),
    _dc2deta(declareProperty<Real>("dc2deta"))

{
}

void
SubConcentration::initQpStatefulProperties()
{
  _c1[_qp] = 0.6;
  _c2[_qp] = 0.4;

  _c_mat[_qp] = (1 - _h[_qp]) * _c1[_qp] + _h[_qp] * c2[_qp];

  _eta_mat[_qp] = _h[_qp];
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

  init_err[0] = 400 * log(1 - old_ci_Newton[1]) - 400 * log(1 - old_cold_ci_Newton[0]) +
                400 * log(old_ci_Newton[0]) - 400 * log(old_ci_Newton[1]) - 67;
  init_err[1] = _c[_qp] - old_ci_Newton[1] * _h[_qp] + old_ci_Newton[0] * (_h[_qp] - 1);
  init_err_norm = std::sqrt(Utility::pow<2>(init_err[0]) + Utility::pow<2>(init_err[1]));

  for (unsigned int nloop = 0; nloop < _maxiter; ++nloop)
  {
    // new_ci[0] = old_ci_Newton[0] -
    //             ((old_ci_Newton[0] * (old_ci_Newton[0] - 1) *
    //               (_c[_qp] - old_ci_Newton[1] * _h[_qp] + old_ci_Newton[0] * (_h[_qp] - 1))) /
    //                  (old_ci_Newton[0] - old_ci_Newton[0] * _h[_qp] + old_ci_Newton[1] * _h[_qp]
    //                  +
    //                   Utility::pow<2>(old_ci_Newton[0]) * _h[_qp] -
    //                   Utility::pow<2>(old_ci_Newton[1]) * _h[_qp]
    //                   - Utility::pow<2>(old_ci_Newton[0])) -
    //              (old_ci_Newton[0] * old_ci_Newton[1] * _h[_qp] * (old_ci_Newton[0] - 1) *
    //              (old_ci_Newton[1] - 1) *
    //               (400 * log(1 - old_ci_Newton[0]) - 400 * log(1 - old_ci_Newton[1]) - 400 *
    //               log(old_ci_Newton[0]) +
    //                400 * log(old_ci_Newton[1]) + 4279.99)) /
    //                  (400 * (old_ci_Newton[0] - old_ci_Newton[0] * _h[_qp] + old_ci_Newton[1] *
    //                  _h[_qp] +
    //                          Utility::pow<2>(old_ci_Newton[0]) * _h[_qp] -
    //                          Utility::pow<2>(old_ci_Newton[1]) * _h[_qp] -
    //                          Utility::pow<2>(old_ci_Newton[0]))));

    _c1[_qp] =
        old_ci_Newton[0] -
        (-(840000 * _c[_qp] - 840000 * old_ci_Newton[0] + 1283997 * _h[_qp] -
           960000 * _c[_qp] * old_ci_Newton[1] + 960000 * old_ci_Newton[0] * old_ci_Newton[1] +
           960000 * _c[_qp] * Utility::pow<2>(old_ci_Newton[1]) -
           960000 * old_ci_Newton[0] * Utility::pow<2>(old_ci_Newton[1]) +
           480000 * Utility::pow<2>(old_ci_Newton[0]) * _h[_qp] -
           320000 * Utility::pow<3>(old_ci_Newton[0]) * _h[_qp] +
           480000 * Utility::pow<2>(old_ci_Newton[1]) * _h[_qp] -
           640000 * Utility::pow<3>(old_ci_Newton[1]) * _h[_qp] -
           960000 * old_ci_Newton[0] * old_ci_Newton[1] * _h[_qp] +
           960000 * old_ci_Newton[0] * Utility::pow<2>(old_ci_Newton[1]) * _h[_qp]) /
         (120000 * (8 * old_ci_Newton[1] * _h[_qp] - 8 * old_ci_Newton[0] * _h[_qp] -
                    8 * old_ci_Newton[1] + 8 * Utility::pow<2>(old_ci_Newton[0]) * _h[_qp] -
                    8 * Utility::pow<2>(old_ci_Newton[1]) * _h[_qp] +
                    8 * Utility::pow<2>(old_ci_Newton[1]) + 7)));

    // new_ci[1] = old_ci_Newton[1] -
    //             ((old_ci_Newton[1] * (old_ci_Newton[1] - 1) *
    //               (_c[_qp] - old_ci_Newton[1] * _h[_qp] + old_ci_Newton[0] * (_h[_qp] - 1))) /
    //                  (old_ci_Newton[0] - old_ci_Newton[0] * _h[_qp] + old_ci_Newton[1] * _h[_qp]
    //                  +
    //                   Utility::pow<2>(old_ci_Newton[0]) * _h[_qp] -
    //                   Utility::pow<2>(old_ci_Newton[1]) * _h[_qp]
    //                   - Utility::pow<2>(old_ci_Newton[0])) -
    //              (old_ci_Newton[0] * old_ci_Newton[1] * (old_ci_Newton[0] - 1) *
    //              (old_ci_Newton[1] - 1) * (_h[_qp] - 1) *
    //               (400 * log(1 - old_ci_Newton[0]) - 400 * log(1 - old_ci_Newton[1]) - 400 *
    //               log(old_ci_Newton[0]) +
    //                400 * log(old_ci_Newton[1]) + 4279.99)) /
    //                  (400 * (old_ci_Newton[0] - old_ci_Newton[0] * _h[_qp] + old_ci_Newton[1] *
    //                  _h[_qp] +
    //                          Utility::pow<2>(old_ci_Newton[0]) * _h[_qp] -
    //                          Utility::pow<2>(old_ci_Newton[1]) * _h[_qp] -
    //                          Utility::pow<2>(old_ci_Newton[0]))));

    _c2[_qp] =
        old_ci_Newton[1] -
        (-(840000 * _c[_qp] - 840000 * old_ci_Newton[1] + 1283997 * _h[_qp] -
           960000 * _c[_qp] * old_ci_Newton[0] +
           960000 * _c[_qp] * Utility::pow<2>(old_ci_Newton[0]) -
           480000 * Utility::pow<2>(old_ci_Newton[0]) * _h[_qp] +
           640000 * Utility::pow<3>(old_ci_Newton[0]) * _h[_qp] -
           480000 * Utility::pow<2>(old_ci_Newton[1]) * _h[_qp] +
           320000 * Utility::pow<3>(old_ci_Newton[1]) * _h[_qp] +
           480000 * Utility::pow<2>(old_ci_Newton[0]) - 640000 * Utility::pow<3>(old_ci_Newton[0]) +
           480000 * Utility::pow<2>(old_ci_Newton[1]) - 320000 * Utility::pow<3>(old_ci_Newton[1]) +
           960000 * old_ci_Newton[0] * old_ci_Newton[1] * _h[_qp] -
           960000 * Utility::pow<2>(old_ci_Newton[0]) * old_ci_Newton[1] * _h[_qp] - 1283997) /
         (120000 * (8 * old_ci_Newton[1] * _h[_qp] - 8 * old_ci_Newton[0] * _h[_qp] -
                    8 * old_ci_Newton[1] + 8 * Utility::pow<2>(old_ci_Newton[0]) * _h[_qp] -
                    8 * Utility::pow<2>(old_ci_Newton[1]) * _h[_qp] +
                    8 * Utility::pow<2>(old_ci_Newton[1]) + 7)));

    // new_ci[0] = old_ci_Newton[0] -
    //             ((old_ci_Newton[0] * (old_ci_Newton[0] - 1) *
    //               (_c[_qp] - old_ci_Newton[1] * _h[_qp] + old_ci_Newton[0] * (_h[_qp] - 1))) /
    //                  (old_ci_Newton[0] - old_ci_Newton[0] * _h[_qp] + old_ci_Newton[1] * _h[_qp]
    //                  +
    //                   Utility::pow<2>(old_ci_Newton[0]) * _h[_qp] -
    //                   Utility::pow<2>(old_ci_Newton[1]) * _h[_qp]
    //                   - Utility::pow<2>(old_ci_Newton[0])) -
    //              (old_ci_Newton[0] * old_ci_Newton[1] * _h[_qp] * (old_ci_Newton[0] - 1) *
    //              (old_ci_Newton[1] - 1) *
    //               (400 * log(1 - old_ci_Newton[0]) - 400 * log(1 - old_ci_Newton[1]) - 400 *
    //               log(old_ci_Newton[0]) +
    //                400 * log(old_ci_Newton[1]) + 67)) /
    //                  (400 * (old_ci_Newton[0] - old_ci_Newton[0] * _h[_qp] + old_ci_Newton[1] *
    //                  _h[_qp] +
    //                          Utility::pow<2>(old_ci_Newton[0]) * _h[_qp] -
    //                          Utility::pow<2>(old_ci_Newton[1]) * _h[_qp] -
    //                          Utility::pow<2>(old_ci_Newton[0]))));
    //
    // new_ci[1] = old_ci_Newton[1] -
    //             ((old_ci_Newton[1] * (old_ci_Newton[1] - 1) *
    //               (_c[_qp] - old_ci_Newton[1] * _h[_qp] + old_ci_Newton[0] * (_h[_qp] - 1))) /
    //                  (old_ci_Newton[0] - old_ci_Newton[0] * _h[_qp] + old_ci_Newton[1] * _h[_qp]
    //                  +
    //                   Utility::pow<2>(old_ci_Newton[0]) * _h[_qp] -
    //                   Utility::pow<2>(old_ci_Newton[1]) * _h[_qp]
    //                   - Utility::pow<2>(old_ci_Newton[0])) -
    //              (old_ci_Newton[0] * old_ci_Newton[1] * (old_ci_Newton[0] - 1) *
    //              (old_ci_Newton[1] - 1) * (_h[_qp] - 1) *
    //               (400 * log(1 - old_ci_Newton[0]) - 400 * log(1 - old_ci_Newton[1]) - 400 *
    //               log(old_ci_Newton[0]) +
    //                400 * log(old_ci_Newton[1]) + 67)) /
    //                  (400 * (old_ci_Newton[0] - old_ci_Newton[0] * _h[_qp] + old_ci_Newton[1] *
    //                  _h[_qp] +
    //                          Utility::pow<2>(old_ci_Newton[0]) * _h[_qp] -
    //                          Utility::pow<2>(old_ci_Newton[1]) * _h[_qp] -
    //                          Utility::pow<2>(old_ci_Newton[0]))));

    // abs_err[0] = 400 * log(1 - new_ci[1]) - 400 * log(1 - new_ci[0]) + 400 * log(new_ci[0]) -
    //              400 * log(new_ci[1]) - 4279.99;
    // abs_err[1] = _c[_qp] - new_ci[1] * _h[_qp] + new_ci[0] * (_h[_qp] - 1);
    // two terms of Taylor
    // abs_err[0] = 800 * old_ci_Newton[0] - 800 * old_ci_Newton[1] +
    //               400 * old_ci_Newton[0] * (Utility::pow<2>(old_ci_Newton[0] - 1) -
    //               old_ci_Newton[0] + 2) - 400 * old_ci_Newton[1] *
    //               (Utility::pow<2>(old_ci_Newton[1] - 1) - old_ci_Newton[1] + 2) - 200 *
    //               Utility::pow<2>(old_ci_Newton[0] - 1) + (400 * Utility::pow<3>(old_ci_Newton[0]
    //               - 1)) / 3 + 200 * Utility::pow<2>(old_ci_Newton[1] - 1) - (400 *
    //               Utility::pow<3>(old_ci_Newton[1] - 1)) / 3 + 400 * (old_ci_Newton[0] - 1) *
    //               (Utility::pow<2>(old_ci_Newton[0]) + old_ci_Newton[0] + 1)
    //               - 400 * (old_ci_Newton[1] - 1) * (Utility::pow<2>(old_ci_Newton[1]) +
    //               old_ci_Newton[1] + 1) + 200 * Utility::pow<2>(old_ci_Newton[0]) + (400 *
    //               Utility::pow<3>(old_ci_Newton[0])) / 3 - 200 *
    //               Utility::pow<2>(old_ci_Newton[1]) - (400 * Utility::pow<3>(old_ci_Newton[1])) /
    //               3 - 4279.99;

    // // three terms of Taylor
    // abs_err[0] = 800 * new_ci[0] - 800 * new_ci[1] +
    //              400 * new_ci[0] * (Utility::pow<2>(new_ci[0] - 1) - new_ci[0] + 2) -
    //              400 * new_ci[1] * (Utility::pow<2>(new_ci[1] - 1) - new_ci[1] + 2) -
    //              200 * Utility::pow<2>(new_ci[0] - 1) + (400 * Utility::pow<3>(new_ci[0] - 1)) /
    //              3 + 200 * Utility::pow<2>(new_ci[1] - 1) - (400 * Utility::pow<3>(new_ci[1] -
    //              1)) / 3 + 400 * (new_ci[0] - 1) * (Utility::pow<2>(new_ci[0]) + new_ci[0] + 1) -
    //              400 * (new_ci[1] - 1) * (Utility::pow<2>(new_ci[1]) + new_ci[1] + 1) +
    //              200 * Utility::pow<2>(new_ci[0]) + (400 * Utility::pow<3>(new_ci[0])) / 3 -
    //              200 * Utility::pow<2>(new_ci[1]) - (400 * Utility::pow<3>(new_ci[1])) / 3 -
    //              4279.99;
    // abs_err[1] = _c[_qp] - new_ci[1] * _h[_qp] + new_ci[0] * (_h[_qp] - 1);
    // abs_err_norm = std::sqrt(Utility::pow<2>(abs_err[0]) + Utility::pow<2>(abs_err[1]));

    // calculate dc1/dc and dc1/deta
    _dc1dc[_qp] = (_c1_old[_qp] - _c1[_qp]) / (_c_mat_old[_qp] - _c_mat[_qp]);
    _dc1deta[_qp] = (_c1_old[_qp] - _c1[_qp]) / (_eta_mat_old[_qp] - _eta_mat[_qp]);

    _dc2dc[_qp] = (_c2_old[_qp] - _c2[_qp]) / (_c_mat_old[_qp] - _c_mat[_qp]);
    _dc2deta[_qp] = (_c2_old[_qp] - _c2[_qp]) / (_eta_mat_old[_qp] - _eta_mat[_qp]);

    abs_err[0] = 400 * log(1 - _c2[_qp]) - 400 * log(1 - _c1[_qp]) + 400 * log(_c1[_qp]) -
                 400 * log(_c2[_qp]) - 67;
    abs_err[1] = _c[_qp] - _c2[_qp] * _h[_qp] + _c1[_qp] * (_h[_qp] - 1);
    abs_err_norm = std::sqrt(Utility::pow<2>(abs_err[0]) + Utility::pow<2>(abs_err[1]));

    rel_err_norm = std::abs(abs_err_norm / init_err_norm);

    if (abs_err_norm < _abs_tol)
      break;
    else if (rel_err_norm < _rel_tol)
      break;
  }
}
