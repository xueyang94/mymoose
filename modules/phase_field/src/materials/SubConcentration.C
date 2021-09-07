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

registerMooseObject("PhaseFieldApp", SubConcentration);

InputParameters
SubConcentration::validParams()
{
  // InputParameters params = Material::validParams();
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

{
}

void
SubConcentration::computeQpProperties()
{
  Real n = _eta[_qp];

  // declare and initialize the old ci inside Newton iteration
  // std::vector<Real> old_ci_Newton(2);
  // old_ci_Newton[0] = _c1_old[_qp];
  // old_ci_Newton[1] = _c2_old[_qp];
  // old_ci_Newton[0] = 0.6;
  // old_ci_Newton[1] = 0.4;
  // old_ci_Newton[0] = 0.6;
  // old_ci_Newton[1] = 0.1;
  _c1[_qp] = 0.4;
  _c2[_qp] = 0.6;

  // std::cout << "first_f1 is " << _first_df1[_qp] << std::endl;
  // std::cout << "first_f2 is " << _first_df2[_qp] << std::endl;

  // declare the error vectors and norms of Newton iteration
  std::vector<Real> init_err(2);
  std::vector<Real> abs_err(2);

  Real init_err_norm;
  Real abs_err_norm;
  Real rel_err_norm;

  // compute the initial error norm
  init_err[0] = _first_df1[_qp] - _first_df2[_qp];
  init_err[1] = _c[_qp] - _c2[_qp] * _h[_qp] + _c1[_qp] * (_h[_qp] - 1);
  init_err_norm = std::sqrt(Utility::pow<2>(init_err[0]) + Utility::pow<2>(init_err[1]));

  // Newton iteration
  for (unsigned int nloop = 0; nloop < _maxiter; ++nloop)
  {
    _f1.computePropertiesAtQp(_qp);
    _f2.computePropertiesAtQp(_qp);

    // compute eqn1 and eqn2
    Real eqn1 = _first_df1[_qp] - _first_df2[_qp];
    Real eqn2 = _c[_qp] - _c2[_qp] * _h[_qp] + _c1[_qp] * (_h[_qp] - 1);

    // terms used in the determinant D
    Real deqn1dc1 = _second_df1[_qp];
    Real deqn1dc2 = -_second_df2[_qp];
    Real deqn2dc1 = _h[_qp] - 1;
    Real deqn2dc2 = -_h[_qp];

    // the determinant of the Jacobian matrix
    Real D = deqn1dc1 * deqn2dc2 - deqn1dc2 * deqn2dc1;

    // compute the update of ci (is inv(J)*f_vec)
    Real update[2];

    update[0] = 1 / D * (eqn1 * deqn2dc2 - eqn2 * deqn1dc2);

    update[1] = 1 / D * (-eqn1 * deqn2dc1 + eqn2 * deqn1dc1);

    // compute c1 and c2
    _c1[_qp] -= update[0];

    _c2[_qp] -= update[1];

    // compute the updated absolute Newton error
    abs_err[0] = _first_df1[_qp] - _first_df2[_qp];
    abs_err[1] = _c[_qp] - _c2[_qp] * _h[_qp] + _c1[_qp] * (_h[_qp] - 1);
    abs_err_norm = std::sqrt(Utility::pow<2>(abs_err[0]) + Utility::pow<2>(abs_err[1]));

    // compute the relative Newton error
    rel_err_norm = std::abs(abs_err_norm / init_err_norm);

    // std::cout << "c1 is " << _c1[_qp] << std::endl;
    // std::cout << "c2 is " << _c2[_qp] << std::endl;
    // std::cout << "absolute error is " << abs_err_norm << std::endl;
    // std::cout << "relative error is " << rel_err_norm << std::endl;

    // Newton iteration convergence criterion
    if ((abs_err_norm < _abs_tol) || (rel_err_norm < _rel_tol))
      break;
    // else if (rel_err_norm < _rel_tol)
    //   break;

    if (nloop == (_maxiter - 1))
      mooseError("The SubConcentration Newton iteration exceeds the max iteration.");

    // // update old ci
    // old_ci_Newton[0] = _c1[_qp];
    // old_ci_Newton[1] = _c2[_qp];
  }

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
