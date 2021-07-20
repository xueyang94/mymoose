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
#include <cmath>
#include "libmesh/fparser_ad.hh"

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
  // params.addRequiredParam<MaterialPropertyName>("c_material",
  //                                               "The global concentration as a material
  //                                               property");
  // params.addRequiredParam<MaterialPropertyName>("eta_material",
  //                                               "The order parameter as a material property");
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
  params.addRequiredParam<Real>(
      "plog_tol_value", "The maximum value to start using the Taylor expansion of log functions");
  return params;
}

SubConcentration::SubConcentration(const InputParameters & parameters)
  : Material(parameters),
    _c(coupledValue("global_c")),
    // _c_mat(declareProperty<Real>("c_material")),
    // _c_mat_old(getMaterialPropertyOld<Real>("c_material")), // old
    _eta(coupledValue("eta_name")),
    // _eta_mat(declareProperty<Real>("eta_material")),
    // _eta_mat_old(getMaterialPropertyOld<Real>("eta_material")), // old
    _h(getMaterialProperty<Real>("h_name")),
    _c1(declareProperty<Real>("c1_name")),
    _c2(declareProperty<Real>("c2_name")),
    _c1_old(getMaterialPropertyOld<Real>("c1_name")), // old
    _c2_old(getMaterialPropertyOld<Real>("c2_name")), // old
    _abs_tol(getParam<Real>("absolute_tol_value")),
    _rel_tol(getParam<Real>("relative_tol_value")),
    _maxiter(getParam<Real>("max_iteration")),
    _dc1dc(declareProperty<Real>("dc1dc_name")),
    _dc1deta(declareProperty<Real>("dc1deta_name")),
    _dc2dc(declareProperty<Real>("dc2dc_name")),
    _dc2deta(declareProperty<Real>("dc2deta_name")),
    _tol(getParam<Real>("plog_tol_value"))

{
}

void
SubConcentration::initQpStatefulProperties()
{
  // init the ci property (this will become _c1_old and _c2_old in the first call of
  // computeProperties)
  _c1[_qp] = 0.6;
  _c2[_qp] = 0.1;
}

// Taylor expansion defination of plog()
Real
tlog(Real x, Real eps)
{
  return log(eps) + 1 / eps * (x - eps) - 1 / (2 * eps * eps) * (x - eps) * (x - eps) +
         1 / (3 * eps * eps * eps) * (x - eps) * (x - eps) * (x - eps);
}

void
SubConcentration::computeQpProperties()
{
  FunctionParserADBase<Real> fparser;

  // std::string f1 = "x^2 + 3*a^3";
  // std::string f2 = "7*a^2";

  // // Parse the input expression into bytecode
  // fparser.Parse(f1, "x,a");
  //
  // // transform F -> dF/dx
  // fparser.AutoDiff("x");
  //
  // // run optimizer to simplify the derivative
  // fparser.Optimize();
  //
  // // evaluate the derivative (method from FParserBase<Real>)
  // Real params[2] = {0.1, 1.7};
  // std::cout << fparser.Eval(params) << std::endl;

  // fparser.Parse(f2, "a");
  // fparser.AutoDiff("a");
  // fparser.Optimize();
  // // Real n = 0.1;
  // Real m = 1.7;
  // Real pp[1] = {m};
  // std::cout << fparser.Eval(pp) << std::endl;

  /////////////////////////////////////////////////////////////////////////////////

  Real n = _eta[_qp];

  // old ci inside Newton iteration
  std::vector<Real> old_ci_Newton(2);
  old_ci_Newton[0] = _c1_old[_qp];
  old_ci_Newton[1] = _c2_old[_qp];

  std::string f1;
  std::string f2;

  // f1 and f2 equations based on the range of c1 and c2
  if (old_ci_Newton[0] < _tol && old_ci_Newton[1] < _tol)
  {

    f1 = "20*c1 + 300*(1 - c1) + 400*(c1*tlog(c1) + (1 - c1)*log(1 - c1))";
    f2 = "4000*c2 + 0.01*(1 - c2) + 400*(c2*tlog(c2) + (1 - c2)*log(1 - c2))";
  }
  else if (old_ci_Newton[0] >= _tol && old_ci_Newton[0] < 1 && old_ci_Newton[1] < _tol)
  {

    f1 = "20*c1 + 300*(1 - c1) + 400*(c1*log(c1) + (1 - c1)*log(1 - c1))";
    f2 = "4000*c2 + 0.01*(1 - c2) + 400*(c2*tlog(c2) + (1 - c2)*log(1 - c2))";
  }
  else if (old_ci_Newton[0] >= 1 && old_ci_Newton[1] < _tol)
  {

    f1 = "20*c1 + 300*(1 - c1) + 400*(c1*log(c1) + (1 - c1)*tlog(1 - c1))";
    f2 = "4000*c2 + 0.01*(1 - c2) + 400*(c2*tlog(c2) + (1 - c2)*log(1 - c2))";
  }
  else if (old_ci_Newton[0] < _tol && old_ci_Newton[1] >= _tol && old_ci_Newton[1] < 1)
  {

    f1 = "20*c1 + 300*(1 - c1) + 400*(c1*tlog(c1) + (1 - c1)*log(1 - c1))";
    f2 = "4000*c2 + 0.01*(1 - c2) + 400*(c2*log(c2) + (1 - c2)*log(1 - c2))";
  }
  else if (old_ci_Newton[0] >= _tol && old_ci_Newton[0] < 1 && old_ci_Newton[1] >= _tol &&
           old_ci_Newton[1] < 1)
  {

    f1 = "20*c1 + 300*(1 - c1) + 400*(c1*log(c1) + (1 - c1)*log(1 - c1))";
    f2 = "4000*c2 + 0.01*(1 - c2) + 400*(c2*log(c2) + (1 - c2)*log(1 - c2))";
  }
  else if (old_ci_Newton[0] >= 1 && old_ci_Newton[1] >= _tol && old_ci_Newton[1] < 1)
  {

    f1 = "20*c1 + 300*(1 - c1) + 400*(c1*log(c1) + (1 - c1)*tlog(1 - c1))";
    f2 = "4000*c2 + 0.01*(1 - c2) + 400*(c2*log(c2) + (1 - c2)*log(1 - c2))";
  }
  else if (old_ci_Newton[0] < _tol && old_ci_Newton[1] >= 1)
  {

    f1 = "20*c1 + 300*(1 - c1) + 400*(c1*tlog(c1) + (1 - c1)*log(1 - c1))";
    f2 = "4000*c2 + 0.01*(1 - c2) + 400*(c2*log(c2) + (1 - c2)*tlog(1 - c2))";
  }
  else if (old_ci_Newton[0] >= _tol && old_ci_Newton[0] < 1 && old_ci_Newton[1] >= 1)
  {

    f1 = "20*c1 + 300*(1 - c1) + 400*(c1*log(c1) + (1 - c1)*log(1 - c1))";
    f2 = "4000*c2 + 0.01*(1 - c2) + 400*(c2*log(c2) + (1 - c2)*tlog(1 - c2))";
  }
  else if (old_ci_Newton[0] >= 1 && old_ci_Newton[1] >= 1)
  {

    f1 = "20*c1 + 300*(1 - c1) + 400*(c1*log(c1) + (1 - c1)*tlog(1 - c1))";
    f2 = "4000*c2 + 0.01*(1 - c2) + 400*(c2*log(c2) + (1 - c2)*tlog(1 - c2))";
  }

  Real params[1]; // declare the params used in substitution of symbolic functions

  // compute df1dc1_init
  fparser.Parse(f1, "c1");
  fparser.AutoDiff("c1");
  fparser.Optimize();
  params[0] = {old_ci_Newton[0]};
  Real df1dc1_init = fparser.Eval(params);

  // compute df2dc2_init
  fparser.Parse(f2, "c2");
  fparser.AutoDiff("c2");
  fparser.Optimize();
  params[0] = {old_ci_Newton[1]};
  Real df2dc2_init = fparser.Eval(params);

  // declarations of the error vectors and norms of Newton iteration
  std::vector<Real> init_err(2);
  std::vector<Real> abs_err(2);

  Real init_err_norm;
  Real abs_err_norm;
  Real rel_err_norm;

  // compute the initial error norm
  init_err[0] = df1dc1_init - df2dc2_init;
  init_err[1] = _c[_qp] - old_ci_Newton[1] * _h[_qp] + old_ci_Newton[0] * (_h[_qp] - 1);
  init_err_norm = std::sqrt(Utility::pow<2>(init_err[0]) + Utility::pow<2>(init_err[1]));

  Real count;
  count = 0;

  // declare and initialize the first and second derivatives of the phase energies
  Real df1dc1_Jacob;
  Real d2c1_Jacob;
  Real df2dc2_Jacob;
  Real d2c2_Jacob;

  for (unsigned int nloop = 0; nloop < _maxiter; ++nloop)
  {

    // compute df1dc1_Jacob
    fparser.Parse(f1, "c1");
    fparser.AutoDiff("c1");
    fparser.Optimize();
    params[0] = {old_ci_Newton[0]};
    df1dc1_Jacob = fparser.Eval(params);

    // compute second derivative d2c1_Jacob
    fparser.Parse(f1, "c1");
    fparser.AutoDiff("c1");
    fparser.AutoDiff("c1");
    fparser.Optimize();
    params[0] = {old_ci_Newton[0]};
    d2c1_Jacob = fparser.Eval(params);

    // compute df2dc2_Jacob
    fparser.Parse(f2, "c2");
    fparser.AutoDiff("c2");
    fparser.Optimize();
    params[0] = {old_ci_Newton[1]};
    df2dc2_Jacob = fparser.Eval(params);

    // compute second derivative d2c2_Jacob
    fparser.Parse(f2, "c2");
    fparser.AutoDiff("c2");
    fparser.AutoDiff("c2");
    fparser.Optimize();
    params[0] = {old_ci_Newton[1]};
    d2c2_Jacob = fparser.Eval(params);

    // compute eqn1 and eqn2
    Real eqn1;
    Real eqn2;

    eqn1 = df1dc1_Jacob - df2dc2_Jacob;
    eqn2 = _c[_qp] - old_ci_Newton[1] * _h[_qp] + old_ci_Newton[0] * (_h[_qp] - 1);

    // terms used in the determinant
    Real deqn1dc1 = d2c1_Jacob;
    Real deqn1dc2 = -d2c2_Jacob;
    Real deqn2dc1 = _h[_qp] - 1;
    Real deqn2dc2 = -_h[_qp];

    // The determinant of the Jacobian matrix
    Real D = 1 / (deqn1dc1 * deqn2dc2 - deqn1dc2 * deqn2dc1);

    // compute the update of ci (inv(J)*f_vec)
    Real update[2];

    update[0] = D * (eqn1 * deqn2dc2 - eqn2 * deqn1dc2);

    update[1] = D * (-eqn1 * deqn2dc1 + eqn2 * deqn1dc1);

    // compute c1 and c2
    _c1[_qp] = old_ci_Newton[0] - update[0];

    _c2[_qp] = old_ci_Newton[1] - update[1];

    // compute df1dc1_new and df2dc2_new for calculating the updated absolute error
    fparser.Parse(f1, "c1");
    fparser.AutoDiff("c1");
    fparser.Optimize();
    params[0] = {_c1[_qp]};
    Real df1dc1_new = fparser.Eval(params);

    fparser.Parse(f2, "c2");
    fparser.AutoDiff("c2");
    fparser.Optimize();
    params[0] = {_c2[_qp]};
    Real df2dc2_new = fparser.Eval(params);

    // compute the absolute Newton error
    abs_err[0] = df1dc1_new - df2dc2_new;
    abs_err[1] = _c[_qp] - _c2[_qp] * _h[_qp] + _c1[_qp] * (_h[_qp] - 1);
    abs_err_norm = std::sqrt(Utility::pow<2>(abs_err[0]) + Utility::pow<2>(abs_err[1]));

    // compute the relative Newton error
    rel_err_norm = std::abs(abs_err_norm / init_err_norm);

    // std::cout << "current count is " << count << "
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    //           << std::endl;
    count += 1;

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

  // compute dc1dc, dc2dc, dc1deta, and dc2deta
  // ////////////////////////////////////////////////////////////////////////////////////////
  _dc1dc[_qp] = d2c2_Jacob / (d2c2_Jacob + _h[_qp] * (d2c1_Jacob - d2c2_Jacob));
  _dc2dc[_qp] = d2c1_Jacob / (d2c2_Jacob + _h[_qp] * (d2c1_Jacob - d2c2_Jacob));
  _dc1deta[_qp] = d2c2_Jacob * (_c1[_qp] - _c2[_qp]) * (30.0 * n * n * (n * n - 2.0 * n + 1.0)) /
                  (d2c2_Jacob + _h[_qp] * (d2c1_Jacob - d2c2_Jacob));
  _dc2deta[_qp] = d2c1_Jacob * (_c1[_qp] - _c2[_qp]) * (30.0 * n * n * (n * n - 2.0 * n + 1.0)) /
                  (d2c2_Jacob + _h[_qp] * (d2c1_Jacob - d2c2_Jacob));
}
