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
  InputParameters params = Material::validParams();
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
  params.addRequiredParam<MaterialPropertyName>("df1dc1_name",
                                                "The first derivative of f1 w.r.t. c1");
  params.addRequiredParam<MaterialPropertyName>("df2dc2_name",
                                                "The first derivative of f2 w.r.t. c2");
  params.addRequiredParam<MaterialPropertyName>("d2f1dc1_name",
                                                "The second derivative of f1 w.r.t. c1");
  params.addRequiredParam<MaterialPropertyName>("d2f2dc2_name",
                                                "The second derivative of f2 w.r.t. c2");
  params.addRequiredParam<Real>(
      "plog_tol_value", "The maximum value to start using the Taylor expansion of log functions");
  params.addParam<Real>("dt", "Time step");
  params.addParam<Real>("visco", "Viscosity");
  return params;
}

SubConcentration::SubConcentration(const InputParameters & parameters)
  : Material(parameters),
    // : DerivativeMaterialInterface<Material>(parameters),
    _c(coupledValue("global_c")),
    // _c_name(getVar("global_c", 0)->name()),
    _eta(coupledValue("eta")),
    // _eta_name(getVar("eta", 0)->name()),
    _h(getMaterialProperty<Real>("h_name")),
    _c1(declareProperty<Real>("c1_name")),
    _c2(declareProperty<Real>("c2_name")),
    _c1_old(getMaterialPropertyOld<Real>("c1_name")), // old
    _c2_old(getMaterialPropertyOld<Real>("c2_name")), // old
    _abs_tol(getParam<Real>("absolute_tol_value")),
    _rel_tol(getParam<Real>("relative_tol_value")),
    _maxiter(getParam<Real>("max_iteration")),
    _dc1dc(declareProperty<Real>("dc1dc_name")),
    _dc2dc(declareProperty<Real>("dc2dc_name")),
    _dc1deta(declareProperty<Real>("dc1deta_name")),
    _dc2deta(declareProperty<Real>("dc2deta_name")),
    // _dc1dc(declarePropertyDerivative<Real>("c1_name", _c_name)),
    // _dc2dc(declarePropertyDerivative<Real>("c2_name", _c_name)),
    // _dc1deta(declarePropertyDerivative<Real>("c1_name", _eta_name)),
    // _dc2deta(declarePropertyDerivative<Real>("c2_name", _eta_name)),
    _first_df1(declareProperty<Real>("df1dc1_name")),
    _first_df2(declareProperty<Real>("df2dc2_name")),
    _second_df1(declareProperty<Real>("d2f1dc1_name")),
    _second_df2(declareProperty<Real>("d2f2dc2_name")),
    _dt(getParam<Real>("dt")),
    _visco(getParam<Real>("visco"))

{
  _fparser1 = std::make_unique<FunctionParserADBase<Real>>();
  _fparser2 = std::make_unique<FunctionParserADBase<Real>>();
  _fparser3 = std::make_unique<FunctionParserADBase<Real>>();
  _fparser4 = std::make_unique<FunctionParserADBase<Real>>();

  // declare bulk energies f1 and f2
  // std::string f1 = "2*c1 + 30*(1 - c1) + 400*(c1*plog(c1, 1e-4) + (1 - c1)*plog(1 - c1, 1e-4))";
  // std::string f2 = "40*c2 + (1 - c2) + 400*(c2*plog(c2, 1e-4) + (1 - c2)*plog(1 - c2, 1e-4))";
  std::string f1 = "20*c1 + 300*(1 - c1) + 400*(c1*plog(c1, 1e-4) + (1 - c1)*plog(1 - c1,1e-4))";
  std::string f2 = "4000*c2 + 0.01*(1 - c2) + 400*(c2*plog(c2, 1e-4) + (1 - c2)*plog(1 - c2,1e-4))";
  // std::string f1 = "1e2*(c1 - 0.3)^2 - 10";
  // std::string f2 = "1e2*(c2 - 0.7)^2";

  // parsed function of df1/dc1
  _fparser1->Parse(f1, "c1");
  _fparser1->AutoDiff("c1");
  _fparser1->Optimize();

  // parsed function of df2/dc2
  _fparser2->Parse(f2, "c2");
  _fparser2->AutoDiff("c2");
  _fparser2->Optimize();

  // parsed function of second derivative of f1 w.r.t. c1
  _fparser3->Parse(f1, "c1");
  _fparser3->AutoDiff("c1");
  _fparser3->AutoDiff("c1");
  _fparser3->Optimize();

  // parsed function of second derivative of f2 w.r.t. c2
  _fparser4->Parse(f2, "c2");
  _fparser4->AutoDiff("c2");
  _fparser4->AutoDiff("c2");
  _fparser4->Optimize();
}

void
SubConcentration::initQpStatefulProperties()
{
  // init the ci property (this will become _c1_old and _c2_old in the first call of
  // computeProperties)
  _c1[_qp] = 0.6;
  _c2[_qp] = 0.1;
  // _c1[_qp] = 0.4;
  // _c2[_qp] = 0.6;
}

void
SubConcentration::computeQpProperties()
{
  FunctionParserADBase<Real> fparser;

  Real n = _eta[_qp];

  // declare and initialize the old ci inside Newton iteration
  std::vector<Real> old_ci_Newton(2);
  old_ci_Newton[0] = _c1_old[_qp];
  old_ci_Newton[1] = _c2_old[_qp];
  // old_ci_Newton[0] = 0.6;
  // old_ci_Newton[1] = 0.4;
  // old_ci_Newton[0] = 0.6;
  // old_ci_Newton[1] = 0.1;

  // declare the params used in substitution of symbolic functions fparser.Eval()
  double params;
  double * p = &params;

  // compute df1dc1_init and df2dc2_init for computing the initial error
  params = old_ci_Newton[0];
  Real df1dc1_init = _fparser1->Eval(p);

  params = old_ci_Newton[1];
  Real df2dc2_init = _fparser2->Eval(p);

  // declare the error vectors and norms of Newton iteration
  std::vector<Real> init_err(2);
  std::vector<Real> abs_err(2);

  Real init_err_norm;
  Real abs_err_norm;
  Real rel_err_norm;

  // compute the initial error norm
  init_err[0] = df1dc1_init - df2dc2_init;
  init_err[1] = _c[_qp] - old_ci_Newton[1] * _h[_qp] + old_ci_Newton[0] * (_h[_qp] - 1);
  init_err_norm = std::sqrt(Utility::pow<2>(init_err[0]) + Utility::pow<2>(init_err[1]));

  // Newton iteration
  for (unsigned int nloop = 0; nloop < _maxiter; ++nloop)
  {

    // compute first_df1 in eqn1
    params = old_ci_Newton[0];
    _first_df1[_qp] = _fparser1->Eval(p);

    // compute second derivative second_df1 in determinant D
    params = old_ci_Newton[0];
    _second_df1[_qp] = _fparser3->Eval(p);

    // compute first_df2 in eqn1
    params = old_ci_Newton[1];
    _first_df2[_qp] = _fparser2->Eval(p);

    // compute second derivative second_df2 in determinant D
    params = old_ci_Newton[1];
    _second_df2[_qp] = _fparser4->Eval(p);

    // compute eqn1 and eqn2
    Real eqn1 = _first_df1[_qp] - _first_df2[_qp];
    Real eqn2 = _c[_qp] - old_ci_Newton[1] * _h[_qp] + old_ci_Newton[0] * (_h[_qp] - 1);

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
    _c1[_qp] = old_ci_Newton[0] - update[0];

    _c2[_qp] = old_ci_Newton[1] - update[1];

    // compute df1dc1_new and df2dc2_new for calculating the updated absolute error
    params = _c1[_qp];
    Real df1dc1_new = _fparser1->Eval(p);

    params = _c2[_qp];
    Real df2dc2_new = _fparser2->Eval(p);

    // compute the updated absolute Newton error
    abs_err[0] = df1dc1_new - df2dc2_new;
    abs_err[1] = _c[_qp] - _c2[_qp] * _h[_qp] + _c1[_qp] * (_h[_qp] - 1);
    abs_err_norm = std::sqrt(Utility::pow<2>(abs_err[0]) + Utility::pow<2>(abs_err[1]));

    // compute the relative Newton error
    rel_err_norm = std::abs(abs_err_norm / init_err_norm);

    // std::cout << "nloop is " << nloop << std::endl;
    // std::cout << "initial error norm is " << init_err_norm << std::endl;
    // std::cout << "absolute error norm is " << abs_err_norm << std::endl;
    // std::cout << "relative error norm is " << rel_err_norm << std::endl;

    // Newton iteration convergence criterion
    if (abs_err_norm < _abs_tol)
      break;
    // else if (rel_err_norm < _rel_tol)
    //   break;

    if (nloop == (_maxiter - 1))
      mooseError("The SubConcentration Newton iteration exceeds the max iteration.");

    // update old ci
    old_ci_Newton[0] = _c1[_qp];
    old_ci_Newton[1] = _c2[_qp];
  }

  // use viscosity to reduce ci interfacial oscillation
  _c1[_qp] = (_c1[_qp] + _visco / _dt * _c1_old[_qp]) / (1 + _visco / _dt);
  _c2[_qp] = (_c2[_qp] + _visco / _dt * _c2_old[_qp]) / (1 + _visco / _dt);

  // compute the updated first and second derivatives in the kernel R and chain rule
  // derivatives
  params = _c1[_qp];
  _first_df1[_qp] = _fparser1->Eval(p);

  params = _c2[_qp];
  _first_df2[_qp] = _fparser2->Eval(p);

  params = _c1[_qp];
  _second_df1[_qp] = _fparser3->Eval(p);

  params = _c2[_qp];
  _second_df2[_qp] = _fparser4->Eval(p);

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
