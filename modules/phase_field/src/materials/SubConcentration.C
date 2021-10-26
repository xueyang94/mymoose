//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SubConcentration.h"
#include "MatrixTools.h"

registerMooseObject("PhaseFieldApp", SubConcentration);

InputParameters
SubConcentration::validParams()
{
  InputParameters params = DerivativeMaterialInterface<Material>::validParams();
  params.addClassDescription(
      "Computes the KKS sub-concentrations by using Newton iteration to solve the equal chemical "
      "potential and concentration conservation equations");
  params.addRequiredCoupledVar("global_c", "The interpolated concentration");
  params.addRequiredCoupledVar("all_etas", "Vector of all order parameters for all phases");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "hj_names", "Names of the switching functions in the same order of the all_etas");
  params.addRequiredParam<std::vector<MaterialPropertyName>>("ci_names", "Phase concentrations");

  params.addRequiredParam<MaterialPropertyName>("c1_name", "c1 name");
  params.addRequiredParam<MaterialPropertyName>("c2_name", "c2 name");
  params.addRequiredParam<MaterialPropertyName>("c3_name", "c3 name");
  // params.addRequiredParam<std::vector<MaterialPropertyName>>("test", "Phase concentrations");

  params.addRequiredParam<std::vector<Real>>("ci_IC",
                                             "Initial values of ci in the same order of ci_names");

  params.addRequiredParam<std::vector<MaterialPropertyName>>("dcidc_names",
                                                             "The first derivative of ci wrt c");

  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "dcidetaj_names",
      "The names of dci/detaj in the order of dc1deta1, dc2deta1, dc3deta1, dc1deta2, dc2deta2, "
      "dc3deta2, etc");

  params.addRequiredParam<MaterialName>("F1_material", "F1");
  params.addRequiredParam<MaterialName>("F2_material", "F2");
  params.addRequiredParam<MaterialName>("F3_material", "F3");
  params.addRequiredParam<MaterialPropertyName>("F1_name", "F1");
  params.addRequiredParam<MaterialPropertyName>("F2_name", "F2");
  params.addRequiredParam<MaterialPropertyName>("F3_name", "F3");

  params.addParam<MaterialPropertyName>(
      "nested_iterations", "The number of nested Newton iterations at each quadrature point");
  params.set<unsigned int>("min_iterations") = 10;
  params.set<unsigned int>("max_iterations") = 1000;
  params.set<Real>("absolute_tolerance") = 1e-13;
  params.set<Real>("relative_tolerance") = 1e-8;
  return params;
}

SubConcentration::SubConcentration(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _c(coupledValue("global_c")),
    _num_eta(coupledComponents("all_etas")),
    _eta_names(coupledNames("all_etas")),
    _hj_names(getParam<std::vector<MaterialPropertyName>>("hj_names")),
    _prop_hj(_num_eta),
    _prop_dhjdetai(_num_eta),

    _ci_names(getParam<std::vector<MaterialPropertyName>>("ci_names")),
    _ci_prop(_num_eta),

    _c1_old(getMaterialPropertyOld<Real>("c1_name")), // old
    _c2_old(getMaterialPropertyOld<Real>("c2_name")), // old
    _c3_old(getMaterialPropertyOld<Real>("c3_name")), // old
    // _ci_old(getMaterialPropertyOld<Real>(getParam<std::vector<MaterialPropertyName>>))("ci_names")),
    // _ci_old(getMaterialPropertyOld<std::vector<Real>>("test")),

    _ci_IC(getParam<std::vector<Real>>("ci_IC")),

    _dcidc_names(getParam<std::vector<MaterialPropertyName>>("dcidc_names")),
    _prop_dcidc(_num_eta),

    _dcidetaj_names(getParam<std::vector<MaterialPropertyName>>("dcidetaj_names")),
    _prop_dcidetaj(_num_eta),

    _f1(getMaterial("F1_material")),
    _f2(getMaterial("F2_material")),
    _f3(getMaterial("F3_material")),

    _c1_name("c1"),
    _c2_name("c2"),
    _c3_name("c3"),
    _first_df1(getMaterialPropertyDerivative<Real>("F1_name", _c1_name)),
    _first_df2(getMaterialPropertyDerivative<Real>("F2_name", _c2_name)),
    _first_df3(getMaterialPropertyDerivative<Real>("F3_name", _c3_name)),
    _second_df1(getMaterialPropertyDerivative<Real>("F1_name", _c1_name, _c1_name)),
    _second_df2(getMaterialPropertyDerivative<Real>("F2_name", _c2_name, _c2_name)),
    _second_df3(getMaterialPropertyDerivative<Real>("F3_name", _c3_name, _c3_name)),

    _iter(declareProperty<Real>("nested_iterations")),
    _abs_tol(getParam<Real>("absolute_tolerance")),
    _rel_tol(getParam<Real>("relative_tolerance")),
    _nested_solve(NestedSolve(parameters))

{
  // error check parameters
  if (_hj_names.size() != 0 && _hj_names.size() != _num_eta)
    paramError("hj_names",
               "Specify either as many entries are eta values or none at all for auto-naming the "
               "hj material properties.");

  if (_ci_names.size() != 0 && _ci_names.size() != _num_eta)
    paramError("ci_names",
               "Specify either as many entries are eta values or none at all for auto-naming the "
               "ci material properties.");

  // declare ci material properties
  for (unsigned int i = 0; i < _num_eta; ++i)
    _ci_prop[i] = &declareProperty<Real>(_ci_names[i]);

  // declare dcidc material properties
  for (unsigned int i = 0; i < _num_eta; ++i)
    _prop_dcidc[i] = &declareProperty<Real>(_dcidc_names[i]);

  for (unsigned int m = 0; m < _num_eta; ++m)
  {
    // _ci_old_prop[m] = &declareProperty<Real>(_ci_old[i]);

    // declare h and dh material properties
    _prop_hj[m] = &getMaterialPropertyByName<Real>(_hj_names[m]);

    _prop_dcidetaj[m].resize(_num_eta);

    // _prop_dhjdetai[m][n] is the derivative of h_j w.r.t. eta_i
    _prop_dhjdetai[m].resize(_num_eta);

    for (unsigned int n = 0; n < _num_eta; ++n)
      _prop_dhjdetai[m][n] = &getMaterialPropertyDerivative<Real>(_hj_names[m], _eta_names[n]);
  }

  // Get dcidetaj indexes by converting the vector of _dcidetaj_names to the matrix of
  // _prop_dcidetaj, so that _prop_dcidetaj[m][n] is dci[m]/detaj[n]
  for (unsigned int i = 0; i < _num_eta * _num_eta; ++i)
  {
    if (i >= 0 && i < _num_eta)
    {
      _prop_dcidetaj[i][0] = &declareProperty<Real>(_dcidetaj_names[i]);
      continue;
    }
    if (i >= _num_eta && i < 2 * _num_eta)
    {
      _prop_dcidetaj[i - _num_eta][1] = &declareProperty<Real>(_dcidetaj_names[i]);
      continue;
    }
    if (i >= 2 * _num_eta && i < _num_eta * _num_eta)
    {
      _prop_dcidetaj[i - 2 * _num_eta][2] = &declareProperty<Real>(_dcidetaj_names[i]);
      continue;
    }
  }
}

void
SubConcentration::initQpStatefulProperties()
{
  (*_ci_prop[0])[_qp] = _ci_IC[0];
  (*_ci_prop[1])[_qp] = _ci_IC[1];
  (*_ci_prop[2])[_qp] = _ci_IC[2];
}

// // This function is also defined in NestedSolve.C (protected)
// // And the inverse() method does not give the inverse of A, strange
// void
// linear(RankTwoTensor A, RealVectorValue & x, RealVectorValue b)
// {
//   x = A.inverse() * b;
// }

void
SubConcentration::computeQpProperties()
{
  NestedSolve::Value<> solution(3); // dynamicly sized vector class from the Eigen library
  // solution << _ci_IC[0], _ci_IC[1], _ci_IC[2];
  solution << _c1_old[_qp], _c2_old[_qp], _c3_old[_qp];
  // solution << (*_ci_old_prop[0])[_qp], (*_ci_old_prop[1])[_qp], (*_ci_old_prop[2])[_qp];
  // solution << _ci_old[_qp][0], _ci_old[_qp][1], _ci_old[_qp][2];

  _nested_solve.setAbsoluteTolerance(_abs_tol);
  _nested_solve.setRelativeTolerance(_rel_tol);

  auto compute = [&](const NestedSolve::Value<> & guess,
                     NestedSolve::Value<> & residual,
                     NestedSolve::Jacobian<> & jacobian) {
    (*_ci_prop[0])[_qp] = guess(0);
    (*_ci_prop[1])[_qp] = guess(1);
    (*_ci_prop[2])[_qp] = guess(2);
    _f1.computePropertiesAtQp(_qp);
    _f2.computePropertiesAtQp(_qp);
    _f3.computePropertiesAtQp(_qp);

    residual(0) = _first_df1[_qp] - _first_df2[_qp];
    residual(1) = _first_df2[_qp] - _first_df3[_qp];
    residual(2) = (*_prop_hj[0])[_qp] * guess(0) + (*_prop_hj[1])[_qp] * guess(1) +
                  (*_prop_hj[2])[_qp] * guess(2) - _c[_qp];

    jacobian(0, 0) = _second_df1[_qp];
    jacobian(0, 1) = -_second_df2[_qp];
    jacobian(0, 2) = 0;
    jacobian(1, 0) = 0;
    jacobian(1, 1) = _second_df2[_qp];
    jacobian(1, 2) = -_second_df3[_qp];
    jacobian(2, 0) = (*_prop_hj[0])[_qp];
    jacobian(2, 1) = (*_prop_hj[1])[_qp];
    jacobian(2, 2) = (*_prop_hj[2])[_qp];
  };

  _nested_solve.nonlinear(solution, compute);
  _iter[_qp] = _nested_solve.getIterations();

  if (_nested_solve.getState() == NestedSolve::State::NOT_CONVERGED)
  {
    std::cout << "Newton iteration did not converge." << std::endl;
  }

  (*_ci_prop[0])[_qp] = solution[0];
  (*_ci_prop[1])[_qp] = solution[1];
  (*_ci_prop[2])[_qp] = solution[2];

  // update f after solving for ci
  _f1.computePropertiesAtQp(_qp);
  _f2.computePropertiesAtQp(_qp);
  _f3.computePropertiesAtQp(_qp);

  ////////////////////////////////////////////////////////////////////////////////////////// compute dc1dc, dc2dc, and dc3dc
  // The matrix A used to compute dcidc and dcidetai are the same as the jacobian matrix

  RealVectorValue x_dcidc;
  RealVectorValue b_dcidc{0, 0, 1};
  std::vector<std::vector<Real>> A(3);
  for (auto & row : A)
    row.resize(3);

  A[0][0] = _second_df1[_qp];
  A[0][1] = -_second_df2[_qp];
  A[0][2] = 0;
  A[1][0] = 0;
  A[1][1] = _second_df2[_qp];
  A[1][2] = -_second_df3[_qp];
  A[2][0] = (*_prop_hj[0])[_qp];
  A[2][1] = (*_prop_hj[1])[_qp];
  A[2][2] = (*_prop_hj[2])[_qp];

  MatrixTools::inverse(A, A);

  for (unsigned int i = 0; i < _num_eta; ++i)
  {
    x_dcidc(i) = A[i][0] * b_dcidc(0) + A[i][1] * b_dcidc(1) + A[i][2] * b_dcidc(2);
  }

  (*_prop_dcidc[0])[_qp] = x_dcidc(0);
  (*_prop_dcidc[1])[_qp] = x_dcidc(1);
  (*_prop_dcidc[2])[_qp] = x_dcidc(2);

  //////////////////////////////////////////////////////////////////////////////////////////// compute dc1deta1, dc2deta1, and dc3deta1
  RealVectorValue x_dcideta1;
  RealVectorValue b_dcideta1{0,
                             0,
                             -(*_prop_dhjdetai[0][0])[_qp] * (*_ci_prop[0])[_qp] -
                                 (*_prop_dhjdetai[1][0])[_qp] * (*_ci_prop[1])[_qp] -
                                 (*_prop_dhjdetai[2][0])[_qp] * (*_ci_prop[2])[_qp]};

  for (unsigned int i = 0; i < _num_eta; ++i)
  {
    x_dcideta1(i) = A[i][0] * b_dcideta1(0) + A[i][1] * b_dcideta1(1) + A[i][2] * b_dcideta1(2);
  }

  (*_prop_dcidetaj[0][0])[_qp] = x_dcideta1(0);
  (*_prop_dcidetaj[1][0])[_qp] = x_dcideta1(1);
  (*_prop_dcidetaj[2][0])[_qp] = x_dcideta1(2);

  //////////////////////////////////////////////////////////////////////////////////////////// compute dc1deta2, dc2deta2, and dc3deta2
  RealVectorValue x_dcideta2;
  RealVectorValue b_dcideta2{0,
                             0,
                             -(*_prop_dhjdetai[0][1])[_qp] * (*_ci_prop[0])[_qp] -
                                 (*_prop_dhjdetai[1][1])[_qp] * (*_ci_prop[1])[_qp] -
                                 (*_prop_dhjdetai[2][1])[_qp] * (*_ci_prop[2])[_qp]};

  for (unsigned int i = 0; i < _num_eta; ++i)
  {
    x_dcideta2(i) = A[i][0] * b_dcideta2(0) + A[i][1] * b_dcideta2(1) + A[i][2] * b_dcideta2(2);
  }

  (*_prop_dcidetaj[0][1])[_qp] = x_dcideta2(0);
  (*_prop_dcidetaj[1][1])[_qp] = x_dcideta2(1);
  (*_prop_dcidetaj[2][1])[_qp] = x_dcideta2(2);

  //////////////////////////////////////////////////////////////////////////////////////////// compute dc1deta3, dc2deta3, and dc3deta3
  RealVectorValue x_dcideta3;
  RealVectorValue b_dcideta3{0,
                             0,
                             -(*_prop_dhjdetai[0][2])[_qp] * (*_ci_prop[0])[_qp] -
                                 (*_prop_dhjdetai[1][2])[_qp] * (*_ci_prop[1])[_qp] -
                                 (*_prop_dhjdetai[2][2])[_qp] * (*_ci_prop[2])[_qp]};

  for (unsigned int i = 0; i < _num_eta; ++i)
  {
    x_dcideta3(i) = A[i][0] * b_dcideta3(0) + A[i][1] * b_dcideta3(1) + A[i][2] * b_dcideta3(2);
  }

  (*_prop_dcidetaj[0][2])[_qp] = x_dcideta3(0);
  (*_prop_dcidetaj[1][2])[_qp] = x_dcideta3(1);
  (*_prop_dcidetaj[2][2])[_qp] = x_dcideta3(2);
}
