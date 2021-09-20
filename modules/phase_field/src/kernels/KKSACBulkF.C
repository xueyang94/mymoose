//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KKSACBulkF.h"

registerMooseObject("PhaseFieldApp", KKSACBulkF);

InputParameters
KKSACBulkF::validParams()
{
  // InputParameters params = Kernel::validParams();
  InputParameters params = DerivativeMaterialInterface<Kernel>::validParams();
  params.addClassDescription("KKS model kernel (part 1 of 2) for the Bulk Allen-Cahn. This "
                             "includes all terms NOT dependent on chemical potential.");
  params.addRequiredParam<MaterialPropertyName>("dc1dc_name", "The name of dc1/dc");
  params.addRequiredParam<MaterialPropertyName>("dc2dc_name", "The name of dc2/dc");
  params.addRequiredParam<MaterialPropertyName>("dc1deta_name", "The name of dc1/deta");
  params.addRequiredParam<MaterialPropertyName>("dc2deta_name", "The name of dc2/deta");
  params.addRequiredParam<MaterialPropertyName>("F1_name",
                                                "The name of the bulk energy of phase 1");
  params.addRequiredParam<MaterialPropertyName>("F2_name",
                                                "The name of the bulk energy of phase 2");
  params.addRequiredParam<MaterialPropertyName>("L_name", "The name of the Allen-Cahn mobility");
  params.addRequiredParam<Real>("barrier_height", "Double well height parameter");
  params.addRequiredCoupledVar("w",
                               "Chemical potential non-linear helper variable for the split solve");
  return params;
}

KKSACBulkF::KKSACBulkF(const InputParameters & parameters)
  // : Kernel(parameters),
  : DerivativeMaterialInterface<Kernel>(parameters),
    _dc1dc(getMaterialProperty<Real>("dc1dc_name")),
    _dc2dc(getMaterialProperty<Real>("dc2dc_name")),
    _dc1deta(getMaterialProperty<Real>("dc1deta_name")),
    _dc2deta(getMaterialProperty<Real>("dc2deta_name")),
    _f1(getMaterialProperty<Real>("F1_name")),
    _f2(getMaterialProperty<Real>("F2_name")),
    _c1_name("c1"),
    _c2_name("c2"),
    _first_df1(getMaterialPropertyDerivative<Real>("F1_name", _c1_name)),
    _first_df2(getMaterialPropertyDerivative<Real>("F2_name", _c2_name)),
    _L(getMaterialProperty<Real>("L_name")),
    _m(getParam<Real>("barrier_height")),
    _w_var(coupled("w"))
{
}

Real
KKSACBulkF::computeQpResidual()
{
  Real n = _u[_qp];

  return _L[_qp] *
         (-30.0 * n * n * (n * n - 2.0 * n + 1.0) * (_f1[_qp] - _f2[_qp]) +
          _m * 2.0 * n * (n - 1.0) * (2.0 * n - 1.0)) *
         _test[_i][_qp];
}

Real
KKSACBulkF::computeQpJacobian()
{
  Real n = _u[_qp];

  return _L[_qp] *
         (-(n * (120.0 * n * n - 180.0 * n + 60.0) * (_f1[_qp] - _f2[_qp]) +
            30.0 * n * n * (n * n - 2.0 * n + 1.0) *
                (_first_df1[_qp] * _dc1deta[_qp] - _first_df2[_qp] * _dc2deta[_qp])) +
          _m * (12.0 * (n * n - n) + 2.0)) *
         _phi[_j][_qp] * _test[_i][_qp];

  // return _L[_qp] *
  //        (-(n * (120.0 * n * n - 180.0 * n + 60.0) * (_f1[_qp] - _f2[_qp]) +
  //           30.0 * n * n * (n * n - 2.0 * n + 1.0) *
  //               (_first_df1[_qp] * _dc1deta[_qp] - _first_df2[_qp] * _dc2deta[_qp])) +
  //         _m * (12.0 * (n * n - n) + 2.0)) *
  //        _test[_i][_qp];
}

Real
KKSACBulkF::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real n = _u[_qp];

  // treat w variable explicitly
  if (jvar == _w_var)
    return 0.0;

  // c is the coupled variable
  // return _L[_qp] *
  //        (-30.0 * n * n * (n * n - 2.0 * n + 1.0) *
  //         (_first_df1[_qp] * _dc1dc[_qp] - _first_df2[_qp] * _dc2dc[_qp])) *
  //        _phi[_j][_qp] * _test[_i][_qp];
  return _L[_qp] *
         (-30.0 * n * n * (n * n - 2.0 * n + 1.0) *
          (_first_df1[_qp] * _dc1dc[_qp] - _first_df2[_qp] * _dc2dc[_qp])) *
         _test[_i][_qp];
}
