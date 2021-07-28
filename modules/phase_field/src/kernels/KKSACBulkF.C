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
  InputParameters params = Kernel::validParams();
  params.addClassDescription("KKS model kernel (part 1 of 2) for the Bulk Allen-Cahn. This "
                             "includes all terms NOT dependent on chemical potential.");
  params.addRequiredParam<MaterialPropertyName>("L_name", "The name of the Allen-Cahn mobility");
  params.addRequiredParam<MaterialPropertyName>("f1_name", "The name of the first bulk energy");
  params.addRequiredParam<MaterialPropertyName>("f2_name", "The name of the second bulk energy");
  params.addRequiredParam<Real>("barrier_height", "Double well height parameter");
  params.addRequiredCoupledVar("w", "Chemical potential");
  return params;
}

KKSACBulkF::KKSACBulkF(const InputParameters & parameters)
  : Kernel(parameters),
    _L(getMaterialProperty<Real>("L_name")),
    _f1(getMaterialProperty<Real>("f1_name")),
    _f2(getMaterialProperty<Real>("f2_name")),
    _m(getParam<Real>("barrier_height")),
    _w_var(coupled("w")),
    _w(coupledValue("w"))
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
         (-n * (120.0 * n * n - 180.0 * n + 60.0) * (_f1[_qp] - _f2[_qp]) +
          _m * (12.0 * (n * n - n) + 2.0)) *
         _phi[_j][_qp] * _test[_i][_qp];
}

Real
KKSACBulkF::computeQpOffDiagJacobian(unsigned int jvar)
{
  // treat w variable explicitly
  if (jvar == _w_var)
    return 0.0;

  // c is the coupled variable
  return 0.0;
}
