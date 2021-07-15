//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KKSACBulkC.h"

registerMooseObject("PhaseFieldApp", KKSACBulkC);

InputParameters
KKSACBulkC::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("KKS model kernel (part 2 of 2) for the Bulk Allen-Cahn. This "
                             "includes all terms dependent on chemical potential.");
  params.addRequiredParam<MaterialPropertyName>("L_name", "The name of the Allen-Cahn mobility");
  params.addRequiredParam<MaterialPropertyName>("c1_name",
                                                "The name of the first sub-concentration");
  params.addRequiredCoupledVar("w", "Chemical potential");
  return params;
}

KKSACBulkC::KKSACBulkC(const InputParameters & parameters)
  : Kernel(parameters),
    _L(getMaterialProperty<Real>("L_name")),
    _c1(getMaterialProperty<Real>("c1_name")),
    _w_var(coupled("w")),
    _w(coupledValue("w"))
{
}

Real
KKSACBulkC::computeQpResidual()
{
  Real n = _u[_qp];

  return _L[_qp] * 30.0 * n * n * (n * n - 2.0 * n + 1.0) * 200 * (_c1[_qp] - 0.3) * (-0.4) *
         _test[_i][_qp];
}

Real
KKSACBulkC::computeQpJacobian()
{
  Real n = _u[_qp];

  return _L[_qp] * (-80) *
         (n * (120.0 * n * n - 180.0 * n + 60.0) * (_c1[_qp] - 0.3) -
          0.4 * Utility::pow<2>(30.0 * n * n * (n * n - 2.0 * n + 1.0))) *
         _phi[_j][_qp] * _test[_i][_qp];
}

Real
KKSACBulkC::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _w_var)
    return 0.0;

  Real n = _u[_qp];
  // c is the coupled variable
  return _L[_qp] * (-80) * 30.0 * n * n * (n * n - 2.0 * n + 1.0) * _phi[_j][_qp] * _test[_i][_qp];
}
