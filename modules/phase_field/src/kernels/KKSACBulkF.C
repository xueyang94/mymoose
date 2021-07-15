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
  params.addRequiredParam<Real>("barrier_height", "Double well height parameter");
  params.addRequiredCoupledVar("eta_name", "The name of the order parameter");
  params.addRequiredCoupledVar("global_c_name", "The name of the global concentration");
  params.addRequiredCoupledVar("w", "Chemical potential");
  return params;
}

KKSACBulkF::KKSACBulkF(const InputParameters & parameters)
  : Kernel(parameters),
    _L(getMaterialProperty<Real>("L_name")),
    _m(getParam<Real>("barrier_height")),
    _eta(coupledValue("eta_name")),
    _c(coupledValue("global_c_name")),
    _w_var(coupled("w")),
    _w(coupledValue("w"))
{
}

Real
KKSACBulkF::computeQpResidual()
{
  return _L[_qp] *
         (_m * 2.0 * _eta[_qp] * (_eta[_qp] - 1.0) * (2.0 * _eta[_qp] - 1.0) -
          30.0 * _eta[_qp] * _eta[_qp] * (_eta[_qp] * _eta[_qp] - 2.0 * _eta[_qp] + 1.0) *
              (80 * _c[_qp] -
               32 * _eta[_qp] * _eta[_qp] * _eta[_qp] *
                   (6.0 * _eta[_qp] * _eta[_qp] - 15.0 * _eta[_qp] + 10.0) +
               100 * Utility::pow<2>(_c[_qp] -
                                     0.4 * _eta[_qp] * _eta[_qp] * _eta[_qp] *
                                         (6.0 * _eta[_qp] * _eta[_qp] - 15.0 * _eta[_qp] + 10.0)) -
               100 * Utility::pow<2>(_c[_qp] -
                                     0.4 * _eta[_qp] * _eta[_qp] * _eta[_qp] *
                                         (6.0 * _eta[_qp] * _eta[_qp] - 15.0 * _eta[_qp] + 10.0) +
                                     0.4) +
               6)) *
         _test[_i][_qp];
}

Real
KKSACBulkF::computeQpJacobian()
{
  Real n = _eta[_qp];
  n = n > 1 ? 1 : (n < 0 ? 0 : n);

  return _L[_qp] *
         (30.0 * Utility::pow<2>(n) * (Utility::pow<2>(n) - 2.0 * n + 1.0) *
              (200.0 *
                   (_c[_qp] -
                    0.4 * Utility::pow<3>(n) * (6.0 * Utility::pow<2>(n) - 15.0 * n + 10.0)) *
                   (0.4 * Utility::pow<3>(n) * (12.0 * n - 15.0) +
                    1.2 * Utility::pow<2>(n) * (6.0 * Utility::pow<2>(n) - 15.0 * n + 10.0)) +
               32.0 * Utility::pow<3>(n) * (12.0 * n - 15.0) -
               200.0 *
                   (0.4 * Utility::pow<3>(n) * (12.0 * n - 15.0) +
                    1.2 * Utility::pow<2>(n) * (6.0 * Utility::pow<2>(n) - 15.0 * n + 10.0)) *
                   (_c[_qp] -
                    0.4 * Utility::pow<3>(n) * (6.0 * Utility::pow<2>(n) - 15.0 * n + 10.0) + 0.4) +
               96.0 * Utility::pow<2>(n) * (6.0 * Utility::pow<2>(n) - 15.0 * n + 10.0)) +
          4.0 * _m * n * (n - 1.0) -
          60.0 * n * (Utility::pow<2>(n) - 2.0 * n + 1.0) *
              (80.0 * _c[_qp] +
               100.0 * Utility::pow<2>(_c[_qp] - 0.4 * Utility::pow<3>(n) *
                                                     (6.0 * Utility::pow<2>(n) - 15.0 * n + 10.0)) -
               100.0 * Utility::pow<2>(_c[_qp] -
                                       0.4 * Utility::pow<3>(n) *
                                           (6.0 * Utility::pow<2>(n) - 15.0 * n + 10.0) +
                                       0.4) -
               32.0 * Utility::pow<3>(n) * (6.0 * Utility::pow<2>(n) - 15.0 * n + 10.0) + 6.0) +
          2.0 * _m * (n - 1.0) * (2.0 * n - 1.0) + 2.0 * _m * n * (2.0 * n - 1.0) -
          30.0 * Utility::pow<2>(n) * (2.0 * n - 2.0) *
              (80.0 * _c[_qp] +
               100.0 * Utility::pow<2>(_c[_qp] - 0.4 * Utility::pow<3>(n) *
                                                     (6.0 * Utility::pow<2>(n) - 15.0 * n + 10.0)) -
               100.0 * Utility::pow<2>(_c[_qp] -
                                       0.4 * Utility::pow<3>(n) *
                                           (6.0 * Utility::pow<2>(n) - 15.0 * n + 10.0) +
                                       0.4) -
               32.0 * Utility::pow<3>(n) * (6.0 * Utility::pow<2>(n) - 15.0 * n + 10.0) + 6.0)) *
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
