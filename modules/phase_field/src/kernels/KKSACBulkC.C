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
  params.addRequiredCoupledVar("eta_name", "The name of the order parameter");
  params.addRequiredParam<MaterialPropertyName>("c1_name", "The name of c1");
  params.addRequiredParam<MaterialPropertyName>("c2_name", "The name of c2");
  params.addRequiredParam<MaterialPropertyName>("dc1dc_name", "The name of dc1/dc");
  params.addRequiredParam<MaterialPropertyName>("dc1deta_name", "The name of dc1/deta");
  params.addRequiredParam<MaterialPropertyName>("dc2dc_name", "The name of dc2/dc");
  params.addRequiredParam<MaterialPropertyName>("dc2deta_name", "The name of dc2/deta");
  params.addRequiredCoupledVar("w",
                               "Chemical potential non-linear helper variable for the split solve");
  return params;
}

KKSACBulkC::KKSACBulkC(const InputParameters & parameters)
  : Kernel(parameters),
    _eta(coupledValue("eta_name")),
    _c1(getMaterialProperty<Real>("c1_name")),
    _c2(getMaterialProperty<Real>("c2_name")),
    _dc1dc(getMaterialProperty<Real>("dc1dc_name")),
    _dc1deta(getMaterialProperty<Real>("dc1deta_name")),
    _dc2dc(getMaterialProperty<Real>("dc2dc_name")),
    _dc2deta(getMaterialProperty<Real>("dc2deta_name")),
    _w_var(coupled("w")),
    _w(coupledValue("w"))
{
}

Real
KKSACBulkC::computeQpResidual()
{
  return (30.0 * _eta[_qp] * _eta[_qp] * (_eta[_qp] * _eta[_qp] - 2.0 * _eta[_qp] + 1.0)) *
         (400 * log(_c1[_qp]) - 400 * log(1 - _c1[_qp]) - 28 - _w[_qp]) * (_c1[_qp] - _c2[_qp]) *
         _test[_i][_qp];
}

Real
KKSACBulkC::computeQpJacobian()
{
  return (_eta[_qp] * (120.0 * _eta[_qp] * _eta[_qp] - 180.0 * _eta[_qp] + 60.0) *
              (400 * log(_c1[_qp]) - 400 * log(1 - _c1[_qp]) - 28) * (_c1[_qp] - _c2[_qp]) +
          30.0 * _eta[_qp] * _eta[_qp] * (_eta[_qp] * _eta[_qp] - 2.0 * _eta[_qp] + 1.0) *
              (_c1[_qp] - _c2[_qp]) * _dc1deta[_qp] * (400 / _c1[_qp] + 400 / (1 - _c1[_qp])) +
          (400 * log(_c1[_qp]) - 400 * log(1 - _c1[_qp]) - 28) * (_dc1deta[_qp] - _dc2deta[_qp])) *
         _phi[_j][_qp] * _test[_i][_qp];
}

Real
KKSACBulkC::computeQpOffDiagJacobian(unsigned int jvar)
{
  // treat w variable explicitly
  if (jvar == _w_var)
    return 0.0;

  // c is the coupled variable
  return 30.0 * _eta[_qp] * _eta[_qp] * (_eta[_qp] * _eta[_qp] - 2.0 * _eta[_qp] + 1.0) *
         ((_c1[_qp] - _c2[_qp]) * _dc1dc[_qp] * (400 / _c1[_qp] + 400 / (1 - _c1[_qp])) +
          (400 * log(_c1[_qp]) - 400 * log(1 - _c1[_qp]) - 28) * (_dc1dc[_qp] - _dc2dc[_qp])) *
         _phi[_j][_qp] * _test[_i][_qp];
}
