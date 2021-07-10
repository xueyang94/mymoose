//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KKSSplitCHCRes.h"

registerMooseObject("PhaseFieldApp", KKSSplitCHCRes);

InputParameters
KKSSplitCHCRes::validParams()
{
  InputParameters params = SplitCHBase::validParams();
  params.addClassDescription(
      "KKS model kernel for the split Bulk Cahn-Hilliard term. This kernel operates on the "
      "physical concentration 'c' as the non-linear variable");
  // params.addRequiredParam<MaterialPropertyName>("A2_name", "dFadca");
  params.addRequiredCoupledVar("eta_name", "The name of the order parameter");
  params.addRequiredCoupledVar("global_c_name", "The name of the global concentration");
  params.addRequiredCoupledVar("w",
                               "Chemical potential non-linear helper variable for the split solve");
  return params;
}

KKSSplitCHCRes::KKSSplitCHCRes(const InputParameters & parameters)
  : DerivativeMaterialInterface<JvarMapKernelInterface<SplitCHBase>>(parameters),
    // _A2(getMaterialProperty<Real>("A2_name")),
    _eta(coupledValue("eta_name")),
    _c(coupledValue("global_c_name")),
    _w_var(coupled("w")),
    _w(coupledValue("w"))

{
}

Real
KKSSplitCHCRes::computeQpResidual()
{
  return (200 * (_c[_qp] -
                 0.4 * _eta[_qp] * _eta[_qp] * _eta[_qp] *
                     (6.0 * _eta[_qp] * _eta[_qp] - 15.0 * _eta[_qp] + 10.0) -
                 0.3) -
          _w[_qp]) *
         _test[_i][_qp];
}

Real
KKSSplitCHCRes::computeQpJacobian()
{
  return 200 * _phi[_j][_qp] * _test[_i][_qp];
}

Real
KKSSplitCHCRes::computeQpOffDiagJacobian(unsigned int jvar)
{
  // treat w variable explicitly
  if (jvar == _w_var)
    return -_phi[_j][_qp] * _test[_i][_qp];

  // eta is the couple variable
  return -80 * (30.0 * _eta[_qp] * _eta[_qp] * (_eta[_qp] * _eta[_qp] - 2.0 * _eta[_qp] + 1.0)) *
         _phi[_j][_qp] * _test[_i][_qp];
}
