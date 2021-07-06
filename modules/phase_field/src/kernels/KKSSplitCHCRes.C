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
  params.addRequiredParam<MaterialPropertyName>("A2_name", "dFadca");
  // params.addCoupledVar("args_a", "Vector of additional arguments to A2_name");
  params.addCoupledVar("args", "args");
  params.addRequiredCoupledVar("w",
                               "Chemical potential non-linear helper variable for the split solve");
  return params;
}

KKSSplitCHCRes::KKSSplitCHCRes(const InputParameters & parameters)
  : DerivativeMaterialInterface<JvarMapKernelInterface<SplitCHBase>>(parameters),
    _A2(getMaterialProperty<Real>("A2_name")),
    // _dA2dc(getMaterialPropertyDerivative<Real>("A2_name", _var.name())), // d(A2)/dc
    // _dA2darg(_n_args),
    _w_var(coupled("w")),
    _w(coupledValue("w"))
{
}
// {
//   // get the second derivative material property
//   for (unsigned int i = 0; i < _n_args; ++i)
//     _dA2darg[i] = &getMaterialPropertyDerivative<Real>("A2_name", i);
// }

// void
// KKSSplitCHCRes::initialSetup()
// {
//   validateNonlinearCoupling<Real>("A2_name");
//   validateDerivativeMaterialPropertyBase<Real>("A2_name");
// }

// Real
// KKSSplitCHCRes::computeQpResidual()
// {
//   Real residual = SplitCHBase::computeQpResidual();
//   residual += -_w[_qp] * _test[_i][_qp];
//
//   return residual;
// }
//
// Real
// KKSSplitCHCRes::computeDFDC(PFFunctionType type)
// {
//   switch (type)
//   {
//     case Residual:
//       return _A2[_qp];
//
//     case Jacobian:
//       return _phi[_j][_qp] * _dA2dc[_qp];
//   }
//
//   mooseError("Invalid type passed in");
// }

Real
KKSSplitCHCRes::computeQpResidual()
{
  return (_A2[_qp] - _w[_qp]) * _test[_i][_qp];
}

// Real
// KKSSplitCHCRes::computeQpJacobian()
// {
//   return _phi[_j][_qp] * _dA2dc[_qp] * _test[_i][_qp];
// }

Real
KKSSplitCHCRes::computeQpJacobian()
{
  return 0.0;
}

Real
KKSSplitCHCRes::computeQpOffDiagJacobian(unsigned int jvar)
{
  // treat w variable explicitly
  // if (jvar == _w_var)
  return -_phi[_j][_qp] * _test[_i][_qp];

  // get the coupled variable jvar is referring to
  // const unsigned int cvar = mapJvarToCvar(jvar);
  // return _phi[_j][_qp] * _test[_i][_qp] * (*_dA2darg[cvar])[_qp];
}
