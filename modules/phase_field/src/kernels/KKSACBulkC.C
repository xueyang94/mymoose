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
  // InputParameters params = KKSACBulkBase::validParams();
  InputParameters params = Kernel::validParams();
  params.addClassDescription("KKS model kernel (part 2 of 2) for the Bulk Allen-Cahn. This "
                             "includes all terms dependent on chemical potential.");
  // params.addRequiredParam<MaterialPropertyName>("A1_name", "The product of dFadca and ca-cb");
  params.addRequiredCoupledVar("eta_name", "The name of the order parameter");
  params.addRequiredCoupledVar("global_c_name", "The name of the global concentration");
  params.addRequiredCoupledVar("w", "Chemical potential");
  // params.addRequiredParam<MaterialPropertyName>("h_name", "The product of dFadca and ca-cb");
  return params;
}

KKSACBulkC::KKSACBulkC(const InputParameters & parameters)
  // : KKSACBulkBase(parameters),
  : Kernel(parameters),
    // _A1(getMaterialProperty<Real>("A1_name")),
    // _prop_dA1(getMaterialPropertyDerivative<Real>("A1_name", _eta_name)),
    // _prop_dA1darg(_n_args),
    // _eta_name(getVar("eta_name", 0)->name()),
    _eta(coupledValue("eta_name")),
    _c(coupledValue("global_c_name")),
    _w_var(coupled("w")),
    _w(coupledValue("w"))
// _dh(getMaterialPropertyDerivative<Real>("h_name", _eta_name))
{
}
// {
//   // // get second partial derivatives wrt ca and other coupled variable
//   for (unsigned int i = 0; i < _n_args; ++i)
//     _prop_dA1darg[i] = &getMaterialPropertyDerivative<Real>("A1_name", i);
// }

// Real
// KKSACBulkC::computeDFDOP(PFFunctionType type)
// {
//   switch (type)
//   {
//     case Residual:
//       return 30.0 * _eta[_qp] * _eta[_qp] * (_eta[_qp] * _eta[_qp] - 2.0 * _eta[_qp] + 1.0) * 80
//       *
//              (0.3 +
//               0.4 * (_eta[_qp] * _eta[_qp] * _eta[_qp] *
//                      (6.0 * _eta[_qp] * _eta[_qp] - 15.0 * _eta[_qp] + 10.0)) -
//               _c[_qp]);
//
//     case Jacobian:
//       return 80 *
//              (0.3 * _eta[_qp] * (120.0 * _eta[_qp] * _eta[_qp] - 180.0 * _eta[_qp] + 60.0) +
//               0.4 * _eta[_qp] * (120.0 * _eta[_qp] * _eta[_qp] - 180.0 * _eta[_qp] + 60.0) *
//                   _eta[_qp] * _eta[_qp] * _eta[_qp] *
//                   (6.0 * _eta[_qp] * _eta[_qp] - 15.0 * _eta[_qp] + 10.0) -
//               _c[_qp] * _eta[_qp] * (120.0 * _eta[_qp] * _eta[_qp] - 180.0 * _eta[_qp] + 60.0) +
//               0.4 * 30.0 * _eta[_qp] * _eta[_qp] * (_eta[_qp] * _eta[_qp] - 2.0 * _eta[_qp]
//               + 1.0) *
//                   30.0 * _eta[_qp] * _eta[_qp] * (_eta[_qp] * _eta[_qp] - 2.0 * _eta[_qp] + 1.0))
//                   *
//              _phi[_j][_qp];
//   }
//
//   mooseError("Invalid type passed in");
// }

Real
KKSACBulkC::computeQpResidual()
{

  // Real n = _eta[_qp];
  // n = n > 1 ? 1 : (n < 0 ? 0 : n);
  // Real h;
  // h = n * n * n * (6.0 * n * n - 15.0 * n + 10.0);
  // Real d1h;
  // d1h = 30.0 * n * n * (n * n - 2.0 * n + 1.0);
  // Real d2h;
  // d2h = n * (120.0 * n * n - 180.0 * n + 60.0);
  //
  // return d1h * 80 * (0.3 + 0.4 * h - _c[_qp]) * _test[_i][_qp];

  return 30.0 * _eta[_qp] * _eta[_qp] * (_eta[_qp] * _eta[_qp] - 2.0 * _eta[_qp] + 1.0) * 80 *
         (0.3 +
          0.4 * (_eta[_qp] * _eta[_qp] * _eta[_qp] *
                 (6.0 * _eta[_qp] * _eta[_qp] - 15.0 * _eta[_qp] + 10.0)) -
          _c[_qp]) *
         _test[_i][_qp];
}

Real
KKSACBulkC::computeQpJacobian()
{
  return 80 *
         (0.3 * _eta[_qp] * (120.0 * _eta[_qp] * _eta[_qp] - 180.0 * _eta[_qp] + 60.0) +
          0.4 * _eta[_qp] * (120.0 * _eta[_qp] * _eta[_qp] - 180.0 * _eta[_qp] + 60.0) * _eta[_qp] *
              _eta[_qp] * _eta[_qp] * (6.0 * _eta[_qp] * _eta[_qp] - 15.0 * _eta[_qp] + 10.0) -
          _c[_qp] * _eta[_qp] * (120.0 * _eta[_qp] * _eta[_qp] - 180.0 * _eta[_qp] + 60.0) +
          0.4 * 30.0 * _eta[_qp] * _eta[_qp] * (_eta[_qp] * _eta[_qp] - 2.0 * _eta[_qp] + 1.0) *
              30.0 * _eta[_qp] * _eta[_qp] * (_eta[_qp] * _eta[_qp] - 2.0 * _eta[_qp] + 1.0)) *
         _phi[_j][_qp] * _test[_i][_qp];
}

Real
KKSACBulkC::computeQpOffDiagJacobian(unsigned int jvar)
{
  // first get dependence of mobility _L on other variables using parent class
  // member function
  // Real res = ACBulk<Real>::computeQpOffDiagJacobian(jvar);

  //  for all other vars get the coupled variable jvar is referring to
  // const unsigned int cvar = mapJvarToCvar(jvar);

  // res += _L[_qp] * _prop_dh[_qp] * (*_prop_dA1darg[cvar])[_qp] * _phi[_j][_qp] * _test[_i][_qp];

  // return res;

  if (jvar == _w_var)
    return 0.0;

  // c is the coupled variable
  return -80 * 30.0 * _eta[_qp] * _eta[_qp] * (_eta[_qp] * _eta[_qp] - 2.0 * _eta[_qp] + 1.0) *
         _phi[_j][_qp] * _test[_i][_qp];
}
