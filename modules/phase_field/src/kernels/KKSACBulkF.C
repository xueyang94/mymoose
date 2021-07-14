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
  params.addRequiredParam<Real>("barrier_height", "Double well height parameter");
  params.addRequiredCoupledVar("eta_name", "The name of the order parameter");
  params.addRequiredParam<MaterialPropertyName>("c1_name", "The name of c1");
  params.addRequiredParam<MaterialPropertyName>("c2_name", "The name of c2");
  params.addRequiredParam<MaterialPropertyName>("dc1dc_name", "The name of dc1/dc");
  params.addRequiredParam<MaterialPropertyName>("dc1deta_name", "The name of dc1/deta");
  params.addRequiredParam<MaterialPropertyName>("dc2dc_name", "The name of dc2/dc");
  params.addRequiredParam<MaterialPropertyName>("dc2deta_name", "The name of dc2/deta");
  params.addRequiredParam<MaterialPropertyName>("L_name", "The name of the Allen-Cahn mobility");
  params.addRequiredParam<MaterialPropertyName>("f1_name",
                                                "The name of the bulk energy of phase 1");
  params.addRequiredParam<MaterialPropertyName>("f2_name",
                                                "The name of the bulk energy of phase 2");
  params.addRequiredCoupledVar("w",
                               "Chemical potential non-linear helper variable for the split solve");
  return params;
}

KKSACBulkF::KKSACBulkF(const InputParameters & parameters)
  : Kernel(parameters),
    _m(getParam<Real>("barrier_height")),
    _eta(coupledValue("eta_name")),
    _c1(getMaterialProperty<Real>("c1_name")),
    _c2(getMaterialProperty<Real>("c2_name")),
    _dc1dc(getMaterialProperty<Real>("dc1dc_name")),
    _dc1deta(getMaterialProperty<Real>("dc1deta_name")),
    _dc2dc(getMaterialProperty<Real>("dc2dc_name")),
    _dc2deta(getMaterialProperty<Real>("dc2deta_name")),
    _L(getMaterialProperty<Real>("L_name")),
    _f1(getMaterialProperty<Real>("f1_name")),
    _f2(getMaterialProperty<Real>("f2_name")),
    _w_var(coupled("w"))
// _w(coupledValue("w"))
{
}

Real
KKSACBulkF::computeQpResidual()
{
  // std::cout << "F" << std::endl;
  return _L[_qp] *
         (-30.0 * _eta[_qp] * _eta[_qp] * (_eta[_qp] * _eta[_qp] - 2.0 * _eta[_qp] + 1.0) *
              (_f1[_qp] - _f2[_qp]) +
          _m * 2.0 * _eta[_qp] * (_eta[_qp] - 1.0) * (2.0 * _eta[_qp] - 1.0)) *
         _test[_i][_qp];
}

Real
KKSACBulkF::computeQpJacobian()
{
  // std::cout << "eta is " << _eta[_qp] << std::endl;
  // std::cout << "c1 is " << _c1[_qp] << std::endl;
  // std::cout << "c2 is " << _c2[_qp] << std::endl;
  // std::cout << "dc1dc is " << _dc1dc[_qp] << std::endl;
  // std::cout << "dc2dc is " << _dc2dc[_qp] << std::endl;
  // std::cout << "dc1deta is " << _dc1deta[_qp] << std::endl;
  // std::cout << "dc2deta is " << _dc2deta[_qp] << std::endl;

  return _L[_qp] *
         (-(_eta[_qp] * (120.0 * _eta[_qp] * _eta[_qp] - 180.0 * _eta[_qp] + 60.0) *
                (_f1[_qp] - _f2[_qp]) +
            30.0 * _eta[_qp] * _eta[_qp] * (_eta[_qp] * _eta[_qp] - 2.0 * _eta[_qp] + 1.0) *
                (_dc1deta[_qp] * (400 * log(_c1[_qp]) - 400 * log(1 - _c1[_qp]) - 28) -
                 _dc2deta[_qp] * (400 * log(_c2[_qp]) - 400 * log(1 - _c2[_qp]) + 39))) +
          _m * (12.0 * (_eta[_qp] * _eta[_qp] - _eta[_qp]) + 2.0)) *
         _phi[_j][_qp] * _test[_i][_qp];
}

Real
KKSACBulkF::computeQpOffDiagJacobian(unsigned int jvar)
{
  // treat w variable explicitly
  if (jvar == _w_var)
    return 0.0;

  // Real df1dc1;
  // Real df2dc2;
  //
  // df1dc1 = 400 * log(_c1[_qp]) - 400 * log(1 - _c1[_qp]) - 28;
  // df2dc2 = 400 * log(_c2[_qp]) - 400 * log(1 - _c2[_qp]) + 39;
  //
  // std::cout << "The df1 is " << df1dc1 << std::endl;
  // std::cout << "The df2 is " << df2dc2 << std::endl;

  return _L[_qp] * (-1) * 30.0 * _eta[_qp] * _eta[_qp] *
         (_eta[_qp] * _eta[_qp] - 2.0 * _eta[_qp] + 1.0) *
         (_dc1dc[_qp] * (400 * log(_c1[_qp]) - 400 * log(1 - _c1[_qp]) - 28) -
          _dc2dc[_qp] * (400 * log(_c2[_qp]) - 400 * log(1 - _c2[_qp]) + 39)) *
         _phi[_j][_qp] * _test[_i][_qp];

  // // // c is the coupled variable
  // return _L[_qp] * (-1) * 30.0 * _eta[_qp] * _eta[_qp] *
  //        (_eta[_qp] * _eta[_qp] - 2.0 * _eta[_qp] + 1.0) *
  //        (400 * log(_c1[_qp]) - 400 * log(1 - _c1[_qp]) - 28) * (_dc1dc[_qp] - _dc2dc[_qp]) *
  //        _phi[_j][_qp] * _test[_i][_qp];
}
