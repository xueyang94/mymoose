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
  InputParameters params = KKSACBulkBase::validParams();
  params.addClassDescription("KKS model kernel (part 2 of 2) for the Bulk Allen-Cahn. This "
                             "includes all terms dependent on chemical potential.");
  params.addRequiredParam<MaterialPropertyName>("A1_name", "The product of dFadca and ca-cb");
  return params;
}

KKSACBulkC::KKSACBulkC(const InputParameters & parameters)
  : KKSACBulkBase(parameters),
    _prop_A1(getMaterialProperty<Real>("A1_name")),
    _prop_dA1(getMaterialPropertyDerivative<Real>("A1_name", _eta_name)),
    _prop_dA1darg(_n_args)
{
  // get second partial derivatives wrt ca and other coupled variable
  for (unsigned int i = 0; i < _n_args; ++i)
    _prop_dA1darg[i] = &getMaterialPropertyDerivative<Real>("A1_name", i);
}

Real
KKSACBulkC::computeDFDOP(PFFunctionType type)
{
  switch (type)
  {
    case Residual:
      return _prop_dh[_qp] * _prop_A1[_qp];

    case Jacobian:
      return _phi[_j][_qp] * (_prop_d2h[_qp] * _prop_A1[_qp] + _prop_dh[_qp] * _prop_dA1[_qp]);
      // return _phi[_j][_qp] * (_prop_d2h[_qp] * _prop_A1[_qp]);
  }

  mooseError("Invalid type passed in");
}

Real
KKSACBulkC::computeQpOffDiagJacobian(unsigned int jvar)
{
  // first get dependence of mobility _L on other variables using parent class
  // member function
  Real res = ACBulk<Real>::computeQpOffDiagJacobian(jvar);

  //  for all other vars get the coupled variable jvar is referring to
  const unsigned int cvar = mapJvarToCvar(jvar);

  res += _L[_qp] * _prop_dh[_qp] * (*_prop_dA1darg[cvar])[_qp] * _phi[_j][_qp] * _test[_i][_qp];

  return res;
}
