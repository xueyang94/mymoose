//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "MaterialAuxBase.h"

/**
 * A base class for the various Material related AuxKernal objects
 */
template <typename T = Real>
class MaterialStdVectorAuxBase : public MaterialAuxBase<std::vector<T>>
{
public:
  static InputParameters validParams();

  MaterialStdVectorAuxBase(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// index of the vecor element
  unsigned int _index;

  // Explicitly declare the origin of the following inherited members
  // https://isocpp.org/wiki/faq/templates#nondependent-name-lookup-members
  using MaterialAuxBase<std::vector<T>>::_qp;
  using MaterialAuxBase<std::vector<T>>::_prop;
};

template <typename T>
InputParameters
MaterialStdVectorAuxBase<T>::validParams()
{
  InputParameters params = MaterialAuxBase<T>::validParams();
  params.addParam<unsigned int>("index", 0, "The index to consider for this kernel");
  return params;
}

template <typename T>
MaterialStdVectorAuxBase<T>::MaterialStdVectorAuxBase(const InputParameters & parameters)
  : MaterialAuxBase<std::vector<T>>(parameters),
    _index(this->template getParam<unsigned int>("index"))
{
}

template <typename T>
Real
MaterialStdVectorAuxBase<T>::computeValue()
{
  mooseAssert(_prop[_qp].size() > _index,
              "MaterialStdVectorRealGradientAux: You chose to extract component "
                  << _index << " but your Material property only has size " << _prop[_qp].size());
  return MaterialAuxBase<std::vector<T>>::computeValue();
}
