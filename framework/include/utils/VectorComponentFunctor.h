//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseFunctor.h"

/**
 * This is essentially a forwarding functor that forwards the spatial and temporal evaluation
 * arguments to the parent vector functor and then returns the result indexed at a given component.
 */
template <typename T>
class VectorComponentFunctor : public Moose::FunctorBase<T>
{
public:
  using typename Moose::FunctorBase<T>::ValueType;
  using typename Moose::FunctorBase<T>::GradientType;
  using typename Moose::FunctorBase<T>::DotType;
  using VectorArg = typename libMesh::TensorTools::IncrementRank<T>::type;
  using VectorFunctor = Moose::FunctorBase<VectorArg>;

  VectorComponentFunctor(const VectorFunctor & vector, const unsigned int component)
    : Moose::FunctorBase<T>(vector.functorName() + "_" + std::to_string(component)),
      _vector(vector),
      _component(component)
  {
  }

  std::pair<bool, const Elem *> isExtrapolatedBoundaryFace(const FaceInfo & fi) const override
  {
    return _vector.isExtrapolatedBoundaryFace(fi);
  }

private:
  /// The parent vector functor
  const VectorFunctor & _vector;

  /// The component at which we'll index the parent vector functor evaluation result
  const unsigned int _component;

  ValueType evaluate(const Moose::ElemArg & elem, const unsigned int state) const override final
  {
    return _vector(elem, state)(_component);
  }

  ValueType evaluate(const Moose::ElemFromFaceArg & elem_from_face,
                     const unsigned int state) const override final
  {
    return _vector(elem_from_face, state)(_component);
  }

  ValueType evaluate(const Moose::FaceArg & face, const unsigned int state) const override final
  {
    return _vector(face, state)(_component);
  }

  ValueType evaluate(const Moose::SingleSidedFaceArg & face,
                     const unsigned int state) const override final
  {
    return _vector(face, state)(_component);
  }

  ValueType evaluate(const Moose::ElemQpArg & elem_qp,
                     const unsigned int state) const override final
  {
    return _vector(elem_qp, state)(_component);
  }

  ValueType evaluate(const Moose::ElemSideQpArg & elem_side_qp,
                     const unsigned int state) const override final
  {
    return _vector(elem_side_qp, state)(_component);
  }

  using Moose::FunctorBase<T>::evaluateGradient;
  GradientType evaluateGradient(const Moose::ElemArg & elem_arg,
                                const unsigned int state) const override final
  {
    return _vector.gradient(elem_arg, state).row(_component);
  }

  GradientType evaluateGradient(const Moose::FaceArg & face,
                                const unsigned int state) const override final
  {
    return _vector.gradient(face, state).row(_component);
  }
};
