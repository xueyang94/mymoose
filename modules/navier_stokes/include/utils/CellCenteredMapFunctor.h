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
#include "GreenGaussGradient.h"
#include "MathFVUtils.h"
#include "libmesh/utility.h"
#include "libmesh/type_tensor.h"
#include "libmesh/compare_types.h"
#include "libmesh/threads.h"

namespace Moose
{
template <typename T, typename T2, typename std::enable_if<ScalarTraits<T>::value, int>::type = 0>
inline TypeVector<typename CompareTypes<T, T2>::supertype>
outer_product(const T & a, const TypeVector<T2> & b)
{
  TypeVector<typename CompareTypes<T, T2>::supertype> ret;
  for (const auto i : make_range(Moose::dim))
    ret(i) = a * b(i);

  return ret;
}

template <typename T, typename T2>
inline TypeTensor<typename CompareTypes<T, T2>::supertype>
outer_product(const TypeVector<T> & a, const TypeVector<T2> & b)
{
  return libMesh::outer_product(a, b);
}
}

/**
 * A functor whose evaluation relies on querying a map where the keys are element ids and the values
 * correspond to the element/cell values. This is a very useful data type for storing the result of
 * (possibly repeated) interpolation and reconstruction
 */
template <typename T, typename Map>
class CellCenteredMapFunctor : public Moose::FunctorBase<T>, public Map
{
public:
  using typename Moose::FunctorBase<T>::ValueType;
  using typename Moose::FunctorBase<T>::GradientType;
  using typename Moose::FunctorBase<T>::DotType;
  using ElemArg = Moose::ElemArg;
  using ElemFromFaceArg = Moose::ElemFromFaceArg;
  using FaceArg = Moose::FaceArg;
  using SingleSidedFaceArg = Moose::SingleSidedFaceArg;
  using ElemQpArg = Moose::ElemQpArg;
  using ElemSideQpArg = Moose::ElemSideQpArg;

  CellCenteredMapFunctor(const MooseMesh & mesh, const std::string & name)
    : Moose::FunctorBase<T>(name), _mesh(mesh)
  {
  }

  CellCenteredMapFunctor(const MooseMesh & mesh,
                         const std::set<SubdomainID> & sub_ids,
                         const std::string & name)
    : Moose::FunctorBase<T>(name),
      _mesh(mesh),
      _sub_ids(sub_ids == mesh.meshSubdomains() ? std::set<SubdomainID>() : sub_ids)
  {
  }

  std::pair<bool, const Elem *> isExtrapolatedBoundaryFace(const FaceInfo & fi) const override;

private:
  /// The mesh that this functor lives on
  const MooseMesh & _mesh;

  /// The subdomain IDs that this functor lives on. If empty, then we consider the functor to live
  /// on all subdomains
  const std::set<SubdomainID> _sub_ids;

  ValueType evaluate(const ElemArg & elem_arg, unsigned int) const override final
  {
    const Elem * const elem = elem_arg.elem;

    try
    {
      return libmesh_map_find(*this, elem->id());
    }
    catch (libMesh::LogicError &)
    {
      if (!_sub_ids.empty() && !_sub_ids.count(elem->subdomain_id()))
        mooseError("Attempted to evaluate CellCenteredMapFunctor '",
                   this->functorName(),
                   "' with an element subdomain id of '",
                   elem->subdomain_id(),
                   "' but that subdomain id is not one of the subdomain ids the functor is "
                   "restricted to.");
      else
        mooseError("Attempted access into CellCenteredMapFunctor '",
                   this->functorName(),
                   "' with a key that does not yet exist in the map. Make sure to fill your "
                   "CellCenteredMapFunctor for all elements you will attempt to access later.");
    }
  }

  ValueType evaluate(const ElemFromFaceArg & elem_from_face, unsigned int) const override
  {
    const auto * elem = elem_from_face.elem;
    if (!elem)
      elem = &elem_from_face.fi->elem();

    const auto elem_value = (*this)(ElemArg({elem, elem_from_face.correct_skewness}));

    // For the non-boundary elements
    if (elem_from_face.elem)
      return elem_value;

    // We must be on a boundary and want to compute a ghost cell value. Let's determine what the
    // boundary value should be
    mooseAssert(!elem_from_face.fi->neighborPtr(),
                "We should be on a boundary face if we got to this logic");
    const auto boundary_face = SingleSidedFaceArg({elem_from_face.fi,
                                                   Moose::FV::LimiterType::CentralDifference,
                                                   true,
                                                   elem_from_face.correct_skewness,
                                                   elem_from_face.sub_id});
    const auto boundary_value = (*this)(boundary_face);
    // Linear interpolation: face_value = (elem_value + neighbor_value) / 2. Note that weights of
    // 1/2 are perfectly appropriate here because we can arbitrarily put our ghost cell centroid
    // anywhere and by convention we locate it such that a line drawn between the ghost cell
    // centroid and the element centroid is perfectly bisected by the face centroid
    return 2. * boundary_value - elem_value;
  }

  ValueType evaluate(const FaceArg & face, unsigned int) const override final
  {
    const auto & fi = *face.fi;
    mooseAssert(face.limiter_type == Moose::FV::LimiterType::CentralDifference,
                "this implementation currently only supports linear interpolations");

    if (!isExtrapolatedBoundaryFace(fi).first)
      return Moose::FV::linearInterpolation(*this, face);

    if (!fi.neighborPtr() || _sub_ids.count(fi.elem().subdomain_id()))
    {
      const auto elem_arg = face.makeElem();
      const auto elem_value = (*this)(elem_arg);
      // Two term expansion
      return elem_value + this->gradient(elem_arg) * (fi.faceCentroid() - fi.elemCentroid());
    }
    else
    {
      const auto neighbor_arg = face.makeNeighbor();
      const auto neighbor_value = (*this)(neighbor_arg);
      // Two term expansion
      return neighbor_value +
             this->gradient(neighbor_arg) * (fi.faceCentroid() - fi.neighborCentroid());
    }
  }

  using Moose::FunctorBase<T>::evaluateGradient;

  GradientType evaluateGradient(const ElemArg & elem_arg, unsigned int) const override final
  {
    return Moose::FV::greenGaussGradient(elem_arg, *this, true, _mesh);
  }

  GradientType evaluateGradient(const FaceArg & face, unsigned int) const override final
  {
    const auto & fi = *face.fi;
    if (!isExtrapolatedBoundaryFace(fi).first)
    {
      const auto elem_arg = face.makeElem();
      const auto elem_gradient = this->gradient(elem_arg);
      const auto neighbor_arg = face.makeNeighbor();
      const auto linear_interp_gradient =
          fi.gC() * elem_gradient + (1 - fi.gC()) * this->gradient(neighbor_arg);
      return linear_interp_gradient +
             Moose::outer_product(((*this)(neighbor_arg) - (*this)(elem_arg)) / fi.dCFMag() -
                                      linear_interp_gradient * fi.eCF(),
                                  fi.eCF());
    }

    // One term expansion
    if (!fi.neighborPtr() || _sub_ids.count(fi.elem().subdomain_id()))
      return this->gradient(face.makeElem());
    else
      return this->gradient(face.makeNeighbor());
  }

  ValueType evaluate(const SingleSidedFaceArg & ssf, unsigned int) const override
  {
    return (*this)(Moose::FV::makeCDFace(*ssf.fi));
  }

  ValueType evaluate(const ElemQpArg &, unsigned int) const override
  {
    mooseError("not implemented");
  }

  ValueType evaluate(const ElemSideQpArg &, unsigned int) const override
  {
    mooseError("not implemented");
  }
};

template <typename T, typename Map>
std::pair<bool, const Elem *>
CellCenteredMapFunctor<T, Map>::isExtrapolatedBoundaryFace(const FaceInfo & fi) const
{
  if (!fi.neighborPtr())
    return std::make_pair(true, &fi.elem());

  const bool defined_on_elem = _sub_ids.empty() || _sub_ids.count(fi.elem().subdomain_id());
  const bool defined_on_neighbor = _sub_ids.empty() || _sub_ids.count(fi.neighbor().subdomain_id());
  const bool extrapolated = (defined_on_elem + defined_on_neighbor) == 1;

  mooseAssert(defined_on_elem || defined_on_neighbor,
              "This shouldn't be called if we aren't defined on either side.");
  const Elem * const ret_elem = defined_on_elem ? &fi.elem() : fi.neighborPtr();
  return std::make_pair(extrapolated, ret_elem);
}
