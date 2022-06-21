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
#include "MathFVUtils.h"
#include "FVUtils.h"
#include "MooseMeshUtils.h"
#include "VectorComponentFunctor.h"
#include "ArrayComponentFunctor.h"
#include "libmesh/elem.h"

namespace Moose
{
namespace FV
{
/**
 * Compute a cell gradient using the method of Green-Gauss
 * @param elem_arg An element argument specifying the current element and whether to perform skew
 * corrections
 * @param functor The functor that will provide information such as cell and face value evaluations
 * necessary to construct the cell gradient
 * @param two_term_boundary_expansion Whether to perform a two-term expansion to compute
 * extrapolated boundary face values. If this is true, then an implicit system has to be solved. If
 * false, then the cell center value will be used as the extrapolated boundary face value
 * @param mesh The mesh on which we are computing the gradient
 * @param face_to_value_cache An optional parameter. If provided, we will add face to
 * face-evaluations computed in this function to the map
 * @return The computed cell gradient
 */
template <typename T, typename Enable = typename std::enable_if<ScalarTraits<T>::value>::type>
VectorValue<T>
greenGaussGradient(const ElemArg & elem_arg,
                   const FunctorBase<T> & functor,
                   const bool two_term_boundary_expansion,
                   const MooseMesh & mesh,
                   std::unordered_map<const FaceInfo *, T> * const face_to_value_cache = nullptr)
{
  mooseAssert(elem_arg.elem, "This should be non-null");
  const auto coord_type = mesh.getCoordSystem(elem_arg.elem->subdomain_id());
  const auto rz_radial_coord = mesh.getAxisymmetricRadialCoord();

  T elem_value = functor(elem_arg);

  // We'll save off the extrapolated boundary faces (ebf) for later assignment to the cache (these
  // are the keys)
  std::vector<const FaceInfo *> ebf_faces;

  try
  {
    VectorValue<T> grad;

    bool volume_set = false;
    Real volume = 0;

    // If we are performing a two term Taylor expansion for extrapolated boundary faces (faces on
    // boundaries that do not have associated Dirichlet conditions), then the element gradient
    // depends on the boundary face value and the boundary face value depends on the element
    // gradient, so we have a system of equations to solve. Here is the system:
    //
    // \nabla \phi_C - \frac{1}{V} \sum_{ebf} \phi_{ebf} \vec{S_f} =
    //   \frac{1}{V} \sum_{of} \phi_{of} \vec{S_f}                       eqn. 1
    //
    // \phi_{ebf} - \vec{d_{Cf}} \cdot \nabla \phi_C = \phi_C            eqn. 2
    //
    // where $C$ refers to the cell centroid, $ebf$ refers to an extrapolated boundary face, $of$
    // refers to "other faces", e.g. non-ebf faces, and $f$ is a general face. $d_{Cf}$ is the
    // vector drawn from the element centroid to the face centroid, and $\vec{S_f}$ is the surface
    // vector, e.g. the face area times the outward facing normal

    // ebf eqns: element gradient coefficients, e.g. eqn. 2, LHS term 2 coefficient
    std::vector<VectorValue<Real>> ebf_grad_coeffs;
    // ebf eqns: rhs b values. These will actually correspond to the elem_value so we can use a
    // pointer and avoid copying. This is the RHS of eqn. 2
    std::vector<const T *> ebf_b;

    // elem grad eqns: ebf coefficients, e.g. eqn. 1, LHS term 2 coefficients
    std::vector<VectorValue<Real>> grad_ebf_coeffs;
    // elem grad eqns: rhs b value, e.g. eqn. 1 RHS
    VectorValue<T> grad_b = 0;

    auto action_functor = [&volume_set,
                           &volume,
                           &elem_value,
                           &elem_arg,
                           &ebf_faces,
                           &ebf_grad_coeffs,
                           &ebf_b,
                           &grad_ebf_coeffs,
                           &grad_b,
                           &functor,
                           two_term_boundary_expansion,
                           coord_type,
                           rz_radial_coord](const Elem & libmesh_dbg_var(functor_elem),
                                            const Elem *,
                                            const FaceInfo * const fi,
                                            const Point & surface_vector,
                                            Real coord,
                                            const bool elem_has_info)
    {
      mooseAssert(fi, "We need a FaceInfo for this action_functor");
      mooseAssert(elem_arg.elem == &functor_elem,
                  "Just a sanity check that the element being passed in is the one we passed out.");

      if (functor.isExtrapolatedBoundaryFace(*fi).first)
      {
        if (two_term_boundary_expansion)
        {
          ebf_faces.push_back(fi);

          // eqn. 2
          ebf_grad_coeffs.push_back(-1. * (elem_has_info
                                               ? (fi->faceCentroid() - fi->elemCentroid())
                                               : (fi->faceCentroid() - fi->neighborCentroid())));
          ebf_b.push_back(&elem_value);

          // eqn. 1
          grad_ebf_coeffs.push_back(-surface_vector);
        }
        else
          // We are doing a one-term expansion for the extrapolated boundary faces, in which case we
          // have no eqn. 2 and we have no second term in the LHS of eqn. 1. Instead we apply the
          // element centroid value as the face value (one-term expansion) in the RHS of eqn. 1
          grad_b += surface_vector * elem_value;
      }
      else
        grad_b += surface_vector * functor(makeCDFace(*fi, elem_arg.correct_skewness));

      if (!volume_set)
      {
        // We use the FaceInfo volumes because those values have been pre-computed and cached.
        // An explicit call to elem->volume() here would incur unnecessary expense
        if (elem_has_info)
        {
          MooseMeshUtils::coordTransformFactor(
              fi->elemCentroid(), coord, coord_type, rz_radial_coord);
          volume = fi->elemVolume() * coord;
        }
        else
        {
          MooseMeshUtils::coordTransformFactor(
              fi->neighborCentroid(), coord, coord_type, rz_radial_coord);
          volume = fi->neighborVolume() * coord;
        }

        volume_set = true;
      }
    };

    Moose::FV::loopOverElemFaceInfo(
        *elem_arg.elem, mesh, action_functor, coord_type, rz_radial_coord);

    mooseAssert(volume_set && volume > 0, "We should have set the volume");
    grad_b /= volume;

    if (coord_type == Moose::CoordinateSystemType::COORD_RZ)
    {
      mooseAssert(rz_radial_coord != libMesh::invalid_uint, "rz_radial_coord must be set");
      grad_b(rz_radial_coord) -= elem_value / elem_arg.elem->vertex_average()(rz_radial_coord);
    }

    mooseAssert(
        coord_type != Moose::CoordinateSystemType::COORD_RSPHERICAL,
        "We have not yet implemented the correct translation from gradient to divergence for "
        "spherical coordinates yet.");

    mooseAssert(
        ebf_faces.size() < UINT_MAX,
        "You've created a mystical element that has more faces than can be held by unsigned "
        "int. I applaud you.");
    const auto num_ebfs = static_cast<unsigned int>(ebf_faces.size());

    // test for simple case
    if (num_ebfs == 0)
      grad = grad_b;
    else
    {
      // We have to solve a system
      const unsigned int sys_dim = Moose::dim + num_ebfs;
      DenseVector<T> x(sys_dim), b(sys_dim);
      DenseMatrix<T> A(sys_dim, sys_dim);

      // Let's make i refer to Moose::dim indices, and j refer to num_ebfs indices

      // eqn. 1
      for (const auto i : make_range(Moose::dim))
      {
        // LHS term 1 coeffs
        A(i, i) = 1;

        // LHS term 2 coeffs
        for (const auto j : make_range(num_ebfs))
          A(i, Moose::dim + j) = grad_ebf_coeffs[j](i) / volume;

        // RHS
        b(i) = grad_b(i);
      }

      // eqn. 2
      for (const auto j : make_range(num_ebfs))
      {
        // LHS term 1 coeffs
        A(Moose::dim + j, Moose::dim + j) = 1;

        // LHS term 2 coeffs
        for (const auto i : make_range(unsigned(Moose::dim)))
          A(Moose::dim + j, i) = ebf_grad_coeffs[j](i);

        // RHS
        b(Moose::dim + j) = *ebf_b[j];
      }

      A.lu_solve(b, x);
      for (const auto i : make_range(Moose::dim))
        grad(i) = x(i);

      // Optionally cache the face value information
      if (face_to_value_cache)
        for (const auto j : make_range(num_ebfs))
          face_to_value_cache->emplace(ebf_faces[j], x(Moose::dim + j));
    }

    return grad;
  }
  catch (libMesh::LogicError &)
  {
    // Retry without two-term
    mooseAssert(two_term_boundary_expansion,
                "I believe we should only get singular systems when two-term boundary expansion is "
                "being used");
    const auto grad = greenGaussGradient(elem_arg, functor, false, mesh, face_to_value_cache);

    // We failed to compute the extrapolated boundary faces with two-term expansion and callers of
    // this method may be relying on those values (e.g. if the caller is
    // getExtrapolatedBoundaryFaceValue) so we populate them here with one-term expansion, e.g. we
    // set the boundary face values to the cell centroid value
    if (face_to_value_cache)
      for (auto * const ebf_face : ebf_faces)
        face_to_value_cache->emplace(ebf_face, elem_value);

    return grad;
  }
}

/**
 * Compute a face gradient from Green-Gauss cell gradients, with orthogonality correction
 * On the boundaries, the boundary element value is used
 * @param face_arg A face argument specifying the current faceand whether to perform skew
 * corrections
 * @param functor The functor that will provide information such as cell and face value evaluations
 * necessary to construct the cell gradient
 * @param two_term_boundary_expansion Whether to perform a two-term expansion to compute
 * extrapolated boundary face values. If this is true, then an implicit system has to be solved. If
 * false, then the cell center value will be used as the extrapolated boundary face value
 * @param mesh The mesh on which we are computing the gradient
 * @param face_to_value_cache An optional parameter. If provided, we will add face to
 * face-evaluations computed in this function to the map
 * @return The computed cell gradient
 */
template <typename T, typename Enable = typename std::enable_if<ScalarTraits<T>::value>::type>
VectorValue<T>
greenGaussGradient(const FaceArg & face_arg,
                   const FunctorBase<T> & functor,
                   const bool two_term_boundary_expansion,
                   const MooseMesh & mesh,
                   std::unordered_map<const FaceInfo *, T> * const face_to_value_cache = nullptr)
{
  mooseAssert(face_arg.fi, "We should have a face info to compute a face gradient");
  const auto & fi = *(face_arg.fi);
  const auto & elem_arg = face_arg.makeElem();
  const auto & neighbor_arg = face_arg.makeNeighbor();
  const auto [is_extrapolated, defined_elem] = functor.isExtrapolatedBoundaryFace(fi);

  if (!is_extrapolated)
  {
    // Compute the gradients in the two cells on both sides of the face
    const auto & grad_elem = greenGaussGradient(
        elem_arg, functor, two_term_boundary_expansion, mesh, face_to_value_cache);
    const auto & grad_neighbor = greenGaussGradient(
        neighbor_arg, functor, two_term_boundary_expansion, mesh, face_to_value_cache);

    VectorValue<T> face_gradient;
    Moose::FV::interpolate(
        Moose::FV::InterpMethod::Average, face_gradient, grad_elem, grad_neighbor, fi, true);

    // Perform orthogonality correction
    const auto & value_elem = functor(elem_arg);
    const auto & value_neighbor = functor(neighbor_arg);

    face_gradient += outer_product(
        fi.eCF(), (value_neighbor - value_elem) / fi.dCFMag() - face_gradient * fi.eCF());
    return face_gradient;
  }

  // One term expansion
  if (&fi.elem() == defined_elem)
    return functor.gradient(elem_arg);
  else
    return functor.gradient(neighbor_arg);
}

template <typename T>
TensorValue<T>
greenGaussGradient(const ElemArg & elem_arg,
                   const Moose::FunctorBase<VectorValue<T>> & functor,
                   const bool two_term_boundary_expansion,
                   const MooseMesh & mesh)
{
  TensorValue<T> ret;
  for (const auto i : make_range(Moose::dim))
  {
    VectorComponentFunctor<T> scalar_functor(functor, i);
    const auto row_gradient =
        greenGaussGradient(elem_arg, scalar_functor, two_term_boundary_expansion, mesh);
    for (const auto j : make_range(unsigned(Moose::dim)))
      ret(i, j) = row_gradient(j);
  }

  return ret;
}

template <typename T>
TensorValue<T>
greenGaussGradient(const FaceArg & face_arg,
                   const Moose::FunctorBase<VectorValue<T>> & functor,
                   const bool two_term_boundary_expansion,
                   const MooseMesh & mesh)
{
  TensorValue<T> ret;
  for (const auto i : make_range(unsigned(Moose::dim)))
  {
    VectorComponentFunctor<T> scalar_functor(functor, i);
    const auto row_gradient =
        greenGaussGradient(face_arg, scalar_functor, two_term_boundary_expansion, mesh);
    for (const auto j : make_range(unsigned(Moose::dim)))
      ret(i, j) = row_gradient(j);
  }

  return ret;
}

template <typename T>
typename Moose::FunctorBase<std::vector<T>>::GradientType
greenGaussGradient(const ElemArg & elem_arg,
                   const Moose::FunctorBase<std::vector<T>> & functor,
                   const bool two_term_boundary_expansion,
                   const MooseMesh & mesh)
{
  // Determine the size of the container
  const auto vals = functor(elem_arg);
  typedef typename Moose::FunctorBase<std::vector<T>>::GradientType GradientType;
  GradientType ret(vals.size());
  for (const auto i : index_range(ret))
  {
    // Note that this can be very inefficient. Within the scalar greenGaussGradient routine we're
    // going to do value type evaluations of the array functor from scalar_functor and we will be
    // discarding all the value type evaluations other than the one corresponding to i
    ArrayComponentFunctor<T, FunctorBase<std::vector<T>>> scalar_functor(functor, i);
    ret[i] = greenGaussGradient(elem_arg, scalar_functor, two_term_boundary_expansion, mesh);
  }

  return ret;
}

template <typename T>
typename Moose::FunctorBase<std::vector<T>>::GradientType
greenGaussGradient(const FaceArg & face_arg,
                   const Moose::FunctorBase<std::vector<T>> & functor,
                   const bool two_term_boundary_expansion,
                   const MooseMesh & mesh)
{
  // Determine the size of the container
  const auto vals = functor(face_arg);
  typedef typename Moose::FunctorBase<std::vector<T>>::GradientType GradientType;
  GradientType ret(vals.size());
  for (const auto i : index_range(ret))
  {
    // Note that this can be very inefficient. Within the scalar greenGaussGradient routine we're
    // going to do value type evaluations of the array functor from scalar_functor and we will be
    // discarding all the value type evaluations other than the one corresponding to i
    ArrayComponentFunctor<T, FunctorBase<std::vector<T>>> scalar_functor(functor, i);
    ret[i] = greenGaussGradient(face_arg, scalar_functor, two_term_boundary_expansion, mesh);
  }

  return ret;
}

template <typename T, std::size_t N>
typename Moose::FunctorBase<std::array<T, N>>::GradientType
greenGaussGradient(const ElemArg & elem_arg,
                   const Moose::FunctorBase<std::array<T, N>> & functor,
                   const bool two_term_boundary_expansion,
                   const MooseMesh & mesh)
{
  typedef typename Moose::FunctorBase<std::array<T, N>>::GradientType GradientType;
  GradientType ret;
  for (const auto i : make_range(N))
  {
    // Note that this can be very inefficient. Within the scalar greenGaussGradient routine we're
    // going to do value type evaluations of the array functor from scalar_functor and we will be
    // discarding all the value type evaluations other than the one corresponding to i
    ArrayComponentFunctor<T, FunctorBase<std::array<T, N>>> scalar_functor(functor, i);
    ret[i] = greenGaussGradient(elem_arg, scalar_functor, two_term_boundary_expansion, mesh);
  }

  return ret;
}

template <typename T, std::size_t N>
typename Moose::FunctorBase<std::array<T, N>>::GradientType
greenGaussGradient(const FaceArg & face_arg,
                   const Moose::FunctorBase<std::array<T, N>> & functor,
                   const bool two_term_boundary_expansion,
                   const MooseMesh & mesh)
{
  typedef typename Moose::FunctorBase<std::array<T, N>>::GradientType GradientType;
  GradientType ret;
  for (const auto i : make_range(N))
  {
    // Note that this can be very inefficient. Within the scalar greenGaussGradient routine we're
    // going to do value type evaluations of the array functor from scalar_functor and we will be
    // discarding all the value type evaluations other than the one corresponding to i
    ArrayComponentFunctor<T, FunctorBase<std::array<T, N>>> scalar_functor(functor, i);
    ret[i] = greenGaussGradient(face_arg, scalar_functor, two_term_boundary_expansion, mesh);
  }

  return ret;
}
}
}
