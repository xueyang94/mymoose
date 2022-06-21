//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralUserObject.h"
#include "TaggingInterface.h"
#include "BlockRestrictable.h"
#include "ADReal.h"
#include "MooseTypes.h"
#include "CellCenteredMapFunctor.h"
#include "VectorComponentFunctor.h"
#include "libmesh/vector_value.h"
#include "libmesh/id_types.h"
#include "libmesh/stored_range.h"
#include <unordered_map>
#include <set>
#include <unordered_set>

class MooseMesh;
class INSFVVelocityVariable;
class INSFVPressureVariable;
namespace libMesh
{
class Elem;
class MeshBase;
}

/**
 * This user-object gathers 'a' (on-diagonal velocity coefficients) data. Having the gathered 'a'
 * data, this object is responsible for the computation of the Rhie-Chow velocity, which can be used
 * in advection kernels and postprocessors. This class also supports computation of an average face
 * velocity although this is generally not encouraged as it will lead to a checkerboard in the
 * pressure field
 */
class INSFVRhieChowInterpolator : public GeneralUserObject,
                                  public TaggingInterface,
                                  public BlockRestrictable
{
public:
  static InputParameters validParams();
  INSFVRhieChowInterpolator(const InputParameters & params);

  /**
   * API that momentum residual objects that have on-diagonals for velocity call
   * @param The element we are adding 'a' coefficient data for
   * @param component The velocity component we are adding 'a' coefficient data for
   * @param value The value of 'a' that we are adding
   */
  void addToA(const libMesh::Elem * elem, unsigned int component, const ADReal & value);

  /**
   * Retrieve a face velocity
   * @param m The velocity interpolation method. This is either RhieChow or Average. RhieChow is
   * recommended as it avoids checkboards in the pressure field
   * @param fi The face that we wish to retrieve the velocity for
   * @param tid The thread ID
   * @return The face velocity
   */
  VectorValue<ADReal>
  getVelocity(Moose::FV::InterpMethod m, const FaceInfo & fi, THREAD_ID tid) const;

  void initialSetup() override;
  void residualSetup() override;
  void meshChanged() override;

  void initialize() override final;
  void execute() override final;
  void finalize() override final;

protected:
  /**
   * A virtual method that allows us to only implement getVelocity once for free and porous flows
   */
  virtual const Moose::FunctorBase<ADReal> & epsilon(THREAD_ID tid) const;

  /**
   * perform the setup of this object
   */
  void insfvSetup();

  /// The \p MooseMesh that this user object operates on
  MooseMesh & _moose_mesh;

  /// The \p libMesh mesh that this object acts on
  const libMesh::MeshBase & _mesh;

  /// The dimension of the mesh, e.g. 3 for hexes and tets, 2 for quads and tris
  const unsigned int _dim;

  /// A functor for computing the (non-RC corrected) velocity
  std::vector<std::unique_ptr<PiecewiseByBlockLambdaFunctor<ADRealVectorValue>>> _vel;

  /// The interpolation method to use for the velocity
  Moose::FV::InterpMethod _velocity_interp_method;

  /// The thread 0 copy of the pressure variable
  INSFVPressureVariable * const _p;

  /// The thread 0 copy of the x-velocity variable
  INSFVVelocityVariable * const _u;

  /// The thread 0 copy of the y-velocity variable (null if the problem is 1D)
  INSFVVelocityVariable * const _v;

  /// The thread 0 copy of the z-velocity variable (null if the problem is not 3D)
  INSFVVelocityVariable * const _w;

  /// All the thread copies of the pressure variable
  std::vector<MooseVariableFVReal *> _ps;

  /// All the thread copies of the x-velocity variable
  std::vector<MooseVariableFVReal *> _us;

  /// All the thread copies of the y-velocity variable
  std::vector<MooseVariableFVReal *> _vs;

  /// All the thread copies of the z-velocity variable
  std::vector<MooseVariableFVReal *> _ws;

  /// All the active and elements local to this process that exist on this object's subdomains
  std::unique_ptr<ConstElemRange> _elem_range;

  /// The subdomain ids this object operates on
  const std::set<SubdomainID> _sub_ids;

  /// A map from element IDs to 'a' coefficient data
  CellCenteredMapFunctor<ADRealVectorValue, std::unordered_map<dof_id_type, ADRealVectorValue>> _a;

  /**
   * @name 'a' component functors
   * These vector component functors are not used anywhere within this class but they can be used
   * for outputting, to auxiliary variables, the on-diagonal 'a' coefficients for use in
   * visualization or transfer to other applications
   */
  ///@{
  /// The x-component of 'a'
  VectorComponentFunctor<ADReal> _ax;

  /// The y-component of 'a'
  VectorComponentFunctor<ADReal> _ay;

  /// The z-component of 'a'
  VectorComponentFunctor<ADReal> _az;
  ///@}

  /// Whether we have performed our initial setup. Ordinarily we would do this in initialSetup but
  /// there are wonky things that happen in other objects initialSetup that affect us, like how
  /// Exodus output gathers all elements to 0 on distributed mesh
  bool _initial_setup_done = false;

private:
  /**
   * Fills the _a_read data member at construction time with the appropriate functors. _a_read will
   * be used later when computing the Rhie-Chow velocity
   */
  void fillARead();

  /// The velocity variable numbers
  std::vector<unsigned int> _var_numbers;

  /// Non-local elements that we should push and pull data for across processes
  std::unordered_set<const Elem *> _elements_to_push_pull;

  /// The nonlinear system
  SystemBase & _sys;

  /// An example datum used to help communicate AD vector information in parallel
  const VectorValue<ADReal> _example;

  /// Mutex that prevents multiple threads from saving into the 'a' coefficients at the same time
  Threads::spin_mutex _a_mutex;

  /// A unity functor used in the epsilon virtual method
  const Moose::ConstantFunctor<ADReal> _unity_functor{1};

  /// A zero functor potentially used in _a_read
  const Moose::ConstantFunctor<ADReal> _zero_functor{0};

  /// A vector sized according to the number of threads that holds the 'a' data we will read from
  /// when computing the Rhie-Chow velocity
  std::vector<const Moose::FunctorBase<VectorValue<ADReal>> *> _a_read;

  /// A vector sized according to the number of threads that holds vector composites of 'a'
  /// component functors. This member is leveraged when advecting velocities are auxiliary variables
  /// and the 'a' data has been transferred from another application
  std::vector<std::unique_ptr<Moose::FunctorBase<VectorValue<ADReal>>>> _a_aux;

  /// Whether 'a' data has been provided by the user. This can happen if we are running in an
  /// application solving precursor advection, and another application has computed the fluid flow
  /// field
  bool _a_data_provided;
};

inline const Moose::FunctorBase<ADReal> & INSFVRhieChowInterpolator::epsilon(THREAD_ID) const
{
  return _unity_functor;
}

inline void
INSFVRhieChowInterpolator::addToA(const Elem * const elem,
                                  const unsigned int component,
                                  const ADReal & value)
{
  Threads::spin_mutex::scoped_lock lock(_a_mutex);

  if (elem->processor_id() != this->processor_id())
    _elements_to_push_pull.insert(elem);

  _a[elem->id()](component) += value;
}
