//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseTypes.h"

#include <set>

class InputParameters;
class MooseObject;
class FEProblemBase;
class MooseMesh;
class MortarData;
class AutomaticMortarGeneration;
class SubProblem;
class Assembly;

namespace libMesh
{
class QBase;
}

/**
 * An interface for accessing mortar mesh data
 *
 * \note mci is shorthand for mortar consumer interface, so \p _mci_problem indicates
 *       the mortar consumer interface's problem
 */
class MortarConsumerInterface
{
public:
  MortarConsumerInterface(const MooseObject * moose_object);

  static InputParameters validParams();

  /**
   * Returns the primary lower dimensional subdomain id
   */
  SubdomainID primarySubdomain() const { return _primary_subdomain_id; }

protected:
  const std::set<SubdomainID> & getHigherDimSubdomainIDs() const
  {
    return _higher_dim_subdomain_ids;
  }
  const std::set<BoundaryID> & getBoundaryIDs() const { return _boundary_ids; }

  /**
   * Retrieve the automatic mortar generation object associated with this constraint
   */
  const AutomaticMortarGeneration & amg() const;

  /**
   * Whether to interpolate the nodal normals (e.g. classic idea of evaluating field at quadrature
   * points). If this is set to false, then non-interpolated nodal normals will be used, and then
   * the _normals member should be indexed with _i instead of _qp
   */
  bool interpolateNormals() const { return _interpolate_normals; }

  /**
   * Set the normals vector
   */
  void setNormals();

  FEProblemBase & _mci_fe_problem;
  SubProblem & _mci_subproblem;
  const THREAD_ID _mci_tid;

  /// Mesh to query for boundary and subdomain ID information
  MooseMesh & _mci_mesh;

  Assembly & _mci_assembly;

  /// A reference to the mortar data object that holds all the mortar
  /// mesh information
  const MortarData & _mortar_data;

  /// Boundary ID for the secondary surface
  const BoundaryID _secondary_id;

  /// Boundary ID for the primary surface
  const BoundaryID _primary_id;

  /// Subdomain ID for the secondary surface
  const SubdomainID _secondary_subdomain_id;

  /// Subdomain ID for the primary surface
  const SubdomainID _primary_subdomain_id;

  /// the union of the secondary and primary boundary ids
  std::set<BoundaryID> _boundary_ids;

  /// the higher dimensional subdomain ids corresponding to the interior parents
  std::set<SubdomainID> _higher_dim_subdomain_ids;

  /// Whether to interpolate the nodal normals
  const bool _interpolate_normals;

  /// The locations of the quadrature points on the interior secondary elements
  const MooseArray<Point> & _phys_points_secondary;

  /// The locations of the quadrature points on the interior primary elements
  const MooseArray<Point> & _phys_points_primary;

  /// The quadrature rule on the mortar segment element
  const libMesh::QBase * const & _qrule_msm;

  /// The arbitrary quadrature rule on the lower dimensional secondary face
  const libMesh::QBase * const & _qrule_face;

  /// The secondary face lower dimensional element (not the mortar element!). The mortar element
  /// lives on the secondary side of the mortar interface and *may* correspond to \p
  /// _lower_secondary_elem under the very specific circumstance that the nodes on the primary side
  /// of the mortar interface exactly project onto the secondary side of the mortar interface. In
  /// general projection of primary nodes will split the face elements on the secondary side of the
  /// interface. It is these split elements that are the mortar segment mesh elements
  Elem const * const & _lower_secondary_elem;

  /// The element Jacobian times weights
  const std::vector<Real> & _JxW_msm;

  /// the normals
  std::vector<Point> _normals;

private:
  // Pointer to automatic mortar generation object to give constraints access to mortar geometry
  const AutomaticMortarGeneration * _amg;
};

inline const AutomaticMortarGeneration &
MortarConsumerInterface::amg() const
{
  mooseAssert(_amg, "this should have been set in the constructor");
  return *_amg;
}
