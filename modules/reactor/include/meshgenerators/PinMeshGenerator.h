//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MeshGenerator.h"

/**
 * Generates a pin-like structure for a square or hex grid with the option to be 2- or 3-D.
 */
class PinMeshGenerator : public MeshGenerator
{
public:
  static InputParameters validParams();

  PinMeshGenerator(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  ///The ReactorMeshParams object that is storing the reactor global information for this reactor geometry mesh
  MeshGeneratorName _reactor_params;

  ///The id number for this pin type
  subdomain_id_type _pin_type;

  ///The face-to-face size of this pin
  Real _pitch;

  ///The number of azimuthal divisions
  unsigned int _num_sectors;

  ///The outer radii of concentric rings in the pin
  std::vector<Real> _ring_radii;

  ///The inner apothem of any surrounding ducts in the pin
  std::vector<Real> _duct_radii;

  ///The number of mesh intervals in a radial division starting from the center
  std::vector<unsigned int> _intervals;

  ///The ID that will be assigned as both the block id and name as well as stored as an extra-element integer
  std::vector<std::vector<subdomain_id_type>> _region_ids;

  ///Whether this mesh should be extruded to 3-D, making it the final structure in the reactor mesh
  bool _extrude;

  ///Whether the center most elements in the pin should be quad elements as opposed to tri elements
  bool _quad_center;

  ///The ultimate name of the mesh for generation. The value of this string depends on whether the mesh is extruded.
  std::string mesh_name;

  /// The final mesh that is generated by the subgenerators;
  /// This mesh is generated by the subgenerators with only element and boundary ids changed.
  std::unique_ptr<MeshBase> * _build_mesh;
};