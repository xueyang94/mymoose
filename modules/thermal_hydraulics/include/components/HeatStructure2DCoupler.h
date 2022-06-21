//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "BoundaryBase.h"
#include "MeshAlignment2D2D.h"

/**
 * Couples boundaries of two 2D heat structures via a heat transfer coefficient
 */
class HeatStructure2DCoupler : public BoundaryBase
{
public:
  HeatStructure2DCoupler(const InputParameters & parameters);

  virtual void addMooseObjects() override;

protected:
  virtual void init() override;
  virtual void check() const override;

  /// Primary and secondary heat structure names
  const std::vector<std::string> _hs_names;
  /// Primary and secondary heat structure boundaries
  const std::vector<BoundaryName> _hs_boundaries;

  /// Mesh alignment
  MeshAlignment2D2D _mesh_alignment;
  /// Flag for each heat structure being HeatStructurePlate
  std::vector<bool> _is_plate;
  /// Flag for each heat structure deriving from HeatStructureCylindricalBase
  std::vector<bool> _is_cylindrical;

public:
  static InputParameters validParams();
};
