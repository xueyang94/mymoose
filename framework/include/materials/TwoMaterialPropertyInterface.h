//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MaterialPropertyInterface.h"

// Forward Declarations
class MaterialData;
class TwoMaterialPropertyInterface : public MaterialPropertyInterface
{
public:
  TwoMaterialPropertyInterface(const MooseObject * moose_object,
                               const std::set<SubdomainID> & blocks_ids,
                               const std::set<BoundaryID> & boundary_ids);

  static InputParameters validParams();

  /**
   * Retrieve the property deduced from the name \p name
   */
  template <typename T>
  const MaterialProperty<T> & getNeighborMaterialProperty(const std::string & name);

  /**
   * Retrieve the property named "name" without any deduction
   */
  template <typename T>
  const MaterialProperty<T> & getNeighborMaterialPropertyByName(const std::string & name);

  /**
   * Retrieve the ADMaterialProperty named "name"
   */
  template <typename T>
  const ADMaterialProperty<T> & getNeighborADMaterialProperty(const std::string & name);

  template <typename T>
  const MaterialProperty<T> & getNeighborMaterialPropertyOld(const std::string & name);

  template <typename T>
  const MaterialProperty<T> & getNeighborMaterialPropertyOlder(const std::string & name);

  /**
   * Retrieve the neighbor material property whether AD or not
   */
  template <typename T, bool is_ad, typename std::enable_if<is_ad, int>::type = 0>
  const ADMaterialProperty<T> & getGenericNeighborMaterialProperty(const std::string & name)
  {
    return getNeighborADMaterialProperty<T>(name);
  }
  template <typename T, bool is_ad, typename std::enable_if<!is_ad, int>::type = 0>
  const MaterialProperty<T> & getGenericNeighborMaterialProperty(const std::string & name)
  {
    return getNeighborMaterialProperty<T>(name);
  }

protected:
  std::shared_ptr<MaterialData> _neighbor_material_data;
};

template <typename T>
const MaterialProperty<T> &
TwoMaterialPropertyInterface::getNeighborMaterialProperty(const std::string & name)
{
  // Check if the supplied parameter is a valid input parameter key
  std::string prop_name = deducePropertyName(name);

  // Check if it's just a constant
  const MaterialProperty<T> * default_property = defaultMaterialProperty<T>(prop_name);
  if (default_property)
    return *default_property;

  return this->getNeighborMaterialPropertyByName<T>(prop_name);
}

template <typename T>
const MaterialProperty<T> &
TwoMaterialPropertyInterface::getNeighborMaterialPropertyByName(const std::string & name)
{
  checkExecutionStage();
  checkMaterialProperty(name);

  // mark property as requested
  markMatPropRequested(name);

  // Update the boolean flag.
  _get_material_property_called = true;

  _material_property_dependencies.insert(_material_data->getPropertyId(name));

  return _neighbor_material_data->getProperty<T>(name);
}

template <typename T>
const ADMaterialProperty<T> &
TwoMaterialPropertyInterface::getNeighborADMaterialProperty(const std::string & name)
{
  // Check if the supplied parameter is a valid input parameter key
  std::string prop_name = deducePropertyName(name);

  // Check if it's just a constant
  const ADMaterialProperty<T> * default_property = defaultADMaterialProperty<T>(prop_name);
  if (default_property)
    return *default_property;
  else
  {
    _material_property_dependencies.insert(_material_data->getPropertyId(prop_name));
    return _neighbor_material_data->getADProperty<T>(prop_name);
  }
}

template <typename T>
const MaterialProperty<T> &
TwoMaterialPropertyInterface::getNeighborMaterialPropertyOld(const std::string & name)
{
  // Check if the supplied parameter is a valid input parameter key
  std::string prop_name = deducePropertyName(name);

  // Check if it's just a constant
  const MaterialProperty<T> * default_property = defaultMaterialProperty<T>(prop_name);
  if (default_property)
    return *default_property;
  else
    return _neighbor_material_data->getPropertyOld<T>(prop_name);
}

template <typename T>
const MaterialProperty<T> &
TwoMaterialPropertyInterface::getNeighborMaterialPropertyOlder(const std::string & name)
{
  // Check if the supplied parameter is a valid input parameter key
  std::string prop_name = deducePropertyName(name);

  // Check if it's just a constant
  const MaterialProperty<T> * default_property = defaultMaterialProperty<T>(prop_name);
  if (default_property)
    return *default_property;
  else
    return _neighbor_material_data->getPropertyOlder<T>(prop_name);
}
