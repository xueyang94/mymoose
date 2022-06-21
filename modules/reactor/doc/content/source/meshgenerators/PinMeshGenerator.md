# PinMeshGenerator

!syntax description /Mesh/PinMeshGenerator

## Overview

This object is designed to be used in the Reactor MeshGenerator workflow, which also consists of [`ReactorMeshParams`](ReactorMeshParams.md), [`AssemblyMeshGenerator`](AssemblyMeshGenerator.md), and [`CoreMeshGenerator`](CoreMeshGenerator.md).

The `PinMeshGenerator` object generates square or hexagonal reactor geometry pin cell structures which may be combined into larger assembly structures using `AssemblyMeshGenerator`. The block IDs, external boundary ID, region IDs (e.g., materials), and reporting IDs (extra element integers identifying unique planes and pins, as described in [`CartesianIDPatternedMeshGenerator`](CartesianIDPatternedMeshGenerator.md) and [`HexIDPatternedMeshGenerator`](HexIDPatternedMeshGenerator.md) are automatically assigned once the user provides some basic information.

This pin may be extruded to three dimensions by setting [!param](/Mesh/PinMeshGenerator/extrude) to 'true', however such extruded pins cannot be used as input to `AssemblyMeshGenerator`. Instead, 2-D pins must be inputted to `AssemblyMeshGenerator` and [!param](/Mesh/AssemblyMeshGenerator/extrude) should be set to 'true' at the `AssemblyMeshGenerator` definition to extrude the assembly to 3-D.


The `PinMeshGenerator` object automates the use and functionality of the [`PolygonConcentricCircleMeshGenerator`](PolygonConcentricCircleMeshGenerator.md) and, if extruding to three dimensions, the [`FancyExtruderGenerator'](FancyExtruderGenerator.md) through the use of the `MeshSubgenerator` functionality and supporting functionality from [`TransformGenerator`](TransformGenerator.md), [`RenameBoundaryGenerator`](RenameBoundaryGenerator.md), and [`PlaneIDMeshGenerator`](PlaneIDMeshGenerator.md). In addition to the functionality of the `MeshGenerators` used, this object also automates block ID assignment and boundary ID and name assignment.

The `PinMeshGenerator` object adopts much of the existing input structure of `PolygonConcentricCircleMeshGenerator`](PolygonConcentricCircleMeshGenerator.md) but uses parameters that are more typical for reactor design.

## Block ID Information

The [!param](/Mesh/PinMeshGenerator/region_ids) parameter provides a map of subdomain_id and region_id values to assign to zones in the pin mesh. Each row in this map corresponds to a single axial layer of the pin and contains individual entries corresponding to the radial zones within the pin, starting from the centermost region and extending radially outward. The number of columns (entries in the row) should be identical to the number of rings + 1 (background region) + number of ducts. The required number of rows is dependent on the number of axial layers in the pin. For 2D pins, a single row of entries should be provided. For 3D pins, multiple rows must be provided (one for each axial layer). For 3D pins, the top row corresponds to the bottom of the pin cell.

!alert! note title=Pin block IDs are modified to match region IDs
It should be noted here that both the extra integer "region_id" and the block ID of the resultant pin elements will be modified to match the same value as specified by [!param](/Mesh/PinMeshGenerator/region_ids).
!alert-end!

The region_ids parameter entries can conveniently be selected to match material ids to be assigned to each region of the problem. Using the same value in multiple entries of the region_id parameter will effectively assign elements in multiple zones to the same subdomain_id and the same region_id. For meshes with all quadrilateral elements, this approach does not present any conflicts. However, if [!param](/Mesh/PinMeshGenerator/quad_center_elements) is set to false, the innermost radial zone of the pin is discretized into triangular elements, and therefore only one radial meshing interval may be defined for this radial zone using the [!param](/Mesh/PinMeshGenerator/mesh_intervals) parameter.  This ensures that a single region id does not correspond to a mesh discretization with both triangle and quadrilateral mesh elements, which will occur if more than one meshing interval is applied to the centermost pin zone with triangular elements. Moreover, the subdomain ID associated with any zone of triangular elements should not be shared with another zone containing quadrilateral elements, otherwise MOOSE will error out.

## Reporting ID Information

The `PinMeshGenerator` object also tags the mesh elements with the extra integer reporting ID named "region_id".

The `PinMeshGenerator` object also automatically tags the mesh with the [!param](/Mesh/PinMeshGenerator/pin_type) using the extra integer name "pin_type_id" and, if extruded, the axial layers using the extra integer name "plane_id".

## Exterior Boundary ID Information

The `PinMeshGenerator` object automatically assigns boundary information derived from the [!param](/Mesh/PinMeshGenerator/pin_type) parameter. The exterior pin boundary is assigned the ID equal to 20000 + the pin type ID and is named "outer_pin_<pin_type_id>" (for example a pin with a pin type ID of 1 will have a boundary ID of 20001 and boundary name of "outer_pin_1").

If the pin is extruded to three dimensions the top-most boundary ID must be assigned using [!param](/Mesh/ReactorMeshParams/top_boundary_id) and will have the name "top", while the bottom-most boundary must be assigned using [!param](/Mesh/ReactorMeshParams/bottom_boundary_id) and will have the name "bottom".

## Example Syntax

!listing modules/reactor/test/tests/meshgenerators/pin_mesh_generator/pin_only.i block=Mesh

!media reactor/meshgenerators/pin_mesh_generator.png style=width:60%;

!syntax parameters /Mesh/PinMeshGenerator

!syntax inputs /Mesh/PinMeshGenerator

!syntax children /Mesh/PinMeshGenerator
