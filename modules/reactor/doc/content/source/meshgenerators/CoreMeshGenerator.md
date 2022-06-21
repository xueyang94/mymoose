# CoreMeshGenerator

!syntax description /Mesh/CoreMeshGenerator

## Overview

This object is designed to be used in the Reactor MeshGenerator workflow, which also consists of [`ReactorMeshParams`](ReactorMeshParams.md), [`PinMeshGenerator`](PinMeshGenerator.md), and [`AssemblyMeshGenerator`](AssemblyMeshGenerator.md).

The `CoreMeshGenerator` object generates core-like reactor geometry structures in either square or hexagonal geometries with block ID assignments and reporting (extra integer) IDs, as described in [`CartesianIDPatternedMeshGenerator`](CartesianIDPatternedMeshGenerator.md and [`HexIDPatternedMeshGenerator`](HexIDPatternedMeshGenerator.md). There is expected to only be a single `CoreMeshGenerator` in a Mesh definition.

This object automates the use and functionality of the [`CartesianIDPatternedMeshGenerator`](CartesianIDPatternedMeshGenerator.md) for cartesian  reactor geometry, [`HexIDPatternedMeshGenerator`](HexIDPatternedMeshGenerator.md) for hexagonal reactor geometry and, if extruding to three dimensions, the [`FancyExtruderGenerator'](FancyExtruderGenerator.md) through the use of the `MeshSubgenerator` functionality and supporting functionality from [`RenameBoundaryGenerator`](RenameBoundaryGenerator.md) and [`PlaneIDMeshGenerator'](PlaneIDMeshGenerator.md). In addition to the functionality of the `MeshGenerators` used, this object also automates boundary ID and name assignment.

In addition to the functionality of [`CartesianIDPatternedMeshGenerator`](CartesianIDPatternedMeshGenerator.md) or [`HexIDPatternedMeshGenerator`](HexIDPatternedMeshGenerator.md), this object allows for the definition of "empty" lattice locations using `MeshSubgenerators`. This is achieved through the use of creating "dummy" assembly meshes via [`CartesianMeshGenerator`](CartesianMeshGenerator.md) or [`HexagonConcentricCircleAdaptiveBoundaryMeshGenerator`](HexagonConcentricCircleAdaptiveBoundaryMeshGenerator.md) respectively. These assemblies are then removed after the core mesh creation via [`BlockDeletionGenerator`](BlockDeletionGenerator.md).

The `CoreMeshGenerator` object adopts much of the existing input structure of patterned MeshGenerators but also adapts to use parameters that are more accessible for reactor design.

## Reporting ID Information

The `CoreMeshGenerator` object automatically tags the mesh, if three dimensional, with the axial layers using the extra integer name "plane_id". The assemblies composing the core are also tagged via [`CartesianIDPatternedMeshGenerator`](CartesianIDPatternedMeshGenerator.md) or [`HexIDPatternedMeshGenerator`](HexIDPatternedMeshGenerator.md), using the "cell" assignment type, with the extra integer name "assembly_id" and any "dummy" assembly (identified via the [!param](/Mesh/CoreMeshGenerator/dummy_assembly_name) parameter) locations excluded.

## Exterior Boundary ID Information

The `CoreMeshGenerator` objects automatically assigns boundary information. The exterior core boundary ID is assigned with the parameter [!param](/Mesh/ReactorMeshParams/radial_boundary_id) and will have the name "outer_core".

If the core is extruded to three dimensions the top-most boundary ID must be assigned using [!param](/Mesh/ReactorMeshParams/top_boundary_id) and will have the name "top", while the bottom-most boundary must be assigned using [!param](/Mesh/ReactorMeshParams/bottom_boundary_id) and will have the name "bottom".

## Example Syntax

!listing modules/reactor/test/tests/meshgenerators/core_mesh_generator/core.i block=Mesh

!media reactor/meshgenerators/core_mesh_generator.png style=width:60%;

!syntax parameters /Mesh/CoreMeshGenerator

!syntax inputs /Mesh/CoreMeshGenerator

!syntax children /Mesh/CoreMeshGenerator
