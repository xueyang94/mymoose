[Mesh]
  [pin1]
    type = PolygonConcentricCircleMeshGenerator
    preserve_volumes = true
    ring_radii = 0.4
    ring_intervals = 1
    background_intervals = 1
    num_sides = 6
    num_sectors_per_side = '2 2 2 2 2 2'
    polygon_size = 0.5
  []
  [pin2]
    type = PolygonConcentricCircleMeshGenerator
    preserve_volumes = true
    ring_radii = 0.4
    ring_intervals = 1
    background_intervals = 1
    num_sides = 6
    num_sectors_per_side = '2 2 2 2 2 2'
    polygon_size = 0.5
  []
  [assembly1]
    type = HexIDPatternedMeshGenerator
    inputs = 'pin1 pin2'
    pattern_boundary = hexagon
    pattern = '  1 0 1;
                0 0 0 0;
               1 0 1 0 1;
                0 0 0 0;
                 1 0 1'
    hexagon_size = 2.6
    duct_sizes = '2.4 2.5'
    duct_intervals = '1 1'
    id_name = 'pin_id'
    assign_type = 'cell'
  []
  [assembly2]
    type = HexIDPatternedMeshGenerator
    inputs = 'pin1 pin2'
    pattern_boundary = hexagon
    pattern = '  0 0 0;
                0 1 1 0;
               0 1 0 1 0;
                0 1 1 0;
                 0 0 0'
    hexagon_size = 2.6
    duct_sizes = '2.4 2.5'
    duct_intervals = '1 1'
    id_name = 'pin_id'
    assign_type = 'cell'
  []
  [core]
    type = HexIDPatternedMeshGenerator
    inputs = 'assembly1 assembly2'
    pattern_boundary = none
    pattern = '1 1;
              1 0 1;
               1 1'
    generate_core_metadata = true
    id_name = 'assembly_id'
    assign_type = 'cell'
  []
[]

[Executioner]
  type = Steady
[]

[Problem]
  solve = false
[]

[AuxVariables]
  [pin_id]
    family = MONOMIAL
    order = CONSTANT
  []
  [assembly_id]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [set_pin_id]
    type = ExtraElementIDAux
    variable = pin_id
    extra_id_name = pin_id
  []
  [set_assembly_id]
    type = ExtraElementIDAux
    variable = assembly_id
    extra_id_name = assembly_id
  []
[]

[Outputs]
  exodus = true
  execute_on = timestep_end
[]
