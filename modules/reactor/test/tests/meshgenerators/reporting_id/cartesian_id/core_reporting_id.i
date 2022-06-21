[Mesh]
  [pin1]
    type = ConcentricCircleMeshGenerator
    num_sectors = 2
    radii = '0.4 0.5'
    rings = '1 1 1'
    has_outer_square = on
    pitch = 1.26
    preserve_volumes = yes
    smoothing_max_it = 3
  []

  [pin2]
    type = ConcentricCircleMeshGenerator
    num_sectors = 2
    radii = '0.3 0.4'
    rings = '1 1 1'
    has_outer_square = on
    pitch = 1.26
    preserve_volumes = yes
    smoothing_max_it = 3
  []

  [assembly1]
    type = CartesianIDPatternedMeshGenerator
    inputs = 'pin1 pin2'
    pattern = ' 1  0  1  0;
                0  1  0  1;
                1  0  1  0;
                0  1  0  1'
    assign_type = 'cell'
    id_name = 'pin_id'
  []

  [assembly2]
    type = CartesianIDPatternedMeshGenerator
    inputs = 'pin1 pin2'
    pattern = ' 0  1  1  0;
                1  0  0  1;
                1  0  0  1;
                0  1  1  0'
    assign_type = 'cell'
    id_name = 'pin_id'
  []

  [core]
    type = CartesianIDPatternedMeshGenerator
    inputs = 'assembly1 assembly2'
    pattern = '0  1;
               1  0'
    assign_type = 'cell'
    id_name = 'assembly_id'
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
