a=1.1

[Mesh]
  [./gen_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 2
    xmax = 3
    ymin = 0
    ymax = 1
    nx = 2
    ny = 2
  [../]
[]

[Problem]
  coord_type = 'RZ'
[]

[Variables]
  [./v]
    family = MONOMIAL
    order = CONSTANT
    fv = true
    initial_condition = 1
  [../]
[]

[FVKernels]
  [./advection]
    type = FVAdvection
    variable = v
    velocity = '${a} ${a} 0'
    advected_interp_method = 'average'
  [../]
  [reaction]
    type = FVReaction
    variable = v
  []
  [body_v]
    type = FVBodyForce
    variable = v
    function = 'forcing'
  []
[]

[FVBCs]
  [advection]
    type = FVAdvectionFunctionBC
    boundary = 'left right top bottom'
    exact_solution = 'exact'
    variable = v
    velocity = '${a} ${a} 0'
    advected_interp_method = 'average'
  []
[]

[Functions]
[exact]
  type = ParsedFunction
  value = 'sin(x)*cos(y)'
[]
[forcing]
  type = ParsedFunction
  value = '-a*sin(x)*sin(y) + sin(x)*cos(y) + (x*a*cos(x)*cos(y) + a*sin(x)*cos(y))/x'
  vars = 'a'
  vals = '${a}'
[]
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -sub_pc_factor_shift_type -sub_pc_type'
  petsc_options_value = 'asm      NONZERO                   lu'
[]

[Outputs]
  exodus = true
[]

[Postprocessors]
  [./error]
    type = ElementL2Error
    variable = v
    function = exact
    outputs = 'console'    execute_on = 'timestep_end'
  [../]
  [h]
    type = AverageElementSize
    outputs = 'console'    execute_on = 'timestep_end'
  []
[]
