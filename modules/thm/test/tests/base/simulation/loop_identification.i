# This test tests the loop identification function, which creates a map of component
# names to a loop name. "Loops" are defined to be sets of components which are
# physically connected - heat exchanger connections do not constitute physical
# connections in this sense. Note that this test is not meant to actually perform
# any physical computations, so dummy values are provided for the required parameters.
#
# The test configuration for this test is the following:
#
# pipe1 -> corechannel:pipe -> pipe2 -> hx:primary -> pipe1
#       j1                  j2       j3                 j4
#
# inlet -> hx:secondary -> outlet
#
# This test uses the command-line option "--print-component-loops" to print out
# the lists of components in each loop, with the desired output being the
# following:
#
# Loop 1:
#
#   corechannel:pipe
#   hx:primary
#   j1
#   j2
#   j3
#   j4
#   pipe1
#   pipe2
#
# Loop 2:
#
#   hx:secondary
#   inlet
#   outlet

[GlobalParams]
  closures = simple

  initial_p = 1e6
  initial_T = 300
  initial_vel = 0

  velocity_relaxation = false
  pressure_relaxation = false
  interface_transfer = false
  wall_mass_transfer = false

  initial_p_liquid = 1e6
  initial_T_liquid = 453
  initial_vel_liquid = 0
  initial_p_vapor = 1e6
  initial_T_vapor = 454
  initial_vel_vapor = 0
  initial_alpha_vapor = 0.5

  spatial_discretization = rdg
  rdg_slope_reconstruction = none
[]

[FluidProperties]
  [./fp_2phase]
    type = StiffenedGasTwoPhaseFluidProperties
  [../]
  [./fp_liquid]
    type = StiffenedGasFluidProperties
    gamma = 2.35
    q = -1167e3
    q_prime = 0
    p_inf = 1.e9
    cv = 1816
  [../]
[]

[HeatStructureMaterials]
  [./hx:wall]
    type = SolidMaterialProperties
    k = 1
    Cp = 1
    rho = 1
  [../]
[]

[Components]
  # PRIMARY LOOP

  [./pipe1]
    type = FlowChannel1Phase
    position = '0 0 0'
    orientation = '1 0 0'
    length = 1
    n_elems = 1
    A = 1
    f = 1
    fp = fp_liquid
  [../]
  [./j1]
    type = JunctionOneToOne
    connections = 'pipe1:out corechannel:in'
  [../]
  [./corechannel]
    type = FlowChannel1Phase
    position = '0 0 0'
    orientation = '1 0 0'
    length = 1
    n_elems = 1
    A = 1
    f = 1
    fp = fp_liquid
  [../]
  [./j2]
    type = JunctionOneToOne
    connections = 'corechannel:out pipe2:in'
  [../]
  [./pipe2]
    type = FlowChannel1Phase
    position = '0 0 0'
    orientation = '1 0 0'
    length = 1
    n_elems = 1
    A = 1
    f = 1
    fp = fp_liquid
  [../]
  [./j3]
    type = JunctionOneToOne
    connections = 'pipe2:out hx:primary:in'
  [../]
  [./hx:primary]
    type = FlowChannel1Phase
    position = '0 1 0'
    orientation = '1 0 0'
    length = 1
    n_elems = 1
    A = 1
    f = 1
    fp = fp_liquid
  [../]
  [./j4]
    type = JunctionOneToOne
    connections = 'hx:primary:out pipe1:in'
  [../]

  # HEAT EXCHANGER

  [./hs]
    type = HeatStructurePlate
    position = '0 0 0'
    orientation = '1 0 0'
    length = 1
    n_elems = 1
    materials = hx:wall
    n_part_elems = 1
    names = 0
    widths = 1
    depth = 1
    initial_T = 300
  [../]
  [./ht_primary]
    type = HeatTransferFromHeatStructure1Phase
    hs = hs
    flow_channel = hx:primary
    hs_side = outer
    Hw = 0
  [../]
  [./ht_secondary]
    type = HeatTransferFromHeatStructure2Phase
    hs = hs
    flow_channel = hx:secondary
    hs_side = inner
    Hw_liquid = 0
    Hw_vapor = 0
  [../]

  # SECONDARY LOOP

  [./inlet]
    type = SolidWall
    input = 'hx:secondary:out'
  [../]
  [./hx:secondary]
    type = FlowChannel2Phase
    position = '0 0 0'
    orientation = '1 0 0'
    length = 1
    n_elems = 1
    A = 1
    f = 1
    f_interface = 1
    fp = fp_2phase
  [../]
  [./outlet]
    type = SolidWall
    input = 'hx:secondary:in'
  [../]
[]

[Problem]
  solve = false
[]

[Executioner]
  type = Transient
  num_steps = 1
[]

[Outputs]
  [./console]
    type = Console
    system_info = ''
    enable = false
  [../]
[]
