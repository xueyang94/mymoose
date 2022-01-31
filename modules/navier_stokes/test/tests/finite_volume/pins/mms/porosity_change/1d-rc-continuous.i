mu=1.5
rho=1.1
advected_interp_method='average'
velocity_interp_method='rc'

[Mesh]
  [mesh]
    type = CartesianMeshGenerator
    dim = 1
    dx = '1 1'
    ix = '15 15'
  []
[]

[GlobalParams]
  rhie_chow_user_object = 'rc'
[]

[UserObjects]
  [rc]
    type = PINSFVRhieChowInterpolator
    u = u
    porosity = porosity
    pressure = pressure
  []
[]

[Problem]
  error_on_jacobian_nonzero_reallocation = true
[]

[Variables]
  [u]
    type = PINSFVSuperficialVelocityVariable
    initial_condition = 1
  []
  [pressure]
    type = INSFVPressureVariable
  []
[]

[AuxVariables]
  [porosity]
    family = MONOMIAL
    order = CONSTANT
    fv = true
  []
[]

[ICs]
  [porosity_continuous]
    type = FunctionIC
    variable = porosity
    function = smooth_jump
  []
[]

[Functions]
  [smooth_jump]
    type = ParsedFunction
    value = '1 - 0.5 * 1 / (1 + exp(-30*(x-1)))'
  []

  # Generated by compute-functions-1d.py
  [exact_u]
    type = ParsedFunction
    value = 'cos((1/2)*x*pi)'
  []
  [forcing_u]
    type = ADParsedFunction
    value = '-mu*(1 - 0.5/(exp(30 - 30*x) + 1))*(-1/4*pi^2*cos((1/2)*x*pi)/(1 - 0.5/(exp(30 - 30*x) + 1)) - 15.0*pi*exp(30 - 30*x)*sin((1/2)*x*pi)/((1 - 0.5/(exp(30 - 30*x) + 1))^2*(exp(30 - 30*x) + 1)^2) - 450.0*exp(30 - 30*x)*cos((1/2)*x*pi)/((1 - 0.5/(exp(30 - 30*x) + 1))^2*(exp(30 - 30*x) + 1)^2) + 900.0*exp(60 - 60*x)*cos((1/2)*x*pi)/((1 - 0.5/(exp(30 - 30*x) + 1))^2*(exp(30 - 30*x) + 1)^3) + 450.0*exp(60 - 60*x)*cos((1/2)*x*pi)/((1 - 0.5/(exp(30 - 30*x) + 1))^3*(exp(30 - 30*x) + 1)^4)) + 15.0*mu*(-1/2*pi*sin((1/2)*x*pi)/(1 - 0.5/(exp(30 - 30*x) + 1)) + 15.0*exp(30 - 30*x)*cos((1/2)*x*pi)/((1 - 0.5/(exp(30 - 30*x) + 1))^2*(exp(30 - 30*x) + 1)^2))*exp(30 - 30*x)/(exp(30 - 30*x) + 1)^2 - pi*rho*sin((1/2)*x*pi)*cos((1/2)*x*pi)/(1 - 0.5/(exp(30 - 30*x) + 1)) + 15.0*rho*exp(30 - 30*x)*cos((1/2)*x*pi)^2/((1 - 0.5/(exp(30 - 30*x) + 1))^2*(exp(30 - 30*x) + 1)^2) + (1 - 0.5/(exp(30 - 30*x) + 1))*cos(x)'
    vars = 'mu rho'
    vals = '${mu} ${rho}'
  []
  [exact_p]
    type = ParsedFunction
    value = 'sin(x)'
  []
  [forcing_p]
    type = ParsedFunction
    value = '-1/2*pi*rho*sin((1/2)*x*pi)'
    vars = 'rho'
    vals = '${rho}'
  []
[]

[FVKernels]
  [mass]
    type = PINSFVMassAdvection
    variable = pressure
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
  []
  [mass_forcing]
    type = FVBodyForce
    variable = pressure
    function = forcing_p
  []

  [u_advection]
    type = PINSFVMomentumAdvection
    variable = u
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
    porosity = porosity
    momentum_component = 'x'
  []
  [u_viscosity]
    type = PINSFVMomentumDiffusion
    variable = u
    mu = ${mu}
    porosity = porosity
    momentum_component = 'x'
  []

  [u_pressure]
    type = PINSFVMomentumPressure
    variable = u
    pressure = pressure
    porosity = porosity
    momentum_component = 'x'
  []

  [u_forcing]
    type = INSFVBodyForce
    variable = u
    functor = forcing_u
    momentum_component = 'x'
  []
[]

[FVBCs]
  [inlet-u]
    type = INSFVInletVelocityBC
    boundary = 'left'
    variable = u
    function = 'exact_u'
  []
  [outlet_p]
    type = INSFVOutletPressureBC
    boundary = 'right'
    variable = pressure
    function = 'exact_p'
  []
[]

[Materials]
  [ins_fv]
    type = INSFVMaterial
    u = 'u'
    pressure = 'pressure'
    rho = ${rho}
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm      200                lu           NONZERO'
  line_search = 'none'

  # ksp_gmres_restart bumped to 200 for linear convergence
  nl_max_its = 100
[]

[Postprocessors]
  [inlet_p]
    type = SideAverageValue
    variable = 'pressure'
    boundary = 'left'
  []
  [outlet-u]
    type = SideIntegralVariablePostprocessor
    variable = u
    boundary = 'right'
  []
  [h]
    type = AverageElementSize
    outputs = 'console csv'
    execute_on = 'timestep_end'
  []
  [L2u]
    type = ElementL2Error
    variable = u
    function = exact_u
    outputs = 'console csv'
    execute_on = 'timestep_end'
  []
  [L2p]
    variable = pressure
    function = exact_p
    type = ElementL2Error
    outputs = 'console csv'
    execute_on = 'timestep_end'
  []
[]

[Outputs]
  exodus = true
  csv = true
[]
