# 1D strain controlled test

[GlobalParams]
  displacements = 'disp_x'
  kernel_large_kinematics = false
  material_large_kinematics = false
  constraint_types = 'strain'
  ndim = 1
  large_kinematics = false
  macro_gradient = hvar
[]

[Mesh]
  [./base]
    type = FileMeshGenerator
    file = '1d.exo'
  [../]
  [./ss]
    type = SideSetsFromPointsGenerator
    input = base
    points = '-1 0 0
               7 0 0'
    new_boundary = 'left right'
  [../]
[]

[Variables]
  [./disp_x]
  [../]
 [./hvar]
    family = SCALAR
    order = FIRST
  [../]
[]

[AuxVariables]
  [./sxx]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./exx]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

[AuxKernels]
  [./sxx]
    type = RankTwoAux
    variable = sxx
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
  [../]
  [./exx]
    type = RankTwoAux
    variable = exx
    rank_two_tensor = mechanical_strain
    index_i = 0
    index_j = 0
  [../]
[]

[UserObjects]
  [./integrator]
    type = HomogenizationConstraintIntegral
    targets = 'func_strain'
    execute_on = 'initial linear'
  [../]
[]

[Kernels]
  [./sdx]
      type = TotalLagrangianStressDivergence
      variable = disp_x
      component = 0
  [../]
[]

[ScalarKernels]
  [./enforce]
    type = HomogenizationConstraintScalarKernel
    variable = hvar
    integrator = integrator        
  [../]
[]

[Functions]
  [./func_stress]
    type = ParsedFunction
    value = '100*t'
  [../]
  [./func_strain]
    type = ParsedFunction
    value = '1.0e-2*t'
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      variable = disp_x
      auto_direction = 'x'
    [../]
  [../]

  [./centerfix_x]
    type = DirichletBC
    boundary = "fixme"
    variable = disp_x
    value = 0
  [../]
[]

[Materials]
  [./elastic_tensor_1]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 100000.0
    poissons_ratio = 0.3
    block = '1'
  [../]
  [./elastic_tensor_2]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 120000.0
    poissons_ratio = 0.21
    block = '2'
  [../]
  [./elastic_tensor_3]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 80000.0
    poissons_ratio = 0.4
    block = '3'
  [../]
  [./elastic_tensor_4]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 76000.0
    poissons_ratio = 0.11
    block = '4'
  [../]
  [./stress_base]
    type = ComputeFiniteStrainElasticStress
  [../]
  [./compute_strain]
    type = CalculateStrainLagrangianKernel
  [../]
  [./wrap_stress]
    type = WrapStressLagrangianKernel
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [./sxx]
    type = ElementAverageValue
    variable = sxx
    execute_on = 'initial timestep_end'
  [../]
  [./exx]
    type = ElementAverageValue
    variable = exx
    execute_on = 'initial timestep_end'
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'newton'
  line_search = default

  automatic_scaling = true

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  l_max_its = 2
  l_tol = 1e-14
  nl_max_its = 15
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8

  start_time = 0.0
  dt = 0.2
  dtmin = 0.2
  end_time = 1.0
[]

[Outputs]
  exodus = false
  csv = true
[]
