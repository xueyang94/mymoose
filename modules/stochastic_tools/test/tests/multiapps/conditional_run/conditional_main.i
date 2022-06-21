[StochasticTools]
[]

[Samplers]
  [cart]
    type = CartesianProduct
    linear_space_items = '1 1 10'
  []
[]

[MultiApps]
  [runner]
    type = SamplerFullSolveMultiApp
    sampler = cart
    input_files = 'sub.i'
    mode = batch-reset
    should_run_reporter = conditional/need_sample
  []
[]

[Transfers]
  [data]
    type = SamplerReporterTransfer
    from_multi_app = runner
    sampler = cart
    from_reporter = 'average/value'
    stochastic_reporter = conditional
  []
[]

[Reporters]
  [conditional]
    type = ConditionalSampleReporter
    sampler = cart
    default_value = 1
    function = 'val >= t'
    sampler_vars = 'val'
    sampler_var_indices = '0'
    parallel_type = ROOT
    execute_on = 'initial timestep_begin'
  []
[]

[Executioner]
  type = Transient
  num_steps = 10
[]

[Outputs]
  execute_on = timestep_end
  [out]
    type = JSON
    execute_system_information_on = none
  []
[]
