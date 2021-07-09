"""
File containing functions used to set and change properties of the network,
before the simulations begin.
"""

"""
    load_defaults()

Loads necessary parameter structures according to their default settings and input constants.

Inputs: None (uses global constants only) \n
Outputs: `... parameter structures ...` - with specified defaults \n
Location: configuration_fns.jl
"""
function load_defaults()

    configuration_options = ConfigurationOptions(rng_seed = RNGSEED,
                                                countfinal = COUNTFINAL,
                                                endtime = ENDTIME,
                                                cmax = CMAX,
                                                workertypes = WORKERTYPES)

    network_params,
    workplace_generation_params = find_network_parameters(WORKERTYPES, CMAX)

    infection_params = InfectionParameters()
    #If there are not 41 sectors, use an arbitrary average for transmission risks
    #(THIS SHOULD NOT BE RELIED ON FOR ACCURATE RESULTS)
    if WORKERTYPES != 41
      infection_params.transrisk_static_work_mean = repeat([mean(infection_params.transrisk_static_work_mean)], WORKERTYPES)
      infection_params.transrisk_dynamic_work_mean = repeat([mean(infection_params.transrisk_dynamic_work_mean)], WORKERTYPES)
      infection_params.transrisk_static_work_sd = repeat([mean(infection_params.transrisk_static_work_sd)], WORKERTYPES)
      infection_params.transrisk_dynamic_work_sd = repeat([mean(infection_params.transrisk_dynamic_work_sd)], WORKERTYPES)
    end
    contact_tracing_params = ContactTracingParameters()

    workplace_closure_params = WorkplaceClosureParameters()

    return configuration_options::ConfigurationOptions,
            network_params::NetworkParameters,
            contact_tracing_params::ContactTracingParameters,
            infection_params::InfectionParameters,
            workplace_generation_params::WorkplaceGenerationParameters,
            workplace_closure_params::WorkplaceClosureParameters
end

"""
    load_configurations(configuration_name::String)

Load configuration variables according to previously defined configuration. New configurations must define n_configs > 0.

Inputs: `configuration_name` - name of configuration to be loaded \n
Output: `configuration_variables` - ConfigurationVariables structure \n
Location: configuration_fns.jl
"""
function load_configurations(configuration_name::String)

    configuration_variables = ConfigurationVariables()

    if configuration_name=="default"
        configuration_variables.n_configs = 10
        configuration_variables.workpercent = [0.3*ones(WORKERTYPES)]
    elseif configuration_name=="synchronised_changedays"
        configuration_variables.toff = [0:4;]
        configuration_variables.sameday = [0]
        configuration_variables.ton = [0]
        configuration_variables.n_configs = length(configuration_variables.toff)
    elseif configuration_name=="variable_consecutive_changedays"
        configuration_variables.toff = [0:4;]
        configuration_variables.sameday = [1]
        configuration_variables.ton = [0]
        configuration_variables.n_configs = length(configuration_variables.toff)
    elseif configuration_name=="variable_nonconsecutive_changedays"
        configuration_variables.toff = [0:4;]
        configuration_variables.sameday = [2]
        configuration_variables.ton = [0]
        configuration_variables.n_configs = length(configuration_variables.toff)
    elseif configuration_name == "groups_off"
        if WORKERTYPES!=4
            error("For groups off, need 4 workertypes")
        else
            configuration_variables.toff = ones(Int64,4)*4
            configuration_variables.ton = zeros(Int64,length(toff))
            configuration_variables.sameday = ones(Int64,length(toff))*0
            configuration_variables.workpercent = [[0,  1,  1,  1]
                                                    [1,  0,  1,  1]
                                                    [1,  1,  0,  1]
                                                    [1,  1,  1,  0]]
            configuration_variables.n_configs = length(configuration_variables.workpercent)
        end
    elseif configuration_name=="weeks_on"
        configuration_variables.toff = [1, 1, 2, 2, 2, 2]
        configuration_variables.ton = [1, 1, 2, 2, 1, 1]
        configuration_variables.sameday = [3, 4, 3, 4, 3, 4]
        configuration_variables.n_configs = length(configuration_variables.toff)
    elseif configuration_name=="amount_backwards_CT"
        configuration_variables.isolation = [1]
        configuration_variables.perform_CT_from_infector = [true]
        configuration_variables.contact_tracing_active = [true]
        configuration_variables.prob_backwards_CT = [0:0.05:0.5;]
        configuration_variables.n_configs = length(configuration_variables.prob_backwards_CT)
    elseif configuration_name=="run_one_run"
        configuration_variables.n_configs = 1
    elseif configuration_name=="workplace_CT_threshold"
        configuration_variables.workplace_closure_active = [true]
        configuration_variables.workplace_CT_threshold = [0:0.05:0.5;]
        configuration_variables.n_configs = length(configuration_variables.workplace_CT_threshold)
    elseif configuration_name=="CT_engagement"
        configuration_variables.isolation = [1]
        configuration_variables.contact_tracing_active = [true]
        configuration_variables.CT_engagement = [0:0.1:1;]
        configuration_variables.n_configs = length(configuration_variables.CT_engagement)
    elseif configuration_name=="dynamic_social_timeframe"
        configuration_variables.network_generation_method_dynamic_social = ["repeated"]
        configuration_variables.group_limit = [6]
        configuration_variables.transrisk_social_mean = [1.]
        configuration_variables.transrisk_social_sd = [0.]
        configuration_variables.dynamic_time_frame = [1,2,3,5,7,10,14]
        configuration_variables.n_configs = length(configuration_variables.dynamic_time_frame)
    elseif configuration_name=="dynamic_social_timeframe_groups"
        configuration_variables.network_generation_method_dynamic_social = ["groups"]
        configuration_variables.group_limit = [6]
        configuration_variables.transrisk_social_mean = [1]
        configuration_variables.transrisk_social_sd = [0]
        configuration_variables.dynamic_time_frame = [1,2,3,5,7,10,14]
        configuration_variables.n_configs = length(configuration_variables.dynamic_time_frame)
    elseif configuration_name=="social_group_size"
        configuration_variables.network_generation_method_dynamic_social = ["repeated"]
        configuration_variables.group_limit = [6]
        configuration_variables.max_contacts_social = [3,6,12,24,48,100]
        configuration_variables.n_configs = length(configuration_variables.max_contacts_social)
    elseif configuration_name=="rule_of_n"
        configuration_variables.network_generation_method_dynamic_social = ["fixed_daily_degree"]
        configuration_variables.group_limit = [2,3,4,6,12]
        configuration_variables.n_configs = length(configuration_variables.group_limit)
    elseif configuration_name=="rule_of_n_weekly"
        configuration_variables.network_generation_method_dynamic_social = ["fixed_daily_degree"]
        configuration_variables.dynamic_time_frame = [7]
        configuration_variables.group_limit = [2,3,4,6,12]
        configuration_variables.n_configs = length(configuration_variables.group_limit)
    elseif configuration_name=="rule_of_n_monthly"
        configuration_variables.network_generation_method_dynamic_social = ["fixed_daily_degree"]
        configuration_variables.dynamic_time_frame = [30]
        configuration_variables.group_limit = [2,3,4,6,12]
        configuration_variables.n_configs = length(configuration_variables.group_limit)
    elseif configuration_name=="lockdown_adherence"
    elseif configuration_name=="change_cmax"
        configuration_variables.cmax = [10000:10000:50000;]
        configuration_variables.n_configs = length(configuration_variables.cmax)
    elseif configuration_name=="CS_workplace_no_control"
        variable_ops = [2, 5, 10, 15, 25, 50, 100, 1000]
        configuration_variables.CS_team_size = variable_ops
        configuration_variables.n_configs = length(variable_ops)
    elseif configuration_name=="RNGseed_svty"
        variable_ops = [100, 200, 300]
        configuration_variables.rng = variable_ops
        configuration_variables.n_configs = length(variable_ops)
    elseif configuration_name=="transscaling_svty"
        variable_ops = [0.6 0.6 0.6 0.6 0.6; 0.7 0.7 0.7 0.7 0.7; 0.8 0.8 0.8 0.8 0.8;
                        0.9 0.9 0.9 0.9 0.9; 1 1 1 1 1]
        configuration_variables.scaling_social = variable_ops[:,1]
        configuration_variables.scaling_work_static = variable_ops[:,2]
        configuration_variables.scaling_work_dynamic = variable_ops[:,3]
        configuration_variables.scaling_household = variable_ops[:,4]
        configuration_variables.scaling_random = variable_ops[:,5]
        configuration_variables.n_configs = length(variable_ops[:,1])
    elseif configuration_name=="adherence_svty"
        variable_ops = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
        configuration_variables.adherence = variable_ops
        configuration_variables.n_configs = length(variable_ops)
    elseif configuration_name=="clustering_svty"
        variable_ops = [0.05 0.5]
        configuration_variables.clustering = variable_ops
        configuration_variables.between_workplace_contact_probs = [[variable_ops[ii,1]] for ii=1:length(variable_ops[:,1])]
        configuration_variables.friend_of_friend_prob = variable_ops[:,2]
        configuration_variables.n_configs = length(variable_ops[:,1])
    elseif configuration_name=="workplace_clustering_svty"
        variable_ops = [0, 0.05, 0.1, 0.2, 0.3, 0.5]
        configuration_variables.between_workplace_contact_probs = [[variable_ops[ii]] for ii=1:length(variable_ops)]
        configuration_variables.n_configs = length(variable_ops)
    elseif configuration_name=="social_clustering_svty"
        variable_ops = [0, 0.25, 0.5, 0.75, 1]
        configuration_variables.friend_of_friend_prob = variable_ops
        configuration_variables.n_configs = length(variable_ops)
    elseif configuration_name=="popsize_svty"
        variable_ops = [5000, 10000, 25000, 50000, 100000]
        configuration_variables.cmax = variable_ops
        configuration_variables.n_configs = length(variable_ops)
    elseif configuration_name=="ER_no_control"
        configuration_variables.n_configs = 1
        configuration_variables.network_generation_method = "ER"
        configuration_variables.network_generation_method_dynamic_social = "ER"
    elseif configuration_name=="CS_workplace_no_control"
        configuration_variables.CS_active_flag = [true]
        configuration_variables.CS_scale_transrisk = [ones(WORKERTYPES)]
        variable_ops = [2, 5, 10, 15, 25, 50, 100, 1000]
        configuration_variables.CS_team_size = variable_ops
        configuration_variables.n_configs = length(variable_ops)
    elseif configuration_name=="workplace_clustering_svty_no_social"
        configuration_variables.transrisk_social_mean = [0]
        configuration_variables.transrisk_social_sd = [0]
        variable_ops = [0, 0.05, 0.1, 0.2, 0.3, 0.5]
        configuration_variables.between_workplace_contact_probs = [[variable_ops[ii]] for ii=1:length(variable_ops)]
        configuration_variables.n_configs = length(variable_ops)
    elseif configuration_name=="workplace_clustering_svty_no_household"
        configuration_variables.transrisk_household_group_mean = [[0,0,0,0]]
        configuration_variables.transrisk_household_group_sd = [[0,0,0,0]]
        variable_ops = [0, 0.05, 0.1, 0.2, 0.3, 0.5]
        configuration_variables.between_workplace_contact_probs = [[variable_ops[ii]] for ii=1:length(variable_ops)]
        configuration_variables.n_configs = length(variable_ops)
    elseif configuration_name=="workplace_clustering_svty_workplace_only"
        configuration_variables.transrisk_social_mean = [0]
        configuration_variables.transrisk_social_sd = [0]
        configuration_variables.transrisk_household_group_mean = [[0,0,0,0]]
        configuration_variables.transrisk_household_group_sd = [[0,0,0,0]]
        variable_ops = [0, 0.05, 0.1, 0.2, 0.3, 0.5]
        configuration_variables.between_workplace_contact_probs = [[variable_ops[ii]] for ii=1:length(variable_ops)]
        configuration_variables.n_configs = length(variable_ops)
    elseif configuration_name=="workplace_clustering_svty_no_dynamic"
        configuration_variables.transrisk_dynamic_work_mean = [zeros(Float64, WORKERTYPES)]
        configuration_variables.transrisk_dynamic_work_sd = [zeros(Float64, WORKERTYPES)]
        variable_ops = [0, 0.05, 0.1, 0.2, 0.3, 0.5]
        configuration_variables.between_workplace_contact_probs = [[variable_ops[ii]] for ii=1:length(variable_ops)]
        configuration_variables.n_configs = length(variable_ops)
    elseif configuration_name=="workplace_clustering_svty_workplace_static_only"
        configuration_variables.transrisk_social_mean = [0]
        configuration_variables.transrisk_social_sd = [0]
        configuration_variables.transrisk_household_group_mean = [[0,0,0,0]]
        configuration_variables.transrisk_household_group_sd = [[0,0,0,0]]
        configuration_variables.transrisk_dynamic_work_mean = [zeros(Float64, WORKERTYPES)]
        configuration_variables.transrisk_dynamic_work_sd = [zeros(Float64, WORKERTYPES)]
        variable_ops = [0, 0.05, 0.1, 0.2, 0.3, 0.5]
        configuration_variables.between_workplace_contact_probs = [[variable_ops[ii]] for ii=1:length(variable_ops)]
        configuration_variables.n_configs = length(variable_ops)
    elseif configuration_name=="workplace_clustering_svty_workplace_static_and_social"
        configuration_variables.transrisk_household_group_mean = [[0,0,0,0]]
        configuration_variables.transrisk_household_group_sd = [[0,0,0,0]]
        configuration_variables.transrisk_dynamic_work_mean = [zeros(Float64, WORKERTYPES)]
        configuration_variables.transrisk_dynamic_work_sd = [zeros(Float64, WORKERTYPES)]
        variable_ops = [0, 0.05, 0.1, 0.2, 0.3, 0.5]
        configuration_variables.between_workplace_contact_probs = [[variable_ops[ii]] for ii=1:length(variable_ops)]
        configuration_variables.n_configs = length(variable_ops)
    elseif configuration_name=="workplace_clustering_svty_workplace_static_and_household"
        configuration_variables.transrisk_social_mean = [0]
        configuration_variables.transrisk_social_sd = [0]
        configuration_variables.transrisk_dynamic_work_mean = [zeros(Float64, WORKERTYPES)]
        configuration_variables.transrisk_dynamic_work_sd = [zeros(Float64, WORKERTYPES)]
        variable_ops = [0, 0.05, 0.1, 0.2, 0.3, 0.5]
        configuration_variables.between_workplace_contact_probs = [[variable_ops[ii]] for ii=1:length(variable_ops)]
        configuration_variables.n_configs = length(variable_ops)
    elseif configuration_name=="workpercent_svty"
        variable_ops = [0.2,0.4,0.6,0.8,1.0]
        workpercent_segment1 = [repeat([variable_ops[i]], WORKERTYPES) for i=1:length(variable_ops)]
        if WORKERTYPES == 41
            workpercent_segment2 = [0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.5, 0.8, 0.8, 0.5, 0.5, 0.8, 0.8, 0.8, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.8, 0.3, 0.3, 0.3, 0.7, 0.5, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8, 0.8, 0.7, 0.3, 0.5, 0.8, 0.8, 0.3]   # proportion of each type of worker returning to work
            push!(workpercent_segment1, workpercent_segment2)
            configuration_variables.workpercent = workpercent_segment1
            configuration_variables.n_configs = length(variable_ops) + 1
        else
            configuration_variables.workpercent = workpercent_segment1
            configuration_variables.n_configs = length(variable_ops)
        end
    elseif configuration_name == "social_easing_config"
        configuration_variables.network_generation_method_dynamic_social = ["household_groups"]
        configuration_variables.workpercent = [0.05.*ones(Float64, WORKERTYPES)]
        configuration_variables.dynamic_time_frame = [1]
        configuration_variables.isolation = [1]
        group_limit = [3:3:15;]
        n_households_per_group = [1:15;]
        configuration_variables.group_limit = repeat(group_limit, inner = length(n_households_per_group))
        configuration_variables.n_households_per_group = repeat(n_households_per_group, outer = length(group_limit))
        configuration_variables.n_configs = length(configuration_variables.n_households_per_group)
    elseif configuration_name == "social_easing_intervention"
        configuration_variables.workpercent = [0.05.*ones(Float64, WORKERTYPES)]
        configuration_variables.isolation = [1]
        configuration_variables.dynamic_time_frame = [1]
        configuration_variables.n_configs = 10
    elseif configuration_name == "social_easing_intervention_large_households"
        configuration_variables.household_size_distribution = [[0.,0.,0.,0.,1.]]
        configuration_variables.workpercent = [0.05.*ones(Float64, WORKERTYPES)]
        configuration_variables.isolation = [1]
        configuration_variables.dynamic_time_frame = [1]
        configuration_variables.n_configs = 10
    elseif configuration_name == "social_easing_config_large_households"
        configuration_variables.household_size_distribution = [[0.,0.,0.,0.,1.]]
        configuration_variables.network_generation_method_dynamic_social = ["household_groups"]
        configuration_variables.workpercent = [0.05.*ones(Float64, WORKERTYPES)]
        configuration_variables.dynamic_time_frame = [1]
        configuration_variables.isolation = [1]
        group_limit = [3:3:15;]
        n_households_per_group = [1:15;]
        configuration_variables.group_limit = repeat(group_limit, inner = length(n_households_per_group))
        configuration_variables.n_households_per_group = repeat(n_households_per_group, outer = length(group_limit))
        configuration_variables.n_configs = length(configuration_variables.n_households_per_group)
    else
        error("Requested configuration not defined. Add configuration to load_configurations function.")
    end

    return configuration_variables::ConfigurationVariables
end

"""
    apply_config_changes!(args)

Changes values contained within parameter structures according to the current configuration.

Inputs: `config_idx` - index of current configuration,
        `... parameter structures ...` \n
Outputs: None \n
Location: configuration_fns.jl
"""
function apply_config_changes!(config_idx::Int64,
                                configuration_variables::ConfigurationVariables,
                                configuration_options::ConfigurationOptions,
                                network_params::NetworkParameters,
                                contact_tracing_params::ContactTracingParameters,
                                infection_params::InfectionParameters,
                                workplace_generation_params::WorkplaceGenerationParameters,
                                workplace_closure_params::WorkplaceClosureParameters)

    all_config_vars = fieldnames(ConfigurationVariables)     # All settings changeable within configuration
    n_config_vars = length(all_config_vars)

    for config_var_idx = 1:n_config_vars
        config_var_name = all_config_vars[config_var_idx]
        if config_var_name != :n_configs
            config_var_obj = getfield(configuration_variables, config_var_name)
            if isassigned(config_var_obj)
                if length(config_var_obj) == 1
                    if hasfield(ConfigurationOptions, config_var_name)
                        setfield!(configuration_options, config_var_name, config_var_obj[1])
                    elseif hasfield(NetworkParameters, config_var_name)
                        setfield!(network_params, config_var_name, config_var_obj[1])
                    elseif hasfield(ContactTracingParameters, config_var_name)
                        setfield!(contact_tracing_params, config_var_name, config_var_obj[1])
                    elseif hasfield(InfectionParameters, config_var_name)
                        setfield!(infection_params, config_var_name, config_var_obj[1])
                    elseif hasfield(WorkplaceGenerationParameters, config_var_name)
                        setfield!(workplace_generation_params, config_var_name, config_var_obj[1])
                    elseif hasfield(WorkplaceClosureParameters, config_var_name)
                        setfield!(workplace_closure_params, config_var_name, config_var_obj[1])
                    else
                        error("Constant configuration variable $config_var_name not found in parameter structures")
                    end
                elseif length(config_var_obj) > 1
                    if length(config_var_obj) != configuration_variables.n_configs
                        error("Mismatch in number of configurations")
                    end
                    if hasfield(ConfigurationOptions, config_var_name)
                        setfield!(configuration_options, config_var_name, config_var_obj[config_idx])
                    elseif hasfield(NetworkParameters, config_var_name)
                        setfield!(network_params, config_var_name, config_var_obj[config_idx])
                    elseif hasfield(ContactTracingParameters, config_var_name)
                        setfield!(contact_tracing_params, config_var_name, config_var_obj[config_idx])
                    elseif hasfield(InfectionParameters, config_var_name)
                        setfield!(infection_params, config_var_name, config_var_obj[config_idx])
                    elseif hasfield(WorkplaceGenerationParameters, config_var_name)
                        setfield!(workplace_generation_params, config_var_name, config_var_obj[config_idx])
                    elseif hasfield(WorkplaceClosureParameters, config_var_name)
                        setfield!(workplace_closure_params, config_var_name, config_var_obj[config_idx])
                    else
                        error("Variable configuration variable $config_var_name not found in parameter structures")
                    end
                end
            end
        end
    end

    return nothing
end

"""
    find_network_parameters(workertypes::Int64, cmax::Int64)

Load default network parameters according to number of sectors and nodes required.

Inputs: `workertypes` - number of sectors,
        `cmax` - number of nodes \n
Outputs: `network_params` - NetworkParameters structure,
         `workplace_generation_params` - WorkplaceGenerationParameters structure \n
Location: configuration_fns.jl
"""
function find_network_parameters(workertypes::Int64, cmax::Int64)

    if workertypes==2
        network_params = NetworkParameters(
            prob_workertype_contact = [0.002,0.002].*(10000/cmax),    # Scale relative to level used for 10,000 nodes
            dd_within_workplace = [5., 10.],     # mean degree distribution of workers within same workplace, for each worker type
            dynamic_conts_mean = [10., 0],         # mean number of dynamic contacts for each worker type
            dynamic_conts_sd = [2., 0])            # sd for number of dynamic contacts for each worker type
        workplace_generation_params = WorkplaceGenerationParameters(
            workforce_proportion = [0.15, 0.85],   # proportion of workforce in each worker type
            workplace_size_mean = [10., 100.],      # mean size of workplace for each worker group
            workplace_size_sd = [2., 10. ])       # sd for size of workplace for each worker group
    elseif workertypes==4
        network_params = NetworkParameters(
            prob_workertype_contact = [0.0001,0.0001,0.0001,0.0001].*(10000/cmax),    # Scale relative to level used for 10,000 nodes
            dd_within_workplace = [10., 10., 10., 10.],     # mean degree distribution of workers within same workplace, for each worker type
            dynamic_conts_mean = [0., 0., 0., 10.],         # mean number of dynamic contacts for each worker type
            dynamic_conts_sd = [0., 0., 0., 2.])            # sd for number of dynamic contacts for each worker type
        workplace_generation_params = WorkplaceGenerationParameters(
            workforce_proportion = [0.030, 0.127, 0.311, 0.532],   # proportion of workforce in each worker type
            workplace_size_mean = [6., 9., 11., 10.],      # mean size of workplace for each worker group
            workplace_size_sd = [14., 57., 122., 109. ])       # sd for size of workplace for each worker group
    elseif workertypes==6
        network_params = NetworkParameters(
            prob_workertype_contact = [0.002,0.002,0.002,0.002,0.002,0.002].*(10000/cmax),    # Scale relative to level used for 10,000 nodes
            dd_within_workplace = [5., 10., 10., 10., 20., 20.],     # mean degree distribution of workers within same workplace, for each worker type
            dynamic_conts_mean = [0., 0., 0., 10., 10., 10.],         # mean number of dynamic contacts for each worker type
            dynamic_conts_sd = [0., 0., 0., 2., 2., 2.])            # sd for number of dynamic contacts for each worker type
        workplace_generation_params = WorkplaceGenerationParameters(
            workforce_proportion = [0.025, 0.105, 0.252, 0.427, 0.072, 0.119],   # proportion of workforce in each worker type
            workplace_size_mean = [6., 9., 11., 10., 65., 33.],      # mean size of workplace for each worker group
            workplace_size_sd = [14., 57., 122., 109., 432., 305. ])       # sd for size of workplace for each worker group
    elseif workertypes==41
        network_params = NetworkParameters(
            prob_workertype_contact = ones(workertypes)*0.002.*(10000/cmax),    # Scale relative to level used for 10,000 nodes
            prob_anyworker_contact = 0.0001*(10000/cmax),    # Scale relative to level used for 10,000 nodes
            prob_social_contact = 0.002*(10000/cmax),    # Scale relative to level used for 10,000 nodes
            prob_random_contact = 0.0001*(10000/cmax),    # Scale relative to level used for 10,000 nodes
            dd_within_workplace = [5., 10., 10., 10., 5., 10., 5., 10., 10., 10., 5., 5., 10., 10., 10., 10., 10., 10., 10., 10., 5., 5., 10., 5., 5., 5., 10., 10., 10., 10., 10., 5., 10., 5., 10., 10., 10., 5., 5., 5., 10.],     # mean degree distribution of workers within same workplace, for each worker type
            dynamic_conts_mean = [0, 0, 0, 0, 2., 0, 5., 2., 10., 10., 0, 10., 10., 10., 0, 0, 0, 0, 0, 0, 5., 5., 0, 5., 0, 2., 0, 0, 2., 10., 2., 5., 10., 5., 10., 10., 0, 5., 5., 5., 0],         # mean number of dynamic contacts for each worker type
            dynamic_conts_sd = [0, 0, 0, 0, 1., 0, 2., 1., 4., 4., 0, 4., 4., 4., 0, 0, 0, 0, 0, 0, 2., 2., 0, 2., 0, 1., 0, 0, 1., 4., 1., 2., 4., 2., 4., 4., 0, 2., 2., 2., 0],            # sd for number of dynamic contacts for each worker type
            workplace_degree_distribution = find_fitted_contact_dds("workplace_fixed",workertypes),
            workplace_dynamic_degree_distribution = find_fitted_contact_dds("workplace_dynamic",workertypes))
        workplace_generation_params = WorkplaceGenerationParameters(
                                            workforce_proportion = [0.0393, 0.0035, 0.0181, 0.0958, 0.0126, 0.0627, 0.0263, 0.0582, 0.0700, 0.0206, 0.0124, 0.0032, 0.0218, 0.0674, 0.0137, 0.0273, 0.0028, 0.0188, 0.0207, 0.0957, 0.0021, 0.0063, 0.0260, 0.0033, 0.0044, 0.0184, 0.0194, 0.0064, 0.0614, 0.0399, 0.0306, 0.0342, 0.0078, 0.0016, 0.0144, 0.0027, 0.0131, 0.0039, 0.0074, 0.0010, 0.0049],   # proportion of workforce in each worker type
                                            workplace_size_mean = [5.3906, 57.6168, 34.7701, 15.5182, 19.2749, 3.7413, 6.9799, 11.5829, 6.8789, 6.1270, 12.8666, 2.9083, 24.7449, 9.9487, 5.6432, 3.3376, 6.4965, 6.3643, 4.2314, 4.1993, 10.8347, 7.2025, 16.9858, 8.0033, 10.1030, 8.6935, 3.3448, 17.3884, 28.3144, 14.6889, 59.8724, 19.5898, 5.0469, 31.3728, 10.4981, 8.1710, 12.1034, 5.7949, 3.5018, 10.8615, 4.0502],      # mean size of workplace for each worker group
                                            workplace_size_sd = [35.6780, 291.5813, 183.5026, 99.9800, 129.5788, 51.6139, 72.9294, 94.4625, 68.6747, 55.8071, 96.7104, 44.4382, 155.5328, 84.0205, 62.4778, 37.7083, 67.0943, 44.9687, 40.6361, 49.4326, 62.3247, 60.7611, 98.3582, 57.2267, 55.1600, 65.2562, 33.1987, 57.6397, 97.3243, 108.6843, 210.1163, 105.7828, 52.6915, 141.2081, 71.0224, 45.4627, 123.3399, 72.5046, 27.0021, 80.4898, 41.5465])     # sd for size of workplace for each worker group
    elseif workertypes==42
        network_params = NetworkParameters(
            prob_workertype_contact = ones(workertypes)*0.002.*(10000/cmax),    # Scale relative to level used for 10,000 nodes
            dd_within_workplace = [5., 10., 10., 10., 5., 10., 5., 10., 10., 10., 5., 5., 10., 10., 10., 10., 10., 10., 10., 10., 5., 5., 10., 5., 5., 5., 10., 10., 10., 10., 10., 5., 10., 5., 10., 10., 10., 5., 5., 5., 10., 10.],     # mean degree distribution of workers within same workplace, for each worker type
            dynamic_conts_mean = [0, 0, 0, 0, 2., 0, 5., 2., 10., 10., 0, 10., 10., 10., 0, 0, 0, 0, 0, 0, 5., 5., 0, 5., 0, 2., 0, 0, 2., 10., 2., 5., 10., 5., 10., 10., 0, 5., 5., 5., 0, 0],         # mean number of dynamic contacts for each worker type
            dynamic_conts_sd = [0, 0, 0, 0, 1., 0, 2., 1., 4., 4., 0, 4., 4., 4., 0, 0, 0, 0, 0, 0, 2., 2., 0, 2., 0, 1., 0, 0, 1., 4., 1., 2., 4., 2., 4., 4., 0, 2., 2., 2., 0, 0])            # sd for number of dynamic contacts for each worker type
        workplace_generation_params = WorkplaceGenerationParameters(
            workforce_proportion = [0.0393, 0.0035, 0.0181, 0.0958, 0.0126, 0.0627, 0.0263, 0.0582, 0.0700, 0.0206, 0.0124, 0.0032, 0.0218, 0.0674, 0.0137, 0.0273, 0.0028, 0.0188, 0.0207, 0.0957, 0.0021, 0.0063, 0.0260, 0.0033, 0.0044, 0.0184, 0.0194, 0.0064, 0.0614, 0.0399, 0.0306, 0.0342, 0.0078, 0.0016, 0.0144, 0.0027, 0.0131, 0.0039, 0.0074, 0.0010, 0.0049, 0.0001],   # proportion of workforce in each worker type
            workplace_size_mean = [5.3906, 57.6168, 34.7701, 15.5182, 19.2749, 3.7413, 6.9799, 11.5829, 6.8789, 6.1270, 12.8666, 2.9083, 24.7449, 9.9487, 5.6432, 3.3376, 6.4965, 6.3643, 4.2314, 4.1993, 10.8347, 7.2025, 16.9858, 8.0033, 10.1030, 8.6935, 3.3448, 17.3884, 28.3144, 14.6889, 59.8724, 19.5898, 5.0469, 31.3728, 10.4981, 8.1710, 12.1034, 5.7949, 3.5018, 10.8615, 4.0502, 79.6000],      # mean size of workplace for each worker group
            workplace_size_sd = [35.6780, 291.5813, 183.5026, 99.9800, 129.5788, 51.6139, 72.9294, 94.4625, 68.6747, 55.8071, 96.7104, 44.4382, 155.5328, 84.0205, 62.4778, 37.7083, 67.0943, 44.9687, 40.6361, 49.4326, 62.3247, 60.7611, 98.3582, 57.2267, 55.1600, 65.2562, 33.1987, 57.6397, 97.3243, 108.6843, 210.1163, 105.7828, 52.6915, 141.2081, 71.0224, 45.4627, 123.3399, 72.5046, 27.0021, 80.4898, 41.5465, 84.8072])       # sd for size of workplace for each worker group
    end

    return network_params, workplace_generation_params
end

"""
    find_fitted_contact_dds(context::String, workertypes::Int64)

Load default static and dynamic workplace contact degree distributions, for each sector.
Only works for 41 sectors (default).

Inputs: `context` - context for contact degree distribution,
        `workertypes` - number of sectors \n
Outputs: `workplace_degree_distribution` or `workplace_dynamic_degree_distribution` - relevant degree distribution (array of distributions) \n
Location: configuration_fns.jl
"""
function find_fitted_contact_dds(context::String, workertypes::Int64)

    if context == "workplace_fixed"
        if workertypes == 41
            work_contacts_meanlog = [1.57, 1.57, 1.57, 1.57, 1.57, 1.57, 1.69, 1.69, 1.69, 1.67, 1.67, 3.48, 1.97, 2.86, 1.73, 1.73, 1.73, 3.1, 1.56, 1.5, 1.98, 1.56, 1.73, 1.97, 1.67, 1.35, 1.73, 1.67, 3.15, 1.96, 1.96, 1.96, 1.67, 1.67, 1.11, 1.67, 1.67, 1.73, 1.5, 1.97, 1.73]
            work_contacts_sdlog = [1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.2, 1.2, 1.2, 0.97, 0.97, 0.53, 1.35, 1.78, 0.91, 0.91, 0.91, 1.44, 0.33, 1.02, 0.92, 0.33, 0.91, 1.35, 1.04, 1.35, 0.91, 1.04, 1.43, 1.25, 1.25, 1.25, 0.49, 0.49, 1.12, 0.49, 1.04, 1.09, 0.11, 1.35, 0.91]

            workplace_degree_distribution = Array{Distribution,1}(undef, length(work_contacts_meanlog))

            for ii = 1:length(workplace_degree_distribution)
                workplace_degree_distribution[ii] = Distributions.LogNormal(work_contacts_meanlog[ii],work_contacts_sdlog[ii])
            end
        end

        return workplace_degree_distribution

    elseif context == "workplace_dynamic"
        if workertypes == 41
            work_contacts_meanlog = [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 2.4, 2.4, 2.4, 3.05, 3.05, 1.76, 1.76, 1.76, 0.96, 0.96, 0.96, 0.69, 1.13, 1.11, 1.85, 1.13, 0.96, 2.85, 1.39, 1.76, 0.96, 1.39, 1.94, 1.54, 1.54, 1.54, 0.73, 0.73, 1.76, 0.73, 1.39, 0.75, 1.63, 1.76, 0.96]
            work_contacts_sdlog = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 1.2, 1.2, 1.2, 0.21, 0.21, 1.42, 1.42, 1.42, 0.99, 0.99, 0.99, 0.57, 0.81, 1.29, 1.31, 0.81, 0.99, 1.06, 1.43, 1.42, 0.99, 1.43, 1.56, 1.12, 1.12, 1.12, 0.72, 0.72, 1.42, 0.72, 1.43, 1.4, 1.59, 1.42, 0.99]

            workplace_dynamic_degree_distribution = Array{Distribution,1}(undef, length(work_contacts_meanlog))

            for ii = 1:length(workplace_dynamic_degree_distribution)
                workplace_dynamic_degree_distribution[ii] = Distributions.LogNormal(work_contacts_meanlog[ii],work_contacts_sdlog[ii])
            end
        end

        return workplace_dynamic_degree_distribution
    end
end

"""
    check_for_errors(configuration_options::ConfigurationOptions)

Check for illogical combinations in specified configuration options.

Inputs: `configuration_options` - ConfigurationOptions structure \n
Outputs: None \n
Location: configuration_fns.jl
"""
function check_for_errors(configuration_options::ConfigurationOptions)

    if (configuration_options.contact_tracing_active == true) &&
        (configuration_options.isolation == 0)
        error("Contact tracing active without isolation")
    end

    return nothing
end
