"""
File containing functions to set up transmission rates within household, workplace,
social and random contexts for each individual
"""

"""
    assign_household_transmit_onegroup!(args)

Allocate transmission risk within household setting for each node, based on single specified normal distribution.
Transmission risks are stored in 'worker_nodes' objects, within NetworkParameters structure.

Inputs: `rng` - random number generator,
        `... parameter structures ...`,
        `household_contacts_per_node` - number of household contacts per node (not used in this function),
        `transrisk_household_group_mean` - mean probability of transmission within household (length must be 1),
        `transrisk_household_group_sd` - SD of probability of transmission within household (length must be 1) \n
Outputs: None \n
Location: assign_transmission_risk_fns.jl
"""
function assign_household_transmit_onegroup!(rng::MersenneTwister,
                                    network_params::NetworkParameters,
                                    configuration_options::ConfigurationOptions,
                                    household_contacts_per_node::Array{Int64,1},
                                    transrisk_household_group_mean::Array{Float64,1},
                                    transrisk_household_group_sd::Array{Float64,1})

    # Unpack parameters
    @unpack worker_nodes = network_params
    @unpack cmax = configuration_options

    # Get number of household groups in use
    n_transrisk_household_group_mean = length(transrisk_household_group_mean)
    n_transrisk_household_group_sd = length(transrisk_household_group_sd)

    # Throw error if more than one group
    if (n_transrisk_household_group_mean != 1)
        error("Should only be a single household group SAR mean estimate, but have found there to be $(n_transrisk_household_group_mean). Please rectify.")
    end

    # Throw error if more than one group
    if (n_transrisk_household_group_sd != 1)
        error("Should only be a single household group SAR standard deviation, but have found there to be $(n_transrisk_household_group_sd). Please rectify.")
    end

    # Construct normal distribution to sample from
    mean_val = transrisk_household_group_mean[1]
    sd_val = transrisk_household_group_sd[1]
    norm_dist = Normal(mean_val,sd_val)

    # Iterate over each individual
    for node_itr = 1:cmax
        worker_nodes[node_itr].transrisk_household = rand(rng,norm_dist)
    end

    network_params.worker_nodes = worker_nodes

    return nothing
end

"""
    assign_household_transmit_household_size!(args)

Allocate transmission risk within household setting for each node, based on specified normal distributions for different household sizes.

Must specify 4 sets of mean/SD, for households of size 2, 3, 4 and 5+. Households of size 1 have no transmission risk.
Transmission risks are stored in 'worker_nodes' objects, within NetworkParameters structure.

Inputs: `rng` - random number generator,
        `... parameter structures ...`,
        `household_contacts_per_node` - number of household contacts per node,
        `transrisk_household_group_mean` - mean probability of transmission within household, depending on household size,
        `transrisk_household_group_sd` - SD of probability of transmission within household, depending on household size \n
Outputs: None \n
Location: assign_transmission_risk_fns.jl
"""
function assign_household_transmit_household_size!(rng::MersenneTwister,
                                    network_params::NetworkParameters,
                                    configuration_options::ConfigurationOptions,
                                    household_contacts_per_node::Array{Int64,1},
                                    transrisk_household_group_mean::Array{Float64,1},
                                    transrisk_household_group_sd::Array{Float64,1})

    # Unpack parameters
    @unpack worker_nodes = network_params
    @unpack cmax = configuration_options

    # Get number of household groups in use
    n_transrisk_household_group_mean = length(transrisk_household_group_mean)
    n_transrisk_household_group_sd = length(transrisk_household_group_sd)

    # Throw error if there is not four groups
    if (n_transrisk_household_group_mean != 4)
        error("Should be four group SAR mean estimates, but have found there to be $(n_transrisk_household_group_mean). Please rectify.")
    end

    # Throw error if more than one group
    if (n_transrisk_household_group_sd != 4)
        error("Should be four group SAR standard deviations, but have found there to be $(n_transrisk_household_group_sd). Please rectify.")
    end

    # Construct normal distribution to sample from
    norm_dists = Normal.(transrisk_household_group_mean,transrisk_household_group_sd)

    # Iterate over each individual
    for node_itr = 1:cmax

        # Check number of people in the household
        hh_size = household_contacts_per_node[node_itr] + 1
        if hh_size == 1
            # Sole person in household. No household contacts, set transmission risk to zero.
            worker_nodes[node_itr].transrisk_household = 0
        elseif hh_size == 2
            # Household size two
            worker_nodes[node_itr].transrisk_household = rand(rng,norm_dists[1])
        elseif hh_size == 3
            # Household size three
            worker_nodes[node_itr].transrisk_household = rand(rng,norm_dists[2])
        elseif hh_size == 4
            # Household size four
            worker_nodes[node_itr].transrisk_household = rand(rng,norm_dists[3])
        else
            # Household size five or more
            worker_nodes[node_itr].transrisk_household = rand(rng,norm_dists[4])
        end
    end

    network_params.worker_nodes = worker_nodes

    return nothing
end

"""
    assign_workplace_static_transmit!(args)

Allocate transmission risk within static work setting for each node, based on specified normal distributions for different work sectors.

Transmission risks are stored in 'worker_nodes' objects, within NetworkParameters structure.

Inputs: `rng` - random number generator,
        `... parameter structures ...`,
        `transrisk_static_work_mean` - mean probability of transmission at workplace (static), depending on sector,
        `transrisk_static_work_sd` - SD of probability of transmission at workplace (static), depending on sector \n
Outputs: None \n
Location: assign_transmission_risk_fns.jl
"""
function assign_workplace_static_transmit!(rng::MersenneTwister,
                                    network_params::NetworkParameters,
                                    configuration_options::ConfigurationOptions,
                                    transrisk_static_work_mean::Array{Float64,1},
                                    transrisk_static_work_sd::Array{Float64,1})

    # Unpack parameters
    @unpack cmax, workertypes = configuration_options

    # Get number of sectors in use
    n_transrisk_mean = length(transrisk_static_work_mean)
    n_transrisk_sd = length(transrisk_static_work_sd)

    # Throw error if dimensions don't match
    if (n_transrisk_mean != workertypes) ||
        (n_transrisk_sd != workertypes)
        error("Dimension mismatch in transrisk parameter arrays. Please rectify.")
    end

    # Iterate over each individual
    for node_itr = 1:cmax
        worker = network_params.worker_nodes[node_itr]
        sector = worker.sector_ID
        mean_val = transrisk_static_work_mean[sector]
        sd_val = transrisk_static_work_sd[sector]
        network_params.worker_nodes[node_itr].transrisk_static_work = rand(rng,Normal(mean_val,sd_val))
    end

    return nothing
end

"""
    assign_workplace_dynamic_transmit!(args)

Allocate transmission risk within dynamic work setting for each node, based on specified normal distributions for different work sectors.

Transmission risks are stored in 'worker_nodes' objects, within NetworkParameters structure.

Inputs: `rng` - random number generator,
        `... parameter structures ...`,
        `transrisk_dynamic_work_mean` - mean probability of transmission at workplace (dynamic), depending on sector,
        `transrisk_dynamic_work_sd` - SD of probability of transmission at workplace (dynamic), depending on sector \n
Outputs: None \n
Location: assign_transmission_risk_fns.jl
"""
function assign_workplace_dynamic_transmit!(rng::MersenneTwister,
                                    network_params::NetworkParameters,
                                    configuration_options::ConfigurationOptions,
                                    transrisk_dynamic_work_mean::Array{Float64,1},
                                    transrisk_dynamic_work_sd::Array{Float64,1})

    # Unpack parameters
    @unpack cmax, workertypes = configuration_options

    # Get number of household groups in use
    n_transrisk_mean = length(transrisk_dynamic_work_mean)
    n_transrisk_sd = length(transrisk_dynamic_work_sd)

    # Throw error if dimensions don't match
    if (n_transrisk_mean != workertypes) ||
        (n_transrisk_sd != workertypes)
        error("Dimension mismatch in transrisk parameter arrays. Please rectify.")
    end

    # Iterate over each individual
    for node_itr = 1:cmax
        worker = network_params.worker_nodes[node_itr]
        sector = worker.sector_ID
        mean_val = transrisk_dynamic_work_mean[sector]
        sd_val = transrisk_dynamic_work_sd[sector]
        network_params.worker_nodes[node_itr].transrisk_dynamic_work = rand(rng,Normal(mean_val,sd_val))
    end

    return nothing
end

"""
    assign_social_transmit!(args)

Allocate transmission risk within social setting for each node, based on single specified normal distribution.

Transmission risks are stored in 'worker_nodes' objects, within NetworkParameters structure.

Inputs: `rng` - random number generator,
        `... parameter structures ...`,
        `transrisk_social_mean` - mean probability of social transmission,
        `transrisk_social_sd` - SD of probability of social transmission \n
Outputs: None \n
Location: assign_transmission_risk_fns.jl
"""
function assign_social_transmit!(rng::MersenneTwister,
                                    network_params::NetworkParameters,
                                    configuration_options::ConfigurationOptions,
                                    transrisk_social_mean::Float64,
                                    transrisk_social_sd::Float64)

    # Unpack parameters
    @unpack cmax = configuration_options

    # Iterate over each individual
    for node_itr = 1:cmax
        network_params.worker_nodes[node_itr].transrisk_social = rand(rng,Normal(transrisk_social_mean,transrisk_social_sd))
    end

    return nothing
end

"""
    assign_random_transmit!(args)

Allocate transmission risk within random setting for each node, based on single specified normal distribution.

Transmission risks are stored in 'worker_nodes' objects, within NetworkParameters structure.

Inputs: `rng` - random number generator,
        `... parameter structures ...`,
        `transrisk_random_mean` - mean probability of random transmission,
        `transrisk_random_sd` - SD of probability of random transmission \n
Outputs: None \n
Location: assign_transmission_risk_fns.jl
"""
function assign_random_transmit!(rng::MersenneTwister,
                                    network_params::NetworkParameters,
                                    configuration_options::ConfigurationOptions,
                                    transrisk_random_mean::Float64,
                                    transrisk_random_sd::Float64)

    # Unpack parameters
    @unpack cmax = configuration_options

    # Iterate over each individual
    for node_itr = 1:cmax
        network_params.worker_nodes[node_itr].transrisk_random = rand(rng,Normal(transrisk_random_mean,transrisk_random_sd))
    end

    return nothing
end
