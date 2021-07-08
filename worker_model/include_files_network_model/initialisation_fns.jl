"""
File containing functions for the initialisation / reinitialisation of network components
within each replicate.
"""

"""
    initialise_settings_and_components_for_replicate!(args)

Reverse changes made to settings and network components during a single replicate, ready to start the next replicate.

Changes include:
- changes made to configuration settings and network structure via intervention are reversed
- node state vectors are reinitialised and infection status seeded
- returned-to-work status of nodes reinitialised
- asymptomatic probability and relative infectiousness are redrawn
- output storage vectors are updated
- transmission risks are redrawn
- reinitialises intervention-dependent parameters (contact tracing, workplace closures, COVID-secure workplaces)

Inputs: `count` - replicate number,
        `intervention_set_itr` - intervention set number,
        `rng` - random number generator,
        `... parameter structures ...`,
        `settings_changed_by_intervention_set` - list of settings changed by intervention during replicate,
        `timed_interventions_active` - boolean flagging if timed interventions are active during current replicates \n
Outputs: None \n
Location: initialisation_fns.jl
"""
function initialise_settings_and_components_for_replicate!(count::Int64,
                                                           intervention_set_itr::Int64,
                                                           rng::MersenneTwister,
                                                           output::SimulationOutputs,
                                                           configuration_options::ConfigurationOptions,
                                                           infection_params::InfectionParameters,
                                                           contact_tracing_params::ContactTracingParameters,
                                                           CT_vars::ContactTracingVariables,
                                                           contacts::ContactStructure,
                                                           network_params::NetworkParameters,
                                                           workplace_generation_params::WorkplaceGenerationParameters,
                                                           workplace_closure_params::WorkplaceClosureParameters,
                                                           states::NodeStates,
                                                           preintervention_settings::InterventionVariables,
                                                           contact_component_store::Array{PreinterventionComponentStore,1},
                                                           settings_changed_by_intervention_set::Array{Symbol,1},
                                                           timed_interventions_active::Bool)

    @unpack cmax, endtime, seed_initial_states_fn = configuration_options

    #  Reset settings if they have been changed
    #  reset / generate node states (except rep_inf_this_timestep),
    #  generate returned_to_work,
    #  seed non-susceptible initial states and add to output
    #  assign transmission risks
    #  reinitialise CT_vars if defined (even if not changed - since involves randomness)
    #  If intervention has changed them, reinitialise:
    #    - workplace params, workplace closure params, CS contacts

    # Reset to pre-intervention conditions
    if (timed_interventions_active == true) && (length(settings_changed_by_intervention_set) > 0)

        settings_changed_by_intervention::Array{Symbol,1} = Symbol[]
        # Reset settings to preintervention states
        change_parameters_by_intervention!(preintervention_settings,
                                         configuration_options,
                                         infection_params,
                                         contact_tracing_params,
                                         network_params,
                                         workplace_generation_params,
                                         workplace_closure_params,
                                         preintervention_settings,
                                         settings_changed_by_intervention_set,
                                         settings_changed_by_intervention)

        # Reset contacts to preintervention setting if changed during replicate
        if length(contact_component_store) > 0
            preintervention_component = contact_component_store[1]
            change_time = preintervention_component.contacts_changed_by_intervention_time
            contacts.workday_social_contacts_by_day[change_time:endtime, :] = preintervention_component.workday_social_contacts_by_day[change_time:endtime, :]
            contacts.nonworkday_social_contacts_by_day[change_time:endtime, :] = preintervention_component.nonworkday_social_contacts_by_day[change_time:endtime, :]
        end
        # Settings changed by intervention remain the same between replicates:
        # No need to reset settings_changed_by_intervention_set,
        #                  preintervention_settings and preintervention_components
    end

    """
    Initialisation / reinitialisation required for every replicate
    """
    # Initialise time series vectors
    initialise_node_states!(states,
                                configuration_options,
                                infection_params,
                                rng)

    # Initialise returned_to_work status of workers
    initialise_returned_to_work!(network_params.worker_nodes,
                                    configuration_options,
                                    workplace_generation_params,
                                    network_params,
                                    rng)


    # Draw asymptomatic probability for current replicate
    infection_params.probasymp = rand(rng, infection_params.probasymp_dist)

    # Draw relative infectiousness of an asymptomatic for current replicate
    infection_params.asymp_trans_scaling = rand(rng, infection_params.asymp_trans_scaling_dist)

    # Sets latent, asymptomatic, symptomatic, recovered nodes
    n_initial_latent::Int64,
    n_initial_asymp::Int64,
    n_initial_symp::Int64,
    n_initial_recovereds::Int64 = seed_initial_states_fn(rng,
                                                        cmax,
                                                        states,
                                                        infection_params,
                                                        configuration_options.recov_propn)


    # Update time series for latent & infecteds after assigning initial
    # infecteds
    output.numlat[1,count,intervention_set_itr] = n_initial_latent
    output.numinf[1,count,intervention_set_itr] = n_initial_asymp + n_initial_symp

    # Update prevalences
    output.prevlat[1,count,intervention_set_itr] = n_initial_latent
    output.prevasymp[1,count,intervention_set_itr] = n_initial_asymp
    output.prevsymp[1,count,intervention_set_itr] = n_initial_symp
    output.prevrec[1,count,intervention_set_itr] = n_initial_recovereds

    # Assign transmission risks
    assign_transrisk_all_settings!(rng,
                                   network_params,
                                   configuration_options,
                                   infection_params,
                                   contacts)

   # If defined, reinitialise contact tracing related variables
   if CT_vars.cmax > 0
       reinitialise_CT_vars!(CT_vars, cmax, rng,
                               contact_tracing_params,
                               states,
                               infection_params)
   end

    """
    Intervention dependent reinitialisation
    """
    # Reinitialise workplace_params
    if length(intersect([:CS_active_flag, :sector_open], settings_changed_by_intervention_set)) > 0
        reinitialise_workplace_params!(network_params.workplace_info,
                                    network_params.sector_open,
                                    configuration_options.CS_active_flag)
    end

    #Reinitialise workplace closure related parameters
    if length(intersect([:workplace_CT_threshold, :workplace_CT_memory], settings_changed_by_intervention_set)) > 0

        reinitialise_workplace_closure_params!(workplace_closure_params,
                                                    configuration_options,
                                                    network_params)
    end

    #Reinitialise CS contacts only if they were generated and changed
    if (network_params.CS_contacts_previously_generated == true) &&
     (network_params.CS_contacts_changed_by_intervention == true)
        contacts.work_contacts_same_workplace_CS = deepcopy(contacts.work_contacts_same_workplace_CS_store)
        contacts.work_contacts_same_workplace_per_node_CS = deepcopy(contacts.work_contacts_same_workplace_per_node_CS_store)
    end

    return nothing
end

"""
    initialise_node_states!(args)

Initialises objects in NodeStates structure ready for replicate.

Inputs: `... parameter structures ...`,
        `rng` - random number generator \n
Outputs: None \n
Location: initialisation_fns.jl
"""
function initialise_node_states!(states::NodeStates,
                                    configuration_options::ConfigurationOptions,
                                    infection_params::InfectionParameters,
                                    rng::MersenneTwister)
    lmul!(0,states.timelat)
    lmul!(0,states.timeinf)
    lmul!(0,states.timesymp)
    lmul!(0,states.asymp)
    lmul!(0,states.lattime)
    lmul!(0,states.hh_isolation)
    lmul!(0,states.delay_adherence)
    lmul!(0,states.acquired_infection)
    lmul!(0,states.hh_in_isolation_array)
    lmul!(0,states.symp_isolation_array)
    lmul!(0,states.CT_isolation_array)
    lmul!(0,states.daily_record_atworkplace)
    lmul!(0,states.daily_record_inisol)
    lmul!(0,states.time_to_symps)
    lmul!(0,states.infected_by)

    # Shuffle adherence heirarchy (changes which nodes will adhere)
    states.adherence_hierarchy = shuffle(rng, 1:configuration_options.cmax)

    # Construct atwork array signifiying when at workplace
    populate_atwork!(rng, states.atwork, configuration_options)

    # Generates hh_isolation, delay_adherence, time_to_symps, lattime
    set_infection_related_times!(rng, states,
                                    configuration_options,
                                    infection_params)
end

"""
    populate_atwork!(args)

Randomly initialises 'atwork' status of nodes (stored in NodeStates structure) according to configuration options.

Inputs: `rng` - random number generator,
        `atwork` - cmax x endtime array indicating work schedule of each node,
        `... parameter structures ...` \n
Outputs: None \n
Location: initialisation_fns.jl
"""
function populate_atwork!(rng::MersenneTwister,
                            atwork::Array{Int64,2},
                            configuration_options::ConfigurationOptions)

    @unpack cmax, endtime, ton, toff, sameday = configuration_options

    # Initialise array to store indicator values of whether individuals are at workplace or not each day
    lmul!(0, atwork)

    # If sameday=0, all workers are at work on the same set of consecutive days.
    # If sameday=1, workers go to work on a random set of consecutive days.
    # If sameday=2, workers go to work on the same number of days, but scattered randomly throughout the week.
    if sameday==0 # all workers go to work on the same set of consecutive days
        # iterate over each node
        for node_itr=1:cmax
            num_reps = floor(Int64,endtime/7) # how many weeks in the simulation
            # iterate over each week
            for reps_itr = 1:num_reps
                for days=ton:toff # all works work days ton to toff
                    atwork[node_itr,(reps_itr-1)*7+(days+1)] = 1
                end
            end
            # and put in the last bit, if there aren't an exact number of repetitions
            for days=ton:toff # all works work days ton to toff
                if (num_reps*7+(days+1))<=endtime
                    atwork[node_itr,num_reps*7+(days+1)] = 1
                end
            end
        end
    elseif sameday==1 # workers go to work on a random set of consecutive days

        # initialise pap
        pap = zeros(7)
        # iterate over each node
        for node_itr=1:cmax
            num_reps = floor(Int64,endtime/7) # how many weeks in the simulation
            dayon = ceil(Int64,7*rand(rng)) # which day do they start work
            lmul!(0,pap) # reinitalise pap
            # stay at work for toff days
            for t_it = 0:toff

                # Get day number of week.
                # uses mod1, rather than mod, so mod1(7,7) returns 7 rather than 0
                day_of_week = mod1(dayon+t_it,7)

                # Give value to pap
                pap[day_of_week] = 1
            end

            # iterate over each week
            for reps_itr = 1:num_reps
                # and put in the days at work in pap
                for day_itr = 1:7
                    atwork[node_itr,(reps_itr-1)*7+day_itr] = pap[day_itr]
                end
            end

            # and put in the last bit, if there aren't an exact number of repetitions
            # and put in the days at work in pap
            for day_itr = 1:7
                if (num_reps*7+day_itr)<=endtime
                    atwork[node_itr,num_reps*7+day_itr] = pap[day_itr]
                end
            end
        end
    elseif sameday==2 # workdays randomly placed throughout the week (but repeat each week)
        # iterate over each node
        pap = zeros(Int64,7)
        for node_itr=1:cmax
            randperm!(rng,pap)  # Get permutation of 1:7
            num_reps = floor(Int64,endtime/7) # how many weeks in the simulation
            # iterate over each week
            for reps_itr = 1:num_reps
                # and put in the first toff days in the random permuation, pap
                for pap_itr = 1:(toff+1)
                    atwork[node_itr,(reps_itr-1)*7+pap[pap_itr]] = 1
                end
            end
            # and put in the last bit, if there aren't an exact number of repetitions
            # and put in the days at work in pap
            for pap_itr = 1:(toff+1)
                if (num_reps*7+pap[pap_itr])<=endtime
                    atwork[node_itr,num_reps*7+pap[pap_itr]] = 1
                end
            end
        end
    elseif sameday==3 # do ton weeks on followed by toff weeks off, all simultaneous
        # iterate over each node
        for node_itr=1:cmax
            num_reps = floor(Int64,endtime/7) # how many weeks in the simulation
            # iterate over each week
            for reps_itr = 1:num_reps
                if mod1(reps_itr,ton+toff) <= ton # using mod1, as mod1(1,1) = 1, whereas mod(1,1) = 0
                    for days=1:5 # only work weekdays
                        atwork[node_itr,(reps_itr-1)*7+days] = 1
                    end
                end
            end
        end
    elseif sameday==4 # do ton weeks on followed by toff weeks off, starting at a random week each
        # initialise pap
        pap = zeros(ton+toff)
        # iterate over each node
        for node_itr=1:cmax
            # which week will they start
            weekon = ceil(Int64,(ton+toff)*rand(rng))
            lmul!(0,pap) # reinitalise pap
            # stay at work for ton weeks
            for t_it = 1:ton
                pap[mod1(weekon+(t_it-1),ton+toff)] = 1
            end
            num_reps = floor(Int64,endtime/7) # how many weeks in the simulation
            # iterate over each week
            for reps_itr = 1:num_reps
                if pap[mod1(reps_itr,ton+toff)] == 1
                    for days=1:5 # only work weekdays
                        atwork[node_itr,(reps_itr-1)*7+days] = 1
                    end
                end
            end
        end
    end

    return nothing
end

"""
    set_infection_related_times!(args)

Randomly initialise infection related values for each node
(time to symptoms, length of latent period, delay to adherence, adhere yes/no)

Inputs: `rng` - random number generator,
        `... parameter structures ...` \n
Outputs: None \n
Location: initialisation_fns.jl
"""
function set_infection_related_times!(rng::MersenneTwister,
                                        states::NodeStates,
                                        configuration_options::ConfigurationOptions,
                                        infection_params::InfectionParameters)

    @unpack cmax, isolation = configuration_options
    @unpack adherence, d_incub, inftime = infection_params
    @unpack adherence_hierarchy = states

    states.time_to_symps .= ceil.(rand(rng,d_incub,cmax)) # time to symptoms
    # (for asymptomatics, the same from a silent start of "symptoms")

    csum_delay_adherence = cumsum(infection_params.delay_adherence_pmf)

    num_adhering = adherence*cmax

    # iterate over nodes to set lattime and hh_isolation
    for node_itr = 1:cmax
        # lattime is the time from infection to infectiousness
        if states.time_to_symps[node_itr]-inftime<1
            states.lattime[node_itr] = 1  # Infectiousness can begin the day after becoming infected
        else
            states.lattime[node_itr] = states.time_to_symps[node_itr]-inftime
        end

        if adherence_hierarchy[node_itr] <= num_adhering # those who adhere will isolate when they get symptoms
            states.hh_isolation[node_itr] = 1 # adherence to household isolation = 1 if adherent, 0 if not.
        end

        # Draw random number
        # Set delay in adherence/symptoms becoming known to household
        # Find interval random number resides in
        states.delay_adherence[node_itr] = draw_sample_from_pmf(csum_delay_adherence,
                                                                rng;
                                                                idx_offset = 1)

    end
end

"""
    initialise_returned_to_work!(args)

Initialise 'returned_to_work' status of workers, according to sector-specific proportion of workers that are present at the workplace ('workpercent').

Proportion of workers at the workplace is applied at a sector level, not workplace level.
Differs to 'atwork', which tracks a workers weekly work / non-workday pattern.

Inputs: `worker_nodes` - array of WorkerParameters structures, containing returned-to-work status of each node,
        `... parameter structures ...`,
        `rng` - random number generator \n
Outputs: None \n
Location: initialisation_fns.jl
"""
function initialise_returned_to_work!(worker_nodes::Array{WorkerParameters,1},
                                        configuration_options::ConfigurationOptions,
                                        workplace_generation_params::WorkplaceGenerationParameters,
                                        network_params::NetworkParameters,
                                        rng::MersenneTwister)

    @unpack cmax, workertypes = configuration_options
    @unpack workpercent = workplace_generation_params
    @unpack workplace_sizes, nodes_by_sector = network_params

    for sector_itr = 1:workertypes
        n_workers = sum(workplace_sizes[sector_itr])

        returned_to_work_hierarchy = shuffle(rng, 1:n_workers)
        workplace_generation_params.returned_to_work_hierarchy[sector_itr] = returned_to_work_hierarchy

        num_returning = ceil(Int64, workpercent[sector_itr]*n_workers)

        for worker_itr = 1:n_workers
            current_node = nodes_by_sector[sector_itr][worker_itr]

            if returned_to_work_hierarchy[worker_itr] <= num_returning
                worker_nodes[current_node].returned_to_work = 1
            else
                worker_nodes[current_node].returned_to_work = 0
            end
        end
    end

    return nothing
end

"""
    assign_transrisk_all_settings!(args)

Randomly generates transmission risk for each node in each setting, according to specified transrisk functions.

Default transrisk functions are:
- assign_household_transmit_household_size! (household)
- assign_workplace_static_transmit! (workplace static)
- assign_workplace_dynamic_transmit! (workplace dynamic)
- assign_social_transmit! (social)
- assign_random_transmit! (random)

Inputs: `rng` - random number generator,
        `... parameter structures ...` \n
Outputs: None \n
Location: initialisation_fns.jl
"""
function assign_transrisk_all_settings!(rng::MersenneTwister,
                                        network_params::NetworkParameters,
                                        configuration_options::ConfigurationOptions,
                                        infection_params::InfectionParameters,
                                        contacts::ContactStructure)

    @unpack transrisk_household_group_mean, transrisk_household_group_sd,
            transrisk_static_work_mean, transrisk_static_work_sd,
            transrisk_dynamic_work_mean, transrisk_dynamic_work_sd,
            transrisk_social_mean, transrisk_social_sd,
            transrisk_random_mean, transrisk_random_sd,
            assign_household_transrisk_fn, assign_workplace_static_transrisk_fn,
            assign_workplace_dynamic_transrisk_fn, assign_social_transrisk_fn,
            assign_random_transrisk_fn = infection_params

    @unpack household_contacts_per_node = contacts


    # Relevant functions are listed in "include_files_network_model/additional_fns.jl"
    assign_household_transrisk_fn(rng,
                                network_params,
                                configuration_options,
                                household_contacts_per_node,
                                transrisk_household_group_mean,
                                transrisk_household_group_sd)

    assign_workplace_static_transrisk_fn(rng,
                                        network_params,
                                        configuration_options,
                                        transrisk_static_work_mean,
                                        transrisk_static_work_sd)

    assign_workplace_dynamic_transrisk_fn(rng,
                                        network_params,
                                        configuration_options,
                                        transrisk_dynamic_work_mean,
                                        transrisk_dynamic_work_sd)

    assign_social_transrisk_fn(rng,
                                network_params,
                                configuration_options,
                                transrisk_social_mean,
                                transrisk_social_sd)

    assign_random_transrisk_fn(rng,
                                network_params,
                                configuration_options,
                                transrisk_random_mean,
                                transrisk_random_sd)
end

"""
    reinitialise_workplace_params!(args)

Reinitialise workplace open and COVID-secure flags to configuration defaults.
Only called if workplace parameters have been changed via intervention.

Inputs: `workplace_info` - array{array} of WorkplaceParameters structures,
        `sector_open` - array indicating whether each sector is open,
        `CS_active_flag` - boolean flagging if COVID-secure workplace measures are in place (applies to all workplaces) \n
Outputs: None \n
Location: initialisation_fns.jl
"""
function reinitialise_workplace_params!(workplace_info::Array{Array{WorkplaceParameters,1},1},
                                        sector_open::Array{Bool,1},
                                        CS_active_flag::Bool)

    # Get number of sectors in use
    n_sectors = length(workplace_info)

    # Iterate over each sector.
    # Within each sector, iterate over each workplace
    # Update field values
    for sector_itr = 1:n_sectors
        n_workplaces = length(workplace_info[sector_itr])
        for workplace_itr = 1:n_workplaces
            workplace_info[sector_itr][workplace_itr].covid_secure = CS_active_flag
            workplace_info[sector_itr][workplace_itr].workplace_open = sector_open[sector_itr]
        end
    end

    return nothing
end

"""
    draw_sample_from_pmf(args)

Randomly draw sample from specified discrete probability distribution.
Used to randomly set infection, adherence and testing waiting times.

Inputs: `csum_pmf` - specified discrete CDF to sample from,
        `rng` - random number generator,
        `idx_offset` - links bin index to desired quantity value \n
Outputs: `val_to_update` - sampled quantity value \n
Location: initialisation_fns.jl
"""
function draw_sample_from_pmf(csum_pmf::Array{Float64,1},
                                rng::MersenneTwister;
                                idx_offset::Int64 = 0)

    # Get number of elements in the pmf
    n_bins = length(csum_pmf)

    # Initialise output value
    val_to_update = 0

    # Draw random number
    # Set delay in adherence/symptoms becoming known to household
    # Find interval random number resides in
    r = rand(rng)
    allocated_flag = false # Intialise allocation flag. Switch to true when allocation done.
    bin_idx = 1   # Current interval being checked
    while (allocated_flag == false)
        if r <= csum_pmf[bin_idx]
            # Assign selected value
            # Subtract idx_offset
            val_to_update = bin_idx - idx_offset

            # Update allocation flag
            allocated_flag = true
        else
            # r does not reside in this interval. Update bin index.
            bin_idx += 1

            # Error check, if not assigned value after checked final bin value
            if bin_idx > n_bins
                error("bin_idx is now $bin_idx. The pmf only has $n_bins bins. Terminating programme.")
            end
        end
    end

    return val_to_update::Int64
end
