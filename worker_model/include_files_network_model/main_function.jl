"""
File containing main function that generates a network according configuration options and runs
multiple simulations of an outbreak across it.
"""

"""
    worker_pattern_network_run(args)

Generates a network for a single configuration and runs countfinal number of replicates of each intervention set,

Inputs: `config_idx` - index of current configuration,
        `... parameter structures ...`,
        `all_intervention_sets` - array{array} of InterventionVariables structures, defining sets of timed interventions,
        `intervention_fns` - array of triggered intervention functions \n
Outputs: `output` - SimulationOutputs structure \n
Location: main_function.jl
"""
function worker_pattern_network_run(config_idx::Int64,
                                        configuration_variables::ConfigurationVariables,
                                        configuration_options::ConfigurationOptions,
                                        network_params::NetworkParameters,
                                        contact_tracing_params::ContactTracingParameters,
                                        infection_params::InfectionParameters,
                                        workplace_generation_params::WorkplaceGenerationParameters,
                                        workplace_closure_params::WorkplaceClosureParameters,
                                        all_intervention_sets::Array{Array{InterventionVariables,1},1};
                                        intervention_fns::Array{Function} = Array{Function}(undef))

    """
    Unpack configuration options
    """
    @unpack rng_seed, cmax, endtime, countfinal, workertypes = configuration_options

    """
    Set the random number generator
    """
    rng = MersenneTwister(rng_seed*config_idx)

    """
    Initialise intervention variables
    """
    # Check if any interventions were specified
    if isassigned(intervention_fns)
        # Number of intervention sets provided is number of rows of intervention_fns
        n_intervention_fns = size(intervention_fns,1)

        triggered_interventions_active = true
    end

    if isassigned(all_intervention_sets)
        n_intervention_sets = length(all_intervention_sets)
        timed_interventions_active = true
    else
        n_intervention_sets = 1
        timed_interventions_active = false
    end

    preintervention_settings = InterventionVariables()
    settings_changed_by_intervention_set = Symbol[]
    contact_component_store = PreinterventionComponentStore[]

    """
    Initialise workplaces and workers
    """
    @time generate_workplaces_and_allocate_workers!(configuration_options,
                                                    workplace_generation_params,
                                                    network_params,
                                                    rng)

    """
    Generate contacts within workplaces & households
    """
    @time contacts::ContactStructure,
        n_households::Int64 = generate_contacts(configuration_options,
                                                network_params,
                                                rng)

    """
    Generate social contacts (workdays and non-workdays)
    """
    contacts.social_meeting_group_sizes = Array{Array{Int64,1},1}(undef, n_intervention_sets+1)

    generate_social_contacts_each_day!(rng,
                                        1,
                                        network_params,
                                        contacts,
                                        configuration_options,
                                        0)

    """
    Generate dynamic contacts
    """
    # Set up dynamic worker contacts, in network_generation_fns.jl

    if network_params.network_generation_method == "ER"
        contacts.dynamic_worker_contacts = generate_dynamic_worker_contacts_ER(rng,
                                                                    cmax,
                                                                    endtime,
                                                                    network_params.worker_nodes,
                                                                    network_params.dynamic_conts_mean,
                                                                    network_params.dynamic_conts_sd)

    elseif network_params.network_generation_method == "configuration"
        contacts.dynamic_worker_contacts = generate_dynamic_worker_contacts_configuration(rng,
                                                                cmax,
                                                                endtime,
                                                                network_params.worker_nodes,
                                                                network_params.workplace_dynamic_degree_distribution,
                                                                network_params.max_contacts_work_dynamic)
    else
        println("Invalid network generation method!")
    end

    """
    Generate random dynamic contacts
    """
    contacts.dynamic_random_contacts = generate_random_contacts(rng,
                                                        cmax,
                                                        endtime,
                                                        network_params.prob_random_contact)

    """
    Initialise storage arrays
    """
    output = SimulationOutputs(endtime=endtime,countfinal=countfinal,cmax=cmax,
                            n_intervention_sets = max(length(all_intervention_sets), 1))

    """
    Initialise node states (reinitialised each replicate)
    """
    states = NodeStates(cmax=cmax, endtime=endtime)

    """
    Initialise contact tracing related variables
    """
    if configuration_options.contact_tracing_active == true
        CT_vars = ContactTracingVariables(cmax=cmax)
    else
        CT_vars = ContactTracingVariables()
    end

    """
    If required, initialise workplace closure variables
    """
    if configuration_options.workplace_closure_active == true
        initialise_workplace_closure_params!(workplace_closure_params,
                                         configuration_options,
                                         network_params)
    end

    """
    Run different intervention sets
    """
    for intervention_set_itr = 1:n_intervention_sets

        """
        Select timed intervention set and check timings
        """
        if timed_interventions_active
            println("Intervention set: $intervention_set_itr")

            intervention_set = all_intervention_sets[intervention_set_itr]

            # Find number and times of specified interventions
            n_interventions = length(intervention_set)
            if n_interventions > 0
                intervention_times = [intervention_set[i].start_time[1] for i=1:n_interventions]
            end
            println("Interventions timed for days: $intervention_times")

            if n_interventions > 1
                time_diff = intervention_times[2:end] .- intervention_times[1:end-1]
                if any(x->x<=0, time_diff)
                    error("Interventions not in chronological order")
                end
            end

            reset_times = Int64[]
            for i = 1:n_interventions
                if isassigned(intervention_set[i].reset_time)
                    push!(reset_times, intervention_set[i].reset_time[1])
                end
            end
            println("Interventions reset on days: $reset_times")

        else
            println("No timed interventions")

            n_interventions = 0
        end


        """
        Run replicates
        """
        # Perform countfinal number of replicates
        for count=1:countfinal

            """
            Initialisation phase
            """

            initialise_settings_and_components_for_replicate!(count,
                                                                intervention_set_itr,
                                                                rng,
                                                                output,
                                                               configuration_options,
                                                               infection_params,
                                                               contact_tracing_params,
                                                               CT_vars,
                                                               contacts,
                                                               network_params,
                                                               workplace_generation_params,
                                                               workplace_closure_params,
                                                               states,
                                                               preintervention_settings,
                                                               contact_component_store,
                                                               settings_changed_by_intervention_set,
                                                               timed_interventions_active)

           """
           Reset intervention trackers between intervention sets
           """
           if (count == 1) && (intervention_set_itr > 1)
               preintervention_settings = InterventionVariables()
               settings_changed_by_intervention_set = Symbol[]
               #contact_component_store = PreinterventionComponentStore[]
           end

            """
            Generate random number arrays
            """
            r_symp_test = rand(rng, cmax, endtime)
            r_backwards_CT = rand(rng, cmax, endtime)

            """
            Check for errors
            """
            check_for_errors(configuration_options)

            """
            Run single replicate
            """
            for time=1:endtime

                # Initial timepoint is for initial conditions
                # Set row to accessed in output arrays for this timestep
                output_time_idx = time + 1

                # if (time == 1) ||  (time == 20) ||  (time == 21) || (time == 50) ||  (time == 51)
                #     println("Time: $time, generation = $(network_params.network_generation_method_dynamic_social)")
                # end

                """
                Implement or reset intervention
                """
                if n_interventions > 0
                    if time ∈ intervention_times
                        intervention_ids = findall(intervention_times.==time)

                        for intervention_itr = 1:length(intervention_ids)
                            current_intervention = intervention_set[intervention_ids[intervention_itr]]

                            affect_intervention!(time,
                                                current_intervention,
                                                intervention_set_itr,
                                                configuration_options,
                                                infection_params,
                                                contact_tracing_params,
                                                CT_vars,
                                                contacts,
                                                network_params,
                                                workplace_generation_params,
                                                states,
                                                workplace_closure_params,
                                                preintervention_settings,
                                                settings_changed_by_intervention_set,
                                                contact_component_store,
                                                rng)
                        end
                        check_for_errors(configuration_options)
                    end
                    if time ∈ reset_times
                        affect_intervention!(time,
                                            preintervention_settings,
                                            intervention_set_itr,
                                            configuration_options,
                                            infection_params,
                                            contact_tracing_params,
                                            CT_vars,
                                            contacts,
                                            network_params,
                                            workplace_generation_params,
                                            states,
                                            workplace_closure_params,
                                            preintervention_settings,
                                            settings_changed_by_intervention_set,
                                            contact_component_store,
                                            rng)

                        check_for_errors(configuration_options)
                    end
                end
                """
                Reinitialise variables at start of timestep
                """
               # Reinitialise timestep specific values
               lmul!(0,states.rep_inf_this_timestep)

               # reinitalise the current workplace_memory slot
               if configuration_options.workplace_closure_active == true
                   workplace_closure_params.WP_memory_slot = mod1(time, workplace_closure_params.workplace_CT_memory)
                   for worktypeID = 1:workertypes
                       # Iterate over each workplace for current sector type
                       n_workplaces = workplace_closure_params.num_workplaces[worktypeID]
                       for workplace_itr = 1:n_workplaces
                           workplace_closure_params.workplace_memory[worktypeID][workplace_itr,workplace_closure_params.WP_memory_slot] = 0
                       end
                   end
               end

                """
                Assign outputs
                """
                # Assign counts in each disease state to array
                output.numlat[output_time_idx,count,intervention_set_itr] = output.numlat[output_time_idx-1,count,intervention_set_itr]
                output.numinf[output_time_idx,count,intervention_set_itr] = output.numinf[output_time_idx-1,count,intervention_set_itr]
                output.numrep[output_time_idx,count,intervention_set_itr] = output.numrep[output_time_idx-1,count,intervention_set_itr]

                """
                Increment counters
                """
                # Increment counters if node is currently in that state.
                increment_counters!(states) # in additional_fns.jl

                """
                Increment infection process
                """
                increment_infection_process!(output_time_idx,
                                                count,
                                                intervention_set_itr,
                                                r_symp_test,
                                                configuration_options,
                                                states,
                                                output,
                                                network_params,
                                                infection_params,
                                                contact_tracing_params,
                                                CT_vars,
                                                contacts,
                                                workplace_closure_params)

                """
                Get daily isolation & atwork status of nodes
                """
                get_daily_isol_atwork_status!(time,
                                            output_time_idx,
                                            configuration_options,
                                            states,
                                            contacts,
                                            network_params)

                """
                Transmit infections

                Structure:
                - Household
                - At work
                    -- Social contacts
                    -- Work contacts (with checks based on covid-secure status)
                    -- Dynamic contacts
                - Not at work
                """
                transmit_infections!(rng,
                                        time,
                                        count,
                                        intervention_set_itr,
                                        configuration_options,
                                        states,
                                        contacts,
                                        network_params,
                                        infection_params,
                                        output)

                """
                Perform contact tracing
                """
                # If in use, enact contact tracing from index cases reported today
                if (configuration_options.contact_tracing_active == true) &&
                    (configuration_options.isolation > 0)
                    perform_contact_tracing!(rng,
                                            time,
                                            output_time_idx,
                                            count,
                                            intervention_set_itr,
                                            r_backwards_CT,
                                            configuration_options,
                                            states,
                                            contacts,
                                            network_params,
                                            infection_params,
                                            contact_tracing_params,
                                            CT_vars,
                                            output)
                end

                """
                Assign prevalence & isolation outputs
                """
                assign_prev_isol_outputs!(output_time_idx,
                                            count,
                                            intervention_set_itr,
                                            configuration_options,
                                            states,
                                            output)

                """
                Reactive workplace closure check
                """
                # close workplaces with too many infections
                if configuration_options.workplace_closure_active==true
                    reactive_workplace_closures!(configuration_options,
                                                        network_params,
                                                        workplace_closure_params)
                end

                """
                Check and run triggered interventions
                """
                # Check if any interventions are triggered
                # Update statuses as needed
                if isassigned(intervention_fns) # Check if any intervetion were specified
                    check_triggered_interventions!(time,
                                                    network_params,
                                                    states,
                                                    output)
                end

            end

            # Find how many nodes were infected by the initial nodes
            initial_nodes = findall(states.infected_by.==-1)
            sum_infections = 0
            output.num_init_infected[count] = zeros(Int64,length(initial_nodes),max(n_intervention_sets,1)) # Initialise output array
            for initial_node_it = 1:length(initial_nodes)
                output.num_init_infected[count][initial_node_it,intervention_set_itr] = output.num_infected[initial_nodes[initial_node_it],count,intervention_set_itr]
                sum_infections+=output.num_init_infected[count][initial_node_it,intervention_set_itr]
            end

            # find mean generation time
            if sum_infections>0
                output.mean_init_generation_time[count,intervention_set_itr] = output.mean_init_generation_time[count,intervention_set_itr]/sum_infections
            end

            # divide number of infections by number of infectors to find Rt
            for time=1:(endtime+1)
                # divide by the number of nodes that were infected (entered the latent state)
                # at time
                if time == 1
                    output.Rt[time,count,intervention_set_itr] = output.Rt[time,count,intervention_set_itr] / output.numlat[time,count,intervention_set_itr]
                else
                    output.Rt[time,count,intervention_set_itr] = output.Rt[time,count,intervention_set_itr] / (output.numlat[time,count,intervention_set_itr]-output.numlat[time-1,count,intervention_set_itr])
                end
            end

            # Print to screen info on run just completed
            println("Run $count, intervention $intervention_set_itr, configuration $config_idx complete.")
        end

    end

    # Calculate summary of social meeting group sizes, if relevant
    if any([isassigned(contacts.social_meeting_group_sizes,ii) for ii=1:(n_intervention_sets+1)]) == true
        for intervention_set_itr = 1:(n_intervention_sets+1)
            if isassigned(contacts.social_meeting_group_sizes,intervention_set_itr)
                output.social_meeting_groups[1, intervention_set_itr] = mean(contacts.social_meeting_group_sizes[intervention_set_itr])
                output.social_meeting_groups[2:8, intervention_set_itr] = quantile(contacts.social_meeting_group_sizes[intervention_set_itr], [0.025,0.05,0.25,0.5,0.75,0.95,0.975])
            end
        end
    end

    # Compute variance in number of infected per node
    # var_num_infected = zeros(countfinal)
    # for count=1:countfinal
    #     var_num_infected[count] = var(output.num_infected[:,count,intervention_set_itr])
    # end

    # Specify what is output from the function
    return output
end
