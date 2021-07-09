"""
File containing functions that are related to contact tracing
"""

"""
    initialise_CT_vars!(args)

Initialise objects in ContactTracingVariables structure according to specified number of nodes.

Only called when contact tracing is switched on via intervention.

Inputs: `cmax` - number of nodes,
        `CT_vars` - ContactTracingVariables structure to store necessary values \n
Outputs: None \n
Location: contact_tracing_fns.jl
"""
function initialise_CT_vars!(cmax::Int64, CT_vars::ContactTracingVariables)

   CT_vars.cmax = cmax

   # The number of days prior to symptoms that each node remembers
   CT_vars.relevant_prev_days_for_CT::Array{Int64,1} = zeros(Int64,cmax)

   # Vector of vectors for storing IDs of those to be contacted in CT
   CT_vars.Inds_to_be_contacted::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}(undef,cmax)

   # # vector tracking symptomatic cases (positive confirmed or untested)
   # Symp_cases_per_household_pos_or_unknown::Array{Int64,1} = zeros(Int64,n_households)

   # Variables for waiting for test results, set at -1 until activated
   CT_vars.Time_to_test_result::Array{Int64,1} = -1*ones(Int64,cmax)

   # Boolean vector to store whether a false negative test result would be returned
   # Populated in each replicate
   CT_vars.Test_result_false_negative::Array{Bool,1} = Array{Bool,1}(undef,cmax)

   # Boolean vector to store if individual will engage with CT or not
   # Populated in each replicate
   CT_vars.Engage_with_CT::Array{Bool,1} = Array{Bool,1}(undef,cmax)

   # Array to keep track of whether an infected recalls their infector
   CT_vars.Recall_infector::Array{Int64,1} = zeros(Int64,cmax)

   # Delay before test result is returned
   CT_vars.CT_delay_until_test_result::Array{Int64,1} = zeros(Int64,cmax)

end

"""
    reinitialise_CT_vars!(args)

Reinitialise objects in ContactTracingVariables structure.

Delay until test result for each node are randomly drawn (called every replicate).
Engagement with test and trace depends on adherence hierarchy.

Inputs: `CT_vars` - ContactTracingVariables structure to store necessary values,
        `cmax` - number of nodes,
        `rng` - random number generator,
        `... parameter structures ...` \n
Outputs: None \n
Location: contact_tracing_fns.jl
"""
function reinitialise_CT_vars!(CT_vars::ContactTracingVariables,
                                cmax::Int64,
                                rng::MersenneTwister,
                                CT_parameters::ContactTracingParameters,
                                states::NodeStates,
                                infection_params::InfectionParameters)

    @unpack CT_days_before_symptom_included, CT_engagement,
                CT_delay_until_test_result_pmf = CT_parameters
    @unpack adherence_hierarchy, delay_adherence = states
    @unpack adherence = infection_params


    # Variables for waiting for test results
    lmul!(0,CT_vars.Time_to_test_result)
    CT_vars.Time_to_test_result .-= 1 # Reset so all values are -1

    # Repopulate Boolean vector stating whether a false negative test result would be returned
    # and the number of days relevant for contact tracing
    lmul!(0,CT_vars.relevant_prev_days_for_CT)

    num_engaging = CT_engagement*adherence*cmax

    for ii = 1:cmax

      # For each worker, initialise CT_vars.Test_result_false_negative as false
      CT_vars.Test_result_false_negative[ii] = false

      # For each worker, check if they engage with contact tracing
      if adherence_hierarchy[ii] <= num_engaging # engage with contact tracing
          CT_vars.Engage_with_CT[ii] = true
      else # do not engage with contact tracing
          CT_vars.Engage_with_CT[ii] = false
      end

      # Get amount of days to be looked back over
      # Have upper bound of 7 days post symp
      # if we put in reporting delay, needs to be above this
      CT_vars.relevant_prev_days_for_CT[ii] = min(CT_days_before_symptom_included + delay_adherence[ii],
                                          CT_days_before_symptom_included + 7)
    end

    csum_test_result_delay = cumsum(CT_delay_until_test_result_pmf)

    # Repopulate time until test result received for each individual
    lmul!(0,CT_vars.CT_delay_until_test_result)
    for node_itr = 1:cmax
                      CT_vars.CT_delay_until_test_result[node_itr] = draw_sample_from_pmf(csum_test_result_delay,
                                                                                rng;
                                                                                idx_offset = 1)
    end

    # Set up vector of vectors for storing IDs of those to be contacted in CT
    CT_vars.Inds_to_be_contacted = Array{Array{Int64,1},1}(undef,cmax)

    # Initialise array to keep track of whether an infected recalls their infector
    lmul!(0,CT_vars.Recall_infector)
end

"""
    perform_contact_tracing!(args)

For nodes becoming symptomatic, check whether they engage with contact tracing and, if so, contact trace previous contacts.

A delay between becoming symptomatic and engaging with contact tracing may also apply.

Inputs: `rng` - random number generator,
        `time` - current time in simulation,
        `output_time_idx` - current index in output array (corresponds to time + 1),
        `count` - number of current replicate,
        `intervention_set_itr` - number of current intervention set,
        `r_backwards_CT` - random number array to determine if backwards contact tracing is used,
        `... parameter structures ...` \n
Outputs: None \n
Location: contact_tracing_fns.jl
"""
function perform_contact_tracing!(rng::MersenneTwister,
                                time::Int64,
                                output_time_idx::Int64,
                                count::Int64,
                                intervention_set_itr::Int64,
                                r_backwards_CT::Array{Float64,2},
                                configuration_options::ConfigurationOptions,
                                states::NodeStates,
                                contacts::ContactStructure,
                                network_params::NetworkParameters,
                                infection_params::InfectionParameters,
                                contact_tracing_params::ContactTracingParameters,
                                CT_vars::ContactTracingVariables,
                                output::SimulationOutputs)

    @unpack cmax, isolation, endtime, perform_CT_from_infector = configuration_options
    @unpack CT_engagement, prob_backwards_CT, CT_caused_isol_limit = contact_tracing_params

    # Store contacts made during day.
    # For those reporting symptoms, start delay to test result (if needed)
    for node_itr = 1:cmax

        # Check engagement with CT
        if CT_vars.Engage_with_CT[node_itr] == true

            # Increment time to test result if currently waiting for that to occur
            if CT_vars.Time_to_test_result[node_itr] >= 0
                CT_vars.Time_to_test_result[node_itr] += 1
            end

            # For current worker, check if would be leaving pre-symptomatic phase
            # and displaying symptoms
            # If so, and will not return a false negative test result, gather traceable contacts
            if (states.rep_inf_this_timestep[node_itr] == 1)

                # Increment test counter
                output.tests_performed[output_time_idx,count,intervention_set_itr] += 1

                # Initialise CT_vars.Time_to_test_result value
                CT_vars.Time_to_test_result[node_itr] = 0

                # Check if the worker is destined to return a false negative result
                # If false negative to be returned, do not need to work out who traceable contacts are
                # Otherwise, gather traceable contacts
                CT_vars.Inds_to_be_contacted[node_itr] = Int64[] # Initialise vector to store contacts
                if (CT_vars.Test_result_false_negative[node_itr] == false)

                    trace_node!(node_itr,
                                time,
                                CT_vars,
                                contacts,
                                contact_tracing_params,
                                network_params,
                                states,
                                configuration_options,
                                rng)

                    # if we are doing "backward contact tracing"
                    # some small chance the infector is included in this
                    # don't try to backwards trace the initial infections
                    if (r_backwards_CT[node_itr,time] < prob_backwards_CT) && (states.infected_by[node_itr] != -1)
                        append!(CT_vars.Inds_to_be_contacted[node_itr],states.infected_by[node_itr])
                            CT_vars.Recall_infector[node_itr] = 1
                    end
                end

                # Remove duplicates in CT_vars.Inds_to_be_contacted[node_itr]
                unique!(CT_vars.Inds_to_be_contacted[node_itr])
            end

            # Check if delay to test result reached
            if CT_vars.Time_to_test_result[node_itr] >= CT_vars.CT_delay_until_test_result[node_itr]

                # Reset the CT_vars.Time_to_test_result counter
                CT_vars.Time_to_test_result[node_itr] = -1

                # If delay time passed, check if test returned a false negative.
                if CT_vars.Test_result_false_negative[node_itr] == true
                    # Increment false negative counter
                    output.test_outcomes[output_time_idx,count,intervention_set_itr,2] += 1
                else
                    # Increment true positive counter
                    output.test_outcomes[output_time_idx,count,intervention_set_itr,1] += 1

                    # If test is positive, contacts told to isolate
                    # Get number of recallable contacts
                    n_recallable_contacts = length(CT_vars.Inds_to_be_contacted[node_itr])
                    output.num_CT[output_time_idx,count,intervention_set_itr] += n_recallable_contacts

                    # Do contacts isolate or not?
                    # Isolation based on adherence to household isolation of that individual
                    for recallable_contact_idx = 1:n_recallable_contacts
                        recallable_contact_ID = CT_vars.Inds_to_be_contacted[node_itr][recallable_contact_idx]

                        # Check if individual will adhere to guidance
                        # If so, they self-isolate
                        if states.hh_isolation[recallable_contact_ID] == 1
                            # Populate contact tracing isolation time tracker array
                            # Time needed to spend in isolation reduced by test result delay time for index case
                            time_in_CT_isol = max(0,CT_caused_isol_limit - CT_vars.CT_delay_until_test_result[node_itr])
                            for time_itr = 1:time_in_CT_isol
                                array_time_idx = time + time_itr
                                if array_time_idx <= (endtime + 1)
                                    states.CT_isolation_array[recallable_contact_ID,array_time_idx] = 1
                                end
                                    # Note, in isolation_array, col 1 corresponds to time 0,
                                    #                           col 2 corresponds to time 1 etc
                                    # So array_time_idx = time + 1, populates array
                                    # in column corresponding to day value of "time"
                            end
                        end
                    end

                    # Perform forwards CT from infector, if infector has been identified
                    # and such a policy is active
                    if perform_CT_from_infector == true
                        if CT_vars.Recall_infector[node_itr] == 1

                            forwardCT_from_infector!(states.infected_by[node_itr],
                                                        CT_vars,
                                                        contacts,
                                                        contact_tracing_params,
                                                        network_params,
                                                        states,
                                                        configuration_options,
                                                        CT_vars.Engage_with_CT,
                                                        states.atwork,
                                                        time,
                                                        count,
                                                        rng)
                        end
                    end
                end
            end
        end
    end
end

"""
    recallable_dynamic_contacts(args)

For a single specified node, collect expected dynamic contacts for a specified day and check if they occurred and are recallable.

Inputs: `worker_ID` - ID of node being contact traced,
        `dynamic_contact_record` - IDs of dynamic contacts for all nodes on all days,
        `dynamic_contacts_recalled_propn` - array defining proportion of contacts recalled x days ago,
        `daily_record_inisol` - indicates isolation status of each node each day,
        `daily_record_atworkplace` - indicates at workplace status of each node each day
        `time_to_check` - time that contacts to be checked took place,
        `prev_day_val` - number of days prior to current time that time_to_check is,
        `rng` - random number generator \n
Outputs: `recallable_contacts` - array of node IDs that have been recalled as contacts \n
Location: contact_tracing_fns.jl
"""
function recallable_dynamic_contacts(worker_ID::Int64,
                                        dynamic_contact_record::Array{Array{Int64,1},2},
                                        dynamic_contacts_recalled_propn::Array{Float64,1},
                                        daily_record_inisol::Array{Int64,2},
                                        daily_record_atworkplace::Array{Int64,2},
                                        time_to_check::Int64,
                                        prev_day_val::Int64,
                                        rng::MersenneTwister)

    # Get proportion of dynamic contacts from prev_day_val days ago to be retained
    if prev_day_val > length(dynamic_contacts_recalled_propn)
        # No contacts can possible be retained.
        # Return an empty vector
        recallable_contacts = Int64[]
    else
        # Get proportion of contacts expected to be retained from input distribution
        threshold = dynamic_contacts_recalled_propn[prev_day_val]

        # Get IDs of those dynamic contacts on required day
        all_dynamic_contacts = dynamic_contact_record[time_to_check,worker_ID]

        # Check if any dynamic contacts were isolating.
        # If so, will not be a recallable contact
        n_possible_dynamic_contacts = length(all_dynamic_contacts)
        recallable_contact_check = zeros(Int64,n_possible_dynamic_contacts)
        for contact_itr = 1:n_possible_dynamic_contacts
            # If not isolating, then contact did occur
            contact_ID = all_dynamic_contacts[contact_itr]
            if daily_record_inisol[time_to_check,contact_ID]==false
                # Contact occurred
                # Check if contact will be remembered
                r = rand(rng)
                if r<threshold
                    recallable_contact_check[contact_itr] = 1
                end
            end
        end

        # Construct vector of IDs of contacts that did actually occur
        n_recallable_contacts = sum(recallable_contact_check)
        recallable_contacts = zeros(Int64,n_recallable_contacts) # Initialise vector to store IDs of recallable contacts
        recallable_contact_itr = 1 # Initialise counter for vector assignment index
        for contact_itr = 1:n_possible_dynamic_contacts
            if recallable_contact_check[contact_itr] == true
                contact_ID = all_dynamic_contacts[contact_itr]
                recallable_contacts[recallable_contact_itr] = contact_ID
                recallable_contact_itr += 1 # Increment vector assignment index
            end
        end
    end

    return recallable_contacts::Array{Int64,1}
end

"""
    get_worker_contacts(args)

For a single specified node, go over usual workday contacts and check if they actually occurred on a given day.

Inputs: `worker_ID` - ID of node being contact traced,
        `possible_worker_contacts` - IDs of usual workplace contacts, to be checked,
        `daily_record_inisol` - indicates isolation status of each node each day,
        `daily_record_atworkplace` - indicates at workplace status of each node each day
        `time_to_check` - time that contacts to be checked took place,
        `prev_day_val` - number of days prior to current time that time_to_check is,
        `rng` - random number generator,
        `network_params` - NetworkParameters structure,
        `other_workplace_flag` - boolean flagging if tracing contacts in same or different workplace \n
Outputs: `workplace_contacts` - array of node IDs that have been confirmed as contacts \n
Location: contact_tracing_fns.jl
"""
function get_worker_contacts(worker_ID::Int64,
                                possible_worker_contacts::Array{Int64,1},
                                daily_record_inisol::Array{Int64,2},
                                daily_record_atworkplace::Array{Int64,2},
                                time_to_check::Int64,
                                prev_day_val::Int64,
                                rng::MersenneTwister,
                                network_params::NetworkParameters;
                                other_workplace_flag::Bool = false)

    # Unpack parameters required
    if other_workplace_flag == true
        @unpack worker_nodes, workplace_info = network_params
    end

    # Check if any contacts were isolating and/or not at workplace that day
    # If so, will not be a contact
    n_possible_worker_contacts = length(possible_worker_contacts)
    contact_occur_check = zeros(Int64,n_possible_worker_contacts)
    for contact_itr = 1:n_possible_worker_contacts

        # If contact not isolating and at workplace, then contact can occur in same workplace.
        # If contact in another workplace, also need to check CS status
        contact_ID = possible_worker_contacts[contact_itr]
        if (daily_record_inisol[time_to_check,contact_ID]==false) &&
            (daily_record_atworkplace[time_to_check,contact_ID]==true)

            # Check if dealing with same workplace contact
            # or contact with other workplaces
            if other_workplace_flag == true
                # Contact in other workplace. Need to check CS status of that workplace

                # Get information on the contacts workplace
                contact_sector_ID = worker_nodes[contact_ID].sector_ID
                contact_workplace_ID = worker_nodes[contact_ID].workplace_ID
                contact_workplace_info = workplace_info[contact_sector_ID][contact_workplace_ID]

                # Get contacts workplace CS status
                # If not a CS setting, contact has taken place
                contact_workplace_CS_bool = contact_workplace_info.covid_secure
                if (contact_workplace_CS_bool==false)
                    contact_occur_check[contact_itr] = 1
                end
            else
                # Contact in same workplace
                contact_occur_check[contact_itr] = 1
            end
        end
    end

    # Construct vector of IDs of contacts that did actually occur
    n_occur_worker_contacts = sum(contact_occur_check)
    workplace_contacts = zeros(Int64,n_occur_worker_contacts) # Initialise vector to store IDs of recallable contacts
    contact_occur_itr = 1 # Initialise counter for vector assignment index
    for contact_itr = 1:n_possible_worker_contacts
        if contact_occur_check[contact_itr] == true
            contact_ID = possible_worker_contacts[contact_itr]
            workplace_contacts[contact_occur_itr] = contact_ID
            contact_occur_itr += 1 # Increment vector assignment index
        end
    end

    return workplace_contacts::Array{Int64,1}
end

"""
    forwardCT_from_infector!(args)

Perform forward contact-tracing from an identified infector, if infector reports symptoms and engages with contact-tracing.

Inputs: `infector_ID` - ID of infector of contact traced node,
        `... parameter structures ...`,
        `engage_with_CT` - array indicating whether each node engages with contact-tracing,
        `atwork` - work schedule for each node,
        `time` - current time in simulation,
        `count` - number of current replicate,
        `rng` - random number generator \n
Outputs: None \n
Location: contact_tracing_fns.jl
"""
function forwardCT_from_infector!(infector_ID::Int64,
                                CT_vars::ContactTracingVariables,
                                contacts::ContactStructure,
                                CT_parameters::ContactTracingParameters,
                                network_params::NetworkParameters,
                                states::NodeStates,
                                configuration_options::ConfigurationOptions,
                                engage_with_CT::Array{Bool,1},
                                atwork::Array{Int64,2},
                                time::Int64,
                                count::Int64,
                                rng::MersenneTwister)

@unpack dynamic_contacts_recalled_propn, social_contacts_recalled_propn, infector_engage_with_CT_prob = CT_parameters
@unpack Inds_to_be_contacted, Test_result_false_negative = CT_vars

    # Already contact traced? If not, do it now
    if (isassigned(Inds_to_be_contacted,infector_ID)) # Checks if infector_ID reported infection.
                                                      # If not, no forward contact tracing will be done from infector
        if !(length(Inds_to_be_contacted[infector_ID])>0) # Inds_to_be_contacted[infector_ID] being empty signifies
                                                         # infector reported symptoms, but returned false neg and/or
                                                         # did not engage in contact tracing
            # check if the infector will engage with CT
            if (((engage_with_CT[infector_ID] == false) && (rand(rng)<infector_engage_with_CT_prob))
                || (engage_with_CT[infector_ID] == true)) &&
                (Test_result_false_negative[infector_ID] == false)
                        # didn't engage before but does now or already willing to engage
                        # and then checks infector has not previously tested negative

                trace_node!(infector_ID,
                            time,
                            CT_vars,
                            contacts,
                            CT_parameters,
                            network_params,
                            states,
                            configuration_options,
                            rng)
            end
        end
    end
end

"""
    trace_node!(args)

For a single specified node, contact trace all contacts for a specified number of days in the past.

Static work contacts and dynamic (work or social) contacts are treated differently. Random contacts are not traceable.

Inputs: `node_itr` - ID of node being contact traced,
        `time` - current time in simulation,
        `... parameter structures ...`,
        `rng` - random number generator \n
Outputs: None \n
Location: contact_tracing_fns.jl
"""
function trace_node!(node_itr::Int64,
                        time::Int64,
                        CT_vars::ContactTracingVariables,
                        contacts::ContactStructure,
                        CT_parameters::ContactTracingParameters,
                        network_params::NetworkParameters,
                        states::NodeStates,
                        configuration_options::ConfigurationOptions,
                        rng::MersenneTwister)

@unpack worker_nodes, workplace_info = network_params
@unpack CS_active_flag = configuration_options

    # Preload workplace contacts and other worker contacts
    if CS_active_flag == true
        possible_same_workplace_contacts = contacts.work_contacts_same_workplace_CS[node_itr]
    else
        possible_same_workplace_contacts = contacts.work_contacts_same_workplace[node_itr]
    end
    possible_other_workplace_contacts = contacts.work_contacts_other_workplace[node_itr]

    # Preload workplace info
    workertype_ID = worker_nodes[node_itr].sector_ID
    workplace_ID = worker_nodes[node_itr].workplace_ID
    current_workplace_info = workplace_info[workertype_ID][workplace_ID]

    # Get contacts that will be contacted
    for time_itr = 1:CT_vars.relevant_prev_days_for_CT[node_itr]
        if time-time_itr > 0 # Can't look back before the simulation started

            # Get previous time being checked
            time_to_check = time-time_itr

            # If at workplace, get contacts made that day
            if (states.daily_record_atworkplace[time_to_check,node_itr]==1)

                # Note, function "get_worker_contacts"
                # resides in contact_tracing_fns.jl
                # Workplace contacts
                workplace_contacts = get_worker_contacts(node_itr,
                                                            possible_same_workplace_contacts,
                                                            states.daily_record_inisol,
                                                            states.daily_record_atworkplace,
                                                            time_to_check,
                                                            time_itr,
                                                            rng,
                                                            network_params)

                # If there are workplace contacts,
                # add to vector tracking traceable contacts
                if !isempty(workplace_contacts)
                    append!(CT_vars.Inds_to_be_contacted[node_itr],workplace_contacts)
                end
                # Workers in a CS settings will not make any non-workplace
                # worker contacts
                if (current_workplace_info.covid_secure == false)
                    # Worker contacts at other workplaces
                    other_workplace_flag_val = true
                    other_workplace_contacts = get_worker_contacts(node_itr,
                                                                possible_other_workplace_contacts,
                                                                states.daily_record_inisol,
                                                                states.daily_record_atworkplace,
                                                                time_to_check,
                                                                time_itr,
                                                                rng,
                                                                network_params,
                                                                other_workplace_flag=other_workplace_flag_val
                                                                )

                    # If there are workday contacts with
                    # inds at other workplaces,
                    # add to vector tracking traceable contacts
                    if !isempty(other_workplace_contacts)
                        append!(CT_vars.Inds_to_be_contacted[node_itr],other_workplace_contacts)
                    end
                end

                # Dynamic contacts as part of job role at workplace
                # Check at workplace on current timestep
                atwork_dynamic_recallable_contacts = recallable_dynamic_contacts(node_itr,
                                                                    contacts.dynamic_worker_contacts,
                                                                    CT_parameters.dynamic_contacts_recalled_propn,
                                                                    states.daily_record_inisol,
                                                                    states.daily_record_atworkplace,
                                                                    time_to_check,
                                                                    time_itr,
                                                                    rng) # in contact_tracing_fns.jl

                # If there are recallable dynamic contacts, add to vector tracking traceable contacts
                if !isempty(atwork_dynamic_recallable_contacts)
                    append!(CT_vars.Inds_to_be_contacted[node_itr],atwork_dynamic_recallable_contacts)
                end
            end

            # Social contacts check
            # Check node not in isolation, and has
            # potentially made social contacts
            if (states.daily_record_inisol[time_to_check,node_itr]==false) &&
                (contacts.social_contacts_per_node[node_itr] > 0)

                # Workday social contacts
                if (states.daily_record_atworkplace[time_to_check,node_itr]==true)
                    workday_social_recallable_contacts = recallable_dynamic_contacts(node_itr,
                                                                        contacts.workday_social_contacts_by_day,
                                                                        CT_parameters.social_contacts_recalled_propn,
                                                                        states.daily_record_inisol,
                                                                        states.daily_record_atworkplace,
                                                                        time_to_check,
                                                                        time_itr,
                                                                        rng) # in contact_tracing_fns.jl

                    if !isempty(workday_social_recallable_contacts)
                        append!(CT_vars.Inds_to_be_contacted[node_itr],workday_social_recallable_contacts)
                    end

                else
                    # Non-workday social contacts
                    nonworkday_social_recallable_contacts = recallable_dynamic_contacts(node_itr,
                                                                    contacts.nonworkday_social_contacts_by_day,
                                                                    CT_parameters.social_contacts_recalled_propn,
                                                                    states.daily_record_inisol,
                                                                    states.daily_record_atworkplace,
                                                                    time_to_check,
                                                                    time_itr,
                                                                    rng) # in contact_tracing_fns.jl

                    # If there are recallable dynamic contacts, add to vector tracking traceable contacts
                    if !isempty(nonworkday_social_recallable_contacts)
                        append!(CT_vars.Inds_to_be_contacted[node_itr],nonworkday_social_recallable_contacts)
                    end
                end
            end
        end
    end
end
