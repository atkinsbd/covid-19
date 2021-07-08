"""
File containing functions related to the transmission of infection and progression through disease states
"""

"""
    increment_counters!(states::NodeStates)

Increase value of time in infection state variables, called each timestep.

Inputs: `states` - NodeStates structure containing disease states of each node \n
Outputs: None \n
Location: infection_process_fns.jl
"""
function increment_counters!(states::NodeStates)

# Increments time counters
    cmax = length(states.timelat)
    for node_itr = 1:cmax
        if states.timelat[node_itr]>0
            states.timelat[node_itr] += 1
        end

        if states.timeinf[node_itr]>0
            states.timeinf[node_itr] += 1
        end

        if states.timesymp[node_itr]>0
            states.timesymp[node_itr] += 1
        end
    end
    return nothing
end

"""
    increment_infection_process!(args)

Move nodes through disease states according to time spent in current state, called each timestep.

Inputs: `output_time_idx` - current index in output array (corresponds to time + 1),
        `count` - number of current replicate,
        `intervention_set_itr` - number of current intervention set,
        `r_symp_test` - random number array to determine false negatives,
        `... parameter structures ...` \n
Outputs: None \n
Location: infection_process_fns.jl
"""
function increment_infection_process!(output_time_idx::Int64,
                                        count::Int64,
                                        intervention_set_itr::Int64,
                                        r_symp_test::Array{Float64,2},
                                        configuration_options::ConfigurationOptions,
                                        states::NodeStates,
                                        output::SimulationOutputs,
                                        network_params::NetworkParameters,
                                        infection_params::InfectionParameters,
                                        contact_tracing_params::ContactTracingParameters,
                                        CT_vars::ContactTracingVariables,
                                        contacts::ContactStructure,
                                        workplace_closure_params::WorkplaceClosureParameters)

    @unpack cmax, workplace_closure_active, contact_tracing_active, endtime = configuration_options
    @unpack worker_nodes = network_params
    @unpack symp_isoltime, household_isoltime, symptime, inftime = infection_params
    @unpack household_contacts, household_contacts_per_node = contacts
    @unpack WP_memory_slot = workplace_closure_params
    @unpack test_detection_prob_vec = contact_tracing_params

    time = output_time_idx - 1

    # If come to the end of latent time, move to infectious time etc
    for node_itr = 1:cmax
         # if the node has reached the end of latent infection
         if states.timelat[node_itr]>states.lattime[node_itr]
             # move to being infectious
             states.timelat[node_itr] = -1
             states.timeinf[node_itr] = 1

             # Increment time series counts
             output.numinf[output_time_idx,count,intervention_set_itr] += 1
             output.newinf[output_time_idx,count,intervention_set_itr] += 1

             # check if new infected will be asymptomatic
             if states.asymp[node_itr] > 0
                 output.newasymp[output_time_idx,count,intervention_set_itr] += 1
             end

             # check if it is worker that is newly infected
             if worker_nodes[node_itr].returned_to_work==1
                 output.workersinf[output_time_idx,count,intervention_set_itr] += 1

                 # Count workers that are asymptomatically infected
                 if states.asymp[node_itr] > 0
                     output.workersasymp[output_time_idx,count,intervention_set_itr] += 1
                 elseif workplace_closure_active==true
                     # if the newly infected worker is symptomatic, add to
                     # the workplace memory
                     workertype_ID = worker_nodes[node_itr].sector_ID
                     workplace_ID = worker_nodes[node_itr].workplace_ID
                     workplace_closure_params.workplace_memory[workertype_ID][workplace_ID,WP_memory_slot] += 1
                 end
             end
         end

         # Update node disease state time vectors
         if states.timeinf[node_itr]>inftime
             # the node becomes symptomatic (if they develop symptoms)
             states.timeinf[node_itr] = -1
             states.timesymp[node_itr] = 1

             # Increment time series counts
             output.numrep[output_time_idx,count,intervention_set_itr] += 1

             # Check if index case are symptomatic & would have zero adherence delay
             if (states.asymp[node_itr] == 0) && (states.delay_adherence[node_itr] == 0)
                 # Check if infected will isolate
                 if (states.hh_isolation[node_itr] == 1)

                     # Set that the unit has reported infection this timestep
                     states.rep_inf_this_timestep[node_itr] = 1

                     # If testing active and current node engages with test and trace
                     # check if test result will be a false negative
                     if (contact_tracing_active == true) && (CT_vars.Engage_with_CT[node_itr] == true)
                         tot_time_inf = time - states.acquired_infection[node_itr]

                         # Get relevant test sensitivity value based on time since infection
                         if tot_time_inf == 0
                             test_detection_prob = 0.
                         else
                             test_detection_prob = test_detection_prob_vec[tot_time_inf]
                         end

                         # Bernoulli trial to determine if false negative returned
                         if r_symp_test[node_itr, time] < (1 - test_detection_prob)
                             CT_vars.Test_result_false_negative[node_itr] = true
                         end

                         # If testing taking place, check if test result is positive or negative
                         # Set length of isolation
                         if CT_vars.Test_result_false_negative[node_itr] == true
                            # Release from symptomatic isolation once test result received
                            # (if does not exceed usual symptomatic isolation period)
                            length_of_symp_isol = min(CT_vars.CT_delay_until_test_result[node_itr],symp_isoltime)
                            length_of_hh_isol = min(CT_vars.CT_delay_until_test_result[node_itr], household_isoltime)
                         else
                            # Full symptomatic period spent in isolation
                            length_of_symp_isol = symp_isoltime
                            # Household members to spend full time in isolation
                            length_of_hh_isol = household_isoltime
                         end
                    else
                        # Full symptomatic period spent in isolation
                        length_of_symp_isol = symp_isoltime
                        # Household members to spend full time in isolation
                        length_of_hh_isol = household_isoltime
                    end

                     # Isolation due to presence of symptoms
                     for time_itr = 1:length_of_symp_isol
                         array_time_idx = time + time_itr
                         if array_time_idx <= (endtime + 1)
                             states.symp_isolation_array[node_itr,array_time_idx] = 1
                         end
                             # Note, in isolation_array, col 1 corresponds to time 0,
                             #                           col 2 corresponds to time 1 etc
                             # So array_time_idx = time + 1, populates array
                             # in column corresponding to day value of "time"
                     end
                 else
                     # Household members to spend full time in isolation
                     length_of_hh_isol = household_isoltime
                 end

                 # Irrespective of whether index case self-isolates,
                 # adherent members of their household may also isolate.
                 for hh = 1:household_contacts_per_node[node_itr]
                     contact_ID = household_contacts[node_itr][hh]
                     if (states.hh_isolation[contact_ID] == 1)
                         # Populate household isolation time tracker array
                         for time_itr = 1:length_of_hh_isol
                             array_time_idx = time + time_itr
                             if array_time_idx <= endtime+1
                                     states.hh_in_isolation_array[contact_ID,array_time_idx] = 1
                             end
                                 # Note, in isolation_array, col 1 corresponds to time 0,
                                 #                           col 2 corresponds to time 1 etc
                                 # So array_time_idx = time + 1, populates array
                                 # in column corresponding to day value of "time"
                         end
                     end
                 end
             end
         end

         # Check if node, if having a delayed adherence, begins adherence on current day
         if (states.timesymp[node_itr] > 1)&&((states.timesymp[node_itr]-1)==states.delay_adherence[node_itr]) # Condition for node beginning adherence on current day & has been symptomatic for at least one day
             if states.asymp[node_itr] == 0 # Check node is symptomatic and will adhere
                 if states.hh_isolation[node_itr]==1 # Check node will adhere
                     # Set that the unit has reported infection this timestep
                     states.rep_inf_this_timestep[node_itr] = 1

                     # If testing active and node will engage with test and trace
                     # check if test result will be a false negative
                     if (contact_tracing_active == true) && (CT_vars.Engage_with_CT[node_itr] == true)
                         tot_time_inf = time - states.acquired_infection[node_itr]

                         # Get relevant test sensitivity value based on time since infection
                         if tot_time_inf == 0
                             test_detection_prob = 0.
                         else
                             test_detection_prob = test_detection_prob_vec[tot_time_inf]
                         end

                         # Bernoulli trial to determine if false negative returned
                         if r_symp_test[node_itr,time] < (1 - test_detection_prob)
                             CT_vars.Test_result_false_negative[node_itr] = true
                         end

                         # If testing taking place, check if test result is positive or negative
                         # Set length of isolation.
                         # Individual shortens isolation by length of time since
                         # unwell individual began displaying symptoms.
                         if CT_vars.Test_result_false_negative[node_itr] == true
                             # Release from symptomatic isolation once test result received
                             # (if does not exceed remainder of symptomatic isolation period)
                             length_of_symp_isol = min(CT_vars.CT_delay_until_test_result[node_itr],
                                                         symp_isoltime - states.delay_adherence[node_itr])
                             length_of_hh_isol = min(CT_vars.CT_delay_until_test_result[node_itr],
                                                  household_isoltime - states.delay_adherence[node_itr])
                         else
                            # Rest of symptomatic period spent in isolation
                            length_of_symp_isol = max(0,symp_isoltime - states.delay_adherence[node_itr])
                            # Household members to spend full time in isolation
                            length_of_hh_isol = max(0,household_isoltime - states.delay_adherence[node_itr])
                         end
                     else
                         # Rest of symptomatic period spent in isolation
                         length_of_symp_isol = max(0,symp_isoltime - states.delay_adherence[node_itr])
                         # Household members to spend full time in isolation
                         length_of_hh_isol = max(0,household_isoltime - states.delay_adherence[node_itr])
                     end

                     # Isolation due to presence of symptoms
                     for time_itr = 1:length_of_symp_isol
                             # Shorten the isolation period by time already elapsed since
                             # symptom onset
                         array_time_idx = time + time_itr
                         if array_time_idx <= (endtime + 1)
                             states.symp_isolation_array[node_itr,array_time_idx] = 1
                         end
                     end
                 else
                     # Household members to spend full time in isolation
                     length_of_hh_isol = max(0,household_isoltime - states.delay_adherence[node_itr])
                 end

                 # Irrespective of whether index case self-isolates,
                 # adherent members of their household may also isolate.
                 for hh = 1:household_contacts_per_node[node_itr]
                     contact_ID = household_contacts[node_itr][hh]
                     if (states.hh_isolation[contact_ID]==1)
                         # Populate household isolation time tracker array
                         for time_itr = 1:length_of_hh_isol
                             array_time_idx = time + time_itr
                             if array_time_idx <= (endtime + 1)
                                 states.hh_in_isolation_array[contact_ID,array_time_idx] = 1
                             end
                                 # Note, in isolation_array, col 1 corresponds to time 0,
                                 #                           col 2 corresponds to time 1 etc
                                 # So array_time_idx = time + 1, populates array
                                 # in column corresponding to day value of "time"
                         end
                     end
                 end
             end
         end

         # Check if node has reached end of symptom period
         if states.timesymp[node_itr]>symptime
             states.timesymp[node_itr] = -1
         end
    end
end

"""
    get_daily_isol_atwork_status!(args)

Record whether each node is physically present at work or not, or in isolation, on the current day.

Inputs: `time` - current time in simulation,
        `output_time_idx` - current index in output array (corresponds to time + 1),
        `... parameter structures ...` \n
Outputs: None \n
Location: infection_process_fns.jl
"""
function get_daily_isol_atwork_status!(time::Int64,
                                        output_time_idx::Int64,
                                        configuration_options::ConfigurationOptions,
                                        states::NodeStates,
                                        contacts::ContactStructure,
                                        network_params::NetworkParameters)

    @unpack cmax, isolation = configuration_options
    @unpack worker_nodes, workplace_info = network_params

    # record whether nodes are in isolation
    for node_itr = 1:cmax
        # Now record whether node is in isolation on current day
        if (isolation>0) &&
            ((states.hh_in_isolation_array[node_itr,output_time_idx] == 1) ||
            (states.symp_isolation_array[node_itr,output_time_idx] == 1) ||
            (states.CT_isolation_array[node_itr,output_time_idx] == 1))
            # Default value is 0. So only update if node is in isolation for any reason
            states.daily_record_inisol[time,node_itr] = 1
        end

        # Record whether node is at workplace on current day
        # Needs to be returned to work, and a day where at workplace
        # AND not in isolation
        # AND whole workplace not closed
        workertype_ID = worker_nodes[node_itr].sector_ID
        workplace_ID = worker_nodes[node_itr].workplace_ID
        node_workplace_info = workplace_info[workertype_ID][workplace_ID]
        if (worker_nodes[node_itr].returned_to_work==1) &&
            (states.atwork[node_itr,time] == true) &&
            (states.daily_record_inisol[time,node_itr] == false) &&
            (node_workplace_info.workplace_open == true)

            # Default value is 0. So only update if node is at workplace
            states.daily_record_atworkplace[time,node_itr] = 1
        end
    end

end

"""
    transmit_infections!(args)

Iterate over all infectious nodes and transmit infections in all relevant settings.

Inputs: `rng` - random number generator,
        `time` - current time in simulation,
        `count` - number of current replicate,
        `intervention_set_itr` - number of current intervention set,
        `... parameter structures ...` \n
Outputs: None \n
Location: infection_process_fns.jl
"""
function transmit_infections!(rng::MersenneTwister,
                                time::Int64,
                                count::Int64,
                                intervention_set_itr::Int64,
                                configuration_options::ConfigurationOptions,
                                states::NodeStates,
                                contacts::ContactStructure,
                                network_params::NetworkParameters,
                                infection_params::InfectionParameters,
                                output::SimulationOutputs)

    @unpack cmax, CS_active_flag = configuration_options
    @unpack dist_infectivity, iso_trans_scaling,
            asymp_trans_scaling, probasymp, scaling_household,
            scaling_work_static, scaling_work_dynamic,
            scaling_random, scaling_social, CS_scale_transrisk,
            inftime = infection_params
    @unpack worker_nodes = network_params
    @unpack household_contacts, work_contacts_same_workplace_per_node,
            work_contacts_other_workplace_per_node,
            dynamic_random_contacts = contacts

    # Initialise undefined array
    # Used transmit_over! function
    undefined_array = Array{Int64,2}(undef,0,0)

    # Iterate over nodes that may be able to transmit infection
    for node_itr = 1:cmax
        if ((states.timeinf[node_itr]>0) | (states.timesymp[node_itr]>0))
            # Only enter loop if node is capable of transmitting infection

            # find the total time infectious
            if states.timeinf[node_itr]>0
                tot_time_infectious = states.timeinf[node_itr]
            else
                tot_time_infectious = states.timesymp[node_itr]+inftime
            end

            # find the infectiousness
            infectiousness = dist_infectivity[tot_time_infectious]
            current_worker = worker_nodes[node_itr]
            if states.asymp[node_itr]>0 # Asymptomatic
                transtemp_household = current_worker.transrisk_household*scaling_household*infectiousness*asymp_trans_scaling
                transtemp_work_static = current_worker.transrisk_static_work*scaling_work_static*infectiousness*asymp_trans_scaling
                transtemp_work_dynamic = current_worker.transrisk_dynamic_work*scaling_work_dynamic*infectiousness*asymp_trans_scaling
                transtemp_social = current_worker.transrisk_social*scaling_social*infectiousness*asymp_trans_scaling
                transtemp_random = current_worker.transrisk_random*scaling_random*infectiousness*asymp_trans_scaling
            else
                if states.timesymp[node_itr]>0  # symptomatic & less infectious due to cautionary behaviour
                    transtemp_household = current_worker.transrisk_household*scaling_household*infectiousness*iso_trans_scaling
                    transtemp_work_static = current_worker.transrisk_static_work*scaling_work_static*infectiousness*iso_trans_scaling
                    transtemp_work_dynamic = current_worker.transrisk_dynamic_work*scaling_work_dynamic*infectiousness*iso_trans_scaling
                    transtemp_social = current_worker.transrisk_social*scaling_social*infectiousness*iso_trans_scaling
                    transtemp_random = current_worker.transrisk_random*scaling_random*infectiousness*iso_trans_scaling
                else  # infected, unscaled infectiousness
                    transtemp_household = current_worker.transrisk_household*scaling_household*infectiousness
                    transtemp_work_static = current_worker.transrisk_static_work*scaling_work_static*infectiousness
                    transtemp_work_dynamic = current_worker.transrisk_dynamic_work*scaling_work_dynamic*infectiousness
                    transtemp_social = current_worker.transrisk_social*scaling_social*infectiousness
                    transtemp_random = current_worker.transrisk_random*scaling_random*infectiousness
                end
            end

            # Infection check for other household members
            # Transmit over household_contacts[node_itr]
            # checking that contacts are susceptible
            n_hh_contacts = length(household_contacts[node_itr])
            if n_hh_contacts > 0
                   transmit_over!(transtemp_household,output,states,probasymp,rng,time,count,intervention_set_itr,
                                        node_itr,
                                        household_contacts[node_itr],
                                        "household")
            end

            # Check individual is not in household isolation.
            # If not, can see if transmission across cohort, society,
            # dynamic household contacts occured.
            if (states.daily_record_inisol[time,node_itr] == false)
                # Check if node is at workplace or not
                if states.daily_record_atworkplace[time,node_itr] == true
                    if contacts.social_contacts_per_node[node_itr] > 0
                        # Apply contact scaling during transmission if node is adhering
                        if states.hh_isolation[node_itr] == 1
                            social_contact_scaling_node_itr = network_params.social_contact_scaling
                        else
                            social_contact_scaling_node_itr = -1.
                        end
                        # Satisfied condition that node_itr may have social links
                        # transmit over contacts.workday_social_contacts_by_day[time,node_itr]
                        # checking that contacts are susceptible and not isolating
                        transmit_over!(transtemp_social,output,states,probasymp,rng,time,count,intervention_set_itr,
                                node_itr,
                                contacts.workday_social_contacts_by_day[time,node_itr],
                                "social",
                                inisol=states.daily_record_inisol,
                                atwork=undefined_array,
                                social_contact_scaling=social_contact_scaling_node_itr)
                    end

                    # Check if workplace is Covid-Secure
                    # If so, modify the transmission risk and use CS contacts
                    # Also, non-workplace worker contacts made
                    workertype_ID = current_worker.sector_ID
                    workplace_ID = current_worker.workplace_ID
                    current_workplace_info = network_params.workplace_info[workertype_ID][workplace_ID]
                    if (CS_active_flag == true) && (current_workplace_info.covid_secure == true)
                        transtemp_work_static *= CS_scale_transrisk[workertype_ID]
                        transtemp_work_dynamic *= CS_scale_transrisk[workertype_ID]
                        # These are the regular worker contacts in the same workplace
                        # in a CS setting
                        n_work_contacts_same_workplace_CS = contacts.work_contacts_same_workplace_per_node_CS[node_itr]
                        if n_work_contacts_same_workplace_CS > 0
                            # Satisfied condition that node_itr has any workday links
                            # transmit over contacts.work_contacts_same_workplace_CS[time,node_itr]
                            # checking that contacts are susceptible, not isolating and at work
                            transmit_over!(transtemp_work_static,output,states,probasymp,rng,time,count,intervention_set_itr,
                                    node_itr,
                                    contacts.work_contacts_same_workplace_CS[node_itr],
                                    "work",
                                    inisol = states.daily_record_inisol,
                                    atwork = states.daily_record_atworkplace)
                        end
                    else
                        # These are the regular worker contacts in the same workplace
                        n_work_contacts_same_workplace = work_contacts_same_workplace_per_node[node_itr]
                        if n_work_contacts_same_workplace > 0
                            # Satisfied condition that node_itr has any workday links
                            # transmit over contacts.work_contacts_same_workplace[node_itr]
                            # checking that contacts are susceptible, not isolating and at work
                            transmit_over!(transtemp_work_static,output,states,probasymp,rng,time,count,intervention_set_itr,
                                    node_itr,
                                    contacts.work_contacts_same_workplace[node_itr],
                                    "work",
                                    inisol = states.daily_record_inisol,
                                    atwork = states.daily_record_atworkplace)
                        end

                        # These are the regular worker contacts with workers from
                        # other workplaces
                        n_work_contacts_other_workplace = work_contacts_other_workplace_per_node[node_itr]
                        if n_work_contacts_other_workplace > 0
                            # Satisfied condition that node_itr has any workday links
                            # transmit over contacts.work_contacts_other_workplace[node_itr]
                            # checking that contacts are susceptible, not isolating and at work
                            # also check their workplace is not CS - if CS then no close contact occurs.
                            transmit_over_other_workplace!(transtemp_work_static,output,states,probasymp,rng,time,count,intervention_set_itr,
                                                            node_itr,
                                                            contacts.work_contacts_other_workplace[node_itr],
                                                            states.daily_record_inisol,
                                                            states.daily_record_atworkplace,
                                                            network_params)
                        end
                    end

                    # Add dynamic links, if returned to work
                    # transmit over contacts.dynamic_worker_contacts[time,node_itr]
                    # checking that contacts are susceptible and not isolating
                    transmit_over!(transtemp_work_dynamic,output,states,probasymp,rng,time,count,intervention_set_itr,
                            node_itr,
                            contacts.dynamic_worker_contacts[time,node_itr],
                            "work",
                            inisol = states.daily_record_inisol,
                            atwork=undefined_array,
                            dynamic_contact = 1)

                else # otherwise the node is not at work

                    # Add in social contacts on non-work days
                    if contacts.social_contacts_per_node[node_itr] > 0
                        # Apply contact scaling during transmission if node is adhering
                        if states.hh_isolation[node_itr] == 1
                            social_contact_scaling_node_itr = network_params.social_contact_scaling
                        else
                            social_contact_scaling_node_itr = -1.
                        end
                        # Satisfied condition that node_itr may have social links
                        # transmit over contacts.nonworkday_social_contacts_by_day[time,node_itr]
                        # checking that contacts are susceptible and not isolating
                        transmit_over!(transtemp_social,output,states,probasymp,rng,time,count,intervention_set_itr,
                                node_itr,
                                contacts.nonworkday_social_contacts_by_day[time,node_itr],
                                "social",
                                inisol=states.daily_record_inisol,
                                atwork=undefined_array,
                                social_contact_scaling=social_contact_scaling_node_itr)
                    end
                end

                # If not in isolation, also include random daily contacts
                n_random_contacts = length(dynamic_random_contacts[time, node_itr])
                if n_random_contacts > 0
                    # Apply contact scaling during transmission if node is adhering
                    if states.hh_isolation[node_itr] == 1
                        random_contact_scaling_node_itr = network_params.random_contact_scaling
                    else
                        random_contact_scaling_node_itr = -1.
                    end
                   transmit_over!(transtemp_random,output,states,probasymp,rng,time,count,intervention_set_itr,
                                        node_itr,
                                        dynamic_random_contacts[time, node_itr],
                                        "other",
                                        random_contact_scaling=random_contact_scaling_node_itr)
                end
            end
        end
    end
end

"""
    transmit_over!(args)

For a single node and specified contacts, check transmission potential and randomly transmit infection.

Can be applied to all transmission settings apart from contacts from another workplace.

Inputs: `transmission_risk` - relevant transmission probability for current node and context,
        `... parameter structures ...`,
        `probasymp` - probability that new infections will be asymptomatic,
        `rng` - random number generator,
        `time` - current time in simulation,
        `count` - number of current replicate,
        `intervention_set_itr` - number of current intervention set,
        `infecting_by` - ID of node that is infecting others,
        `contacts_to_check` - array of node IDs that may be infected,
        `dynamic_contact` - (optional) flag if contact is dynamic or not,
        `inisol` - (optional) flag if nodes are in isolation,
        `atwork` - (optional) flag if nodes are at work (according to work schedule),
        `transmission_setting` - (optional) current transmission context (e.g. social),
        `social_contact_scaling` - (optional) scaling on social contacts (contacts can be 'deactivated' to reduce number of contacts),
        `random_contact_scaling` - (optional) scaling on random contacts (contacts can be 'deactivated' to reduce number of contacts) \n
Outputs: None \n
Location: infection_process_fns.jl
"""
function transmit_over!(transmission_risk::Float64,
                    output::SimulationOutputs,
                    states::NodeStates,
                    probasymp::Float64,
                    rng::MersenneTwister,
                    time::Int64,
                    count::Int64,
                    intervention_set_itr::Int64,
                    infecting_by::Int64,
                    contacts_to_check::Array{Int64,1},
                    transmission_setting::String;
                    dynamic_contact::Int64=0,
                    inisol::Array{Int64,2} = Array{Int64,2}(undef,0,0),
                    atwork::Array{Int64,2} = Array{Int64,2}(undef,0,0),
                    social_contact_scaling::Float64 = -1.,
                    random_contact_scaling::Float64 = -1.)

    for contact_itr = 1:length(contacts_to_check) # Iterate over each contact
        infecting_to = contacts_to_check[contact_itr]

        contact_active_check = rand(rng)
        # Check if infecting_to is susceptible, not isolating and at work
        # only checks isolation or at work if given the arguments inisol / atwork
        if (states.timelat[infecting_to]==0) &&
            (isassigned(inisol) == false || inisol[time,infecting_to] == false) &&
            (isassigned(atwork) == false || atwork[time,infecting_to] == true) &&
            ((social_contact_scaling < 0) || (contact_active_check <= social_contact_scaling)) &&
            ((random_contact_scaling < 0) || (contact_active_check <= random_contact_scaling))

            if rand(rng) < transmission_risk
                states.timelat[infecting_to] = 1
                states.infected_by[infecting_to] = infecting_by
                output.num_infected[infecting_by,count,intervention_set_itr] += 1
                states.acquired_infection[infecting_to] = time

                # adjust Rt(t) = mean number of infections generated by nodes that were infected at time t
                if states.acquired_infection[infecting_by]>0
                    # Offset time to array indexing. Row 1 of output.Rt is for day 0, Row 2 is for day 1 etc
                    output.Rt[(states.acquired_infection[infecting_by]+1),count,intervention_set_itr] += 1
                end

                # if this was from an initial infection, modify the generation time
                # this sum will be divided by the total number of secondary initial infections
                if states.infected_by[infecting_by]==-1
                    output.mean_init_generation_time[count,intervention_set_itr] += time + states.lattime[infecting_by]
                end

                # Check if infection will be asymptomatic
                if rand(rng) < probasymp
                    states.asymp[infecting_to] = 1
                end

                # Update latent event counter
                output.numlat[time+1,count,intervention_set_itr] += 1

                if transmission_setting == "work"
                    if dynamic_contact==1
                        # Record transmission occurring in work-dynamic setting
                        output.transmission_setting[time+1,count,intervention_set_itr,3] += 1
                        # if this was a dynamic contact, update counter
                        output.dynamic_infection_count[time+1,count,intervention_set_itr] += 1
                    else
                        # Record transmission occurring in work-static setting
                        output.transmission_setting[time+1,count,intervention_set_itr,2] += 1
                    end
                elseif transmission_setting == "social"
                    # Record transmission occurring in social setting
                    output.transmission_setting[time+1,count,intervention_set_itr,1] += 1
                elseif transmission_setting == "household"
                    # Record transmission occurring in household setting
                    output.transmission_setting[time+1,count,intervention_set_itr,4] += 1
                elseif transmission_setting == "other"
                    # Record transmission occurring in other setting
                    output.transmission_setting[time+1,count,intervention_set_itr,5] += 1
                else
                    error("Invalid transmission setting!")
                end
            end
        end
    end
end

"""
    transmit_over_other_workplace!(args)

For a single node and specified contacts, check transmission potential and randomly transmit infection.

Only applied to transmission between contacts from another workplace.

Inputs: `transmission_risk` - relevant transmission probability for current node and context,
        `... parameter structures ...`,
        `probasymp` - probability that new infections will be asymptomatic,
        `rng` - random number generator,
        `time` - current time in simulation,
        `count` - number of current replicate,
        `intervention_set_itr` - number of current intervention set,
        `infecting_by` - ID of node that is infecting others,
        `contacts_to_check` - array of node IDs that may be infected,
        `inisol` - flag if nodes are in isolation,
        `atwork` - flag if nodes are at work (according to work schedule) \n
Outputs: None \n
Location: infection_process_fns.jl
"""
function transmit_over_other_workplace!(transmission_risk::Float64,
                                            output::SimulationOutputs,
                                            states::NodeStates,
                                            probasymp::Float64,
                                            rng::MersenneTwister,
                                            time::Int64,
                                            count::Int64,
                                            intervention_set_itr::Int64,
                                            infecting_by::Int64,
                                            contacts_to_check::Array{Int64,1},
                                            inisol::Array{Int64,2},
                                            atwork::Array{Int64,2},
                                            network_params::NetworkParameters)

    for contact_itr = 1:length(contacts_to_check) # Iterate over each contact
        infecting_to = contacts_to_check[contact_itr]

        # Get information on the contacts workplace
        contact_sector_ID = network_params.worker_nodes[infecting_to].sector_ID
        contact_workplace_ID = network_params.worker_nodes[infecting_to].workplace_ID
        contact_workplace_info = network_params.workplace_info[contact_sector_ID][contact_workplace_ID]

        # Get contacts workplace CS status
        contact_workplace_CS_bool = contact_workplace_info.covid_secure

        # Check if infecting_to is susceptible, not isolating and at work
        # only checks isolation at at work if given the arguments inisol / atwork
        if (states.timelat[infecting_to]==0) &&
            (isassigned(inisol) == false || inisol[time,infecting_to] == false) &&
            (isassigned(atwork) == false || atwork[time,infecting_to] == true) &&
            (contact_workplace_CS_bool==false)

            if rand(rng) < transmission_risk
                states.timelat[infecting_to] = 1
                states.infected_by[infecting_to] = infecting_by
                output.num_infected[infecting_by,count,intervention_set_itr] += 1
                states.acquired_infection[infecting_to] = time

                # adjust Rt(t) = mean number of infections generated by nodes that were infected at time t
                if states.acquired_infection[infecting_by]>0
                    # Offset time to array indexing. Row 1 of output.Rt is for day 0, Row 2 is for day 1 etc
                    output.Rt[(states.acquired_infection[infecting_by]+1),count,intervention_set_itr] += 1
                end

                # if this was from an initial infection, modify the generation time
                # this sum will be divided by the total number of secondary initial infections
                if states.infected_by[infecting_by]==-1
                    output.mean_init_generation_time[count,intervention_set_itr] += time + states.lattime[infecting_by]
                end

                # Check if infection will be asymptomatic
                if rand(rng) < probasymp
                    states.asymp[infecting_to] = 1
                end

                # Update latent event counter
                output.numlat[time+1,count,intervention_set_itr] += 1

                # Record transmission occurring in work-static setting
                output.transmission_setting[time+1,count,intervention_set_itr,2] += 1

            end
        end
    end
end

"""
    assign_prev_isol_outputs!(args)

Store states of nodes in output vectors.

Inputs: `output_time_idx` - current index in output array (corresponds to time + 1),
        `count` - number of current replicate,
        `intervention_set_itr` - number of current intervention set,
        `... parameter structures ...` \n
Outputs: None \n
Location: infection_process_fns.jl
"""
function assign_prev_isol_outputs!(output_time_idx::Int64,
                                    count::Int64,
                                    intervention_set_itr::Int64,
                                    configuration_options::ConfigurationOptions,
                                    states::NodeStates,
                                    output::SimulationOutputs)

    @unpack cmax, contact_tracing_active = configuration_options

    # For this timestep, get number isolating
    # and if they are latent infected or infectious on current timestep
    for node_itr=1:cmax
        # Initialise isolation status flag
        isolating_for_any_reason = false

        # Isolating due to housemate symptoms
        if states.hh_in_isolation_array[node_itr,output_time_idx] == 1
            output.num_household_isolating[output_time_idx,count,intervention_set_itr] += 1
            isolating_for_any_reason = true
        end

        # Isolating due to symptoms
        if states.symp_isolation_array[node_itr,output_time_idx] == 1
            output.num_symp_isolating[output_time_idx,count,intervention_set_itr] += 1
            isolating_for_any_reason = true
        end

        # Isolating as close contact of positive test symptomtic
        if contact_tracing_active == true
            if states.CT_isolation_array[node_itr,output_time_idx] == 1
                output.num_isolating_CTcause[output_time_idx,count,intervention_set_itr] += 1
                isolating_for_any_reason = true
            end
        end

        # Isolating for any reason
        if isolating_for_any_reason == true
            output.num_isolating[output_time_idx,count,intervention_set_itr] += 1
        end

        # Check if latently infected
        if states.timelat[node_itr]>0
            output.prevlat[output_time_idx,count,intervention_set_itr] += 1
        end

        # In presymptomatic infectious period.
        if (states.timeinf[node_itr]>0)
            if states.asymp[node_itr] > 0 # asymptomatic
                output.prevasymp[output_time_idx,count,intervention_set_itr] += 1
            else # will be symptomatic
                output.prevpresymp[output_time_idx,count,intervention_set_itr] += 1
            end
        end

        # After presymp period, check if symptomatic or asymptomatic
        if (states.timesymp[node_itr]>0)
            if states.asymp[node_itr] > 0 # asymptomatic
                output.prevasymp[output_time_idx,count,intervention_set_itr] += 1
            else # symptomatic
                output.prevsymp[output_time_idx,count,intervention_set_itr] += 1
            end
        end

        # Check if recovered
        if states.timesymp[node_itr] == -1
            output.prevrec[output_time_idx,count,intervention_set_itr] += 1
        end
    end
end
