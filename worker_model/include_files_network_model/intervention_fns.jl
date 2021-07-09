"""
File containing functions allowing the implementation of interventions during a simulation.
Interventions can be timed or triggered by a condition.
"""

"""
   affect_intervention!(args)

Implement changes to parameter structures according to given InterventionVariables structure, from current time until end.

Changes are applied in two stages: i) parameters ('settings') are changed, then ii) relevant network components are regenerated if necessary.

Inputs: `time` - current time in simulation,
         `intervention` - InterventionVariables structure defining the intervention to be implemented,
         `intervention_set_itr` - number of current intervention set,
         `... parameter structures ...`,
         `preintervention_settings` - InterventionVariables structure used to reverse current intervention set,
         `settings_changed_by_intervention_set` - array of settings changed by this intervention set,
         `contact_component_store` - storage for preintervention contacts, if intervention affects contact structure,
         `rng` - random number generator \n
Outputs: None \n
Location: intervention_fns.jl
"""
function affect_intervention!(time::Int64,
                              intervention::InterventionVariables,
                              intervention_set_itr::Int64,
                               configuration_options::ConfigurationOptions,
                               infection_params::InfectionParameters,
                               contact_tracing_params::ContactTracingParameters,
                               CT_vars::ContactTracingVariables,
                               contacts::ContactStructure,
                               network_params::NetworkParameters,
                               workplace_generation_params::WorkplaceGenerationParameters,
                               states::NodeStates,
                               workplace_closure_params::WorkplaceClosureParameters,
                               preintervention_settings::InterventionVariables,
                               settings_changed_by_intervention_set::Array{Symbol,1},
                               contact_component_store::Array{PreinterventionComponentStore,1},
                               rng::MersenneTwister)

   settings_changed_by_intervention::Array{Symbol,1} = Symbol[]

   # In order to regenerate some components, store unchanged setting values
   old_symp_isoltime::Int64 = infection_params.symp_isoltime
   old_household_isoltime::Int64 = infection_params.household_isoltime
   old_workplace_CT_memory::Int64 = workplace_closure_params.workplace_CT_memory

   daily_social_contact_setting_names = [:network_generation_method_dynamic_social,
                                    :social_workday_dd, :social_nonworkday_dd,
                                    :group_limit, :dynamic_time_frame,
                                    :n_groups_per_day_distribution,
                                    :n_households_per_group]
   old_daily_social_contact_settings = [getfield(network_params, daily_social_contact_setting_names[ii]) for ii=1:length(daily_social_contact_setting_names)]

   # apply intervention to parameter structures and record preintervention settings
   change_parameters_by_intervention!(intervention,
                                     configuration_options,
                                     infection_params,
                                     contact_tracing_params,
                                     network_params,
                                     workplace_generation_params,
                                     workplace_closure_params,
                                     preintervention_settings,
                                     settings_changed_by_intervention_set,
                                     settings_changed_by_intervention)

   # Error checks for settings changed
   if length(intersect([:CT_engagement, :CT_delay_until_test_result_pmf,
      :CT_days_before_symptom_included, :CT_caused_isol_limit, :perform_CT_from_infector,
      :test_detection_prob_vec, :prob_backwards_CT, :infector_engage_with_CT_prob],
                                                      settings_changed_by_intervention)) > 0
      if configuration_options.contact_tracing_active == false
         error("Attempted to change CT setting(s) whilst CT inactive")
      end
   end
   if length(intersect([:prob_backwards_CT, :infector_engage_with_CT_prob],
                                                      settings_changed_by_intervention)) > 0
      if configuration_options.perform_CT_from_infector == false
         error("Attempted to change backwards CT setting(s) whilst backwards CT inactive")
      end
   end
   if length(intersect([:workplace_CT_threshold, :workplace_CT_memory,
      :time_WC], settings_changed_by_intervention)) > 0
      if configuration_options.workplace_closure_active == false
         error("Attempted to change workplace closure setting(s) whilst workplace closures inactive")
      end
   end
   if length(intersect([:symp_isoltime, :household_isoltime,
      :delay_adherence_pmf, :iso_trans_scaling], settings_changed_by_intervention)) > 0
      if configuration_options.isolation < 1
         error("Attempted to change isolation setting(s) whilst isolation inactive")
      end
   end
   if length(intersect([:CS_scale_transrisk, :CS_team_size], settings_changed_by_intervention)) > 0
      if configuration_options.CS_active_flag == false
         error("Attempted to change CS workplace setting(s) whilst CS workplaces inactive")
      end
   end

    if (:contact_tracing_active ∈ settings_changed_by_intervention) &&
       (configuration_options.contact_tracing_active == true) &&
       (CT_vars.cmax == 0)
         initialise_CT_vars!(configuration_options.cmax, CT_vars)
         reinitialise_CT_vars!(CT_vars, configuration_options.cmax, rng,
                                     contact_tracing_params,
                                     states,
                                     infection_params)
    end

    if (:workplace_closure_active ∈ settings_changed_by_intervention) &&
      (configuration_options.workplace_closure_active == true) &&
      (isassigned(workplace_closure_params.num_workplaces) == false)
         initialise_workplace_closure_params!(workplace_closure_params,
                                             configuration_options,
                                             network_params)
    end

   if length(intersect(settings_changed_by_intervention, [:sameday, :ton, :toff])) > 0
      # Redefine when workers are at work
      populate_atwork!(rng, states.atwork,
                        configuration_options)
   end

   if :isolation ∈ settings_changed_by_intervention
      redefine_isolation_arrays!(time,
                                 configuration_options,
                                 states,
                                 infection_params,
                                 contacts)
   end

   if :workpercent ∈ settings_changed_by_intervention
      redefine_returned_to_work!(configuration_options,
                                 network_params.worker_nodes,
                                 workplace_generation_params,
                                 network_params)
   end

   if :CT_engagement ∈ settings_changed_by_intervention
      redefine_CT_engagement!(configuration_options,
                              contact_tracing_params,
                              CT_vars,
                              states,
                              infection_params)
   end

   if :adherence ∈ settings_changed_by_intervention
      redefine_adherence!(configuration_options.cmax,
                           states,
                           infection_params)
      redefine_CT_engagement!(configuration_options,
                              contact_tracing_params,
                              CT_vars,
                              states,
                              infection_params)
   end

   if :CS_active_flag ∈ settings_changed_by_intervention
      redefine_covid_secure!(rng,
                              configuration_options,
                              network_params,
                              contacts)
   end

   if :CS_team_size ∈ settings_changed_by_intervention
      redefine_CS_team_size!(rng,
                              configuration_options,
                              network_params,
                              contacts)
   end

   if length(intersect(settings_changed_by_intervention,
                        [:network_generation_method_dynamic_social,
                        :social_workday_dd, :social_nonworkday_dd,
                        :group_limit, :dynamic_time_frame,
                        :n_groups_per_day_distribution,
                        :n_households_per_group])) > 0
      redefine_daily_social_contacts!(rng,
                                       time,
                                       intervention_set_itr,
                                       configuration_options,
                                       network_params,
                                       contacts,
                                       contact_component_store,
                                       daily_social_contact_setting_names,
                                       old_daily_social_contact_settings)
   end

   if :sector_open ∈ settings_changed_by_intervention
      redefine_workplace_open!(network_params)
   end

   if :symp_isoltime ∈ settings_changed_by_intervention
      redefine_symp_isolation!(time,
                              old_symp_isoltime,
                              configuration_options,
                              states,
                              infection_params,
                              CT_vars)
   end

   if :household_isoltime ∈ settings_changed_by_intervention
      redefine_household_isolation!(time,
                              old_household_isoltime,
                              configuration_options,
                              states,
                              infection_params,
                              CT_vars)
   end

   if :delay_adherence_pmf ∈ settings_changed_by_intervention
      redefine_delay_adherence!(rng,
                                 configuration_options,
                                 states,
                                 infection_params,
                                 contact_tracing_params,
                                 CT_vars)
   end

   if :CT_delay_until_test_result_pmf ∈ settings_changed_by_intervention
      redefine_CT_delay_until_test_result!(rng,
                                          contact_tracing_params,
                                          CT_vars)
   end

   if :CT_days_before_symptom_included ∈ settings_changed_by_intervention
      redefine_relevant_prev_days_for_CT!(configuration_options,
                                          CT_vars,
                                          contact_tracing_params,
                                          states)
   end

   if :CT_caused_isol_limit ∈ settings_changed_by_intervention
      redefine_CT_isolation!(configuration_options,
                              contact_tracing_params,
                              CT_vars,
                              states)
   end

   if :workplace_CT_threshold ∈ settings_changed_by_intervention
      redefine_workplace_CT_threshold!(configuration_options,
                                 workplace_closure_params,
                                 network_params)
   end

   if :workplace_CT_memory ∈ settings_changed_by_intervention
      redefine_workplace_CT_memory!(configuration_options,
                                 workplace_closure_params,
                                 old_workplace_CT_memory)
   end

   return nothing
end

"""
   change_parameters_by_intervention!(args)

Change settings in parameter structures according to given InterventionVariable structure, and store pre-intervention settings.

Pre-intervention settings are stored as an InterventionVariable structure, which, when applied, changes settings back to their original values.
Pre-intervention settings are the settings before ANY intervention in the simulation, not necessarily the current intervention.

Inputs:  `intervention` - InterventionVariables structure defining the intervention to be implemented,
         `... parameter structures ...`,
         `preintervention_settings` - InterventionVariables structure used to reverse current intervention set,
         `settings_changed_by_intervention_set` - array of settings changed by this intervention set,
         `settings_changed_by_intervention` - array of settings changed by this intervention \n
Outputs: None \n
Location: intervention_fns.jl
"""
function change_parameters_by_intervention!(intervention::InterventionVariables,
                                           configuration_options::ConfigurationOptions,
                                           infection_params::InfectionParameters,
                                           contact_tracing_params::ContactTracingParameters,
                                           network_params::NetworkParameters,
                                           workplace_generation_params::WorkplaceGenerationParameters,
                                           workplace_closure_params::WorkplaceClosureParameters,
                                           preintervention_settings::InterventionVariables,
                                           settings_changed_by_intervention_set::Array{Symbol,1},
                                           settings_changed_by_intervention::Array{Symbol,1})

    all_settings = fieldnames(InterventionVariables)
    n_settings = length(all_settings)

    for setting_idx = 1:n_settings

        setting_name = all_settings[setting_idx]

        if (setting_name != :start_time) && (setting_name != :reset_time)

            intervention_obj = getfield(intervention, setting_name)

            if isassigned(intervention_obj)

                if setting_name ∉ settings_changed_by_intervention_set
                    push!(settings_changed_by_intervention_set, setting_name)
                    push!(settings_changed_by_intervention, setting_name)
                    first_time_change = true
                else
                    push!(settings_changed_by_intervention, setting_name)
                    first_time_change = false
                end

                if hasfield(ConfigurationOptions, setting_name)
                    if first_time_change
                        preintervention_obj = getfield(configuration_options, setting_name)
                        setfield!(preintervention_settings, setting_name, [preintervention_obj])
                    end
                    setfield!(configuration_options, setting_name, intervention_obj[1])
                elseif hasfield(NetworkParameters, setting_name)
                    if first_time_change
                        preintervention_obj = getfield(network_params, setting_name)
                        setfield!(preintervention_settings, setting_name, [preintervention_obj])
                    end
                    setfield!(network_params, setting_name, intervention_obj[1])
                elseif hasfield(ContactTracingParameters, setting_name)
                    if first_time_change
                        preintervention_obj = getfield(contact_tracing_params, setting_name)
                        setfield!(preintervention_settings, setting_name, [preintervention_obj])
                    end
                    setfield!(contact_tracing_params, setting_name, intervention_obj[1])
                elseif hasfield(InfectionParameters, setting_name)
                    if first_time_change
                        preintervention_obj = getfield(infection_params, setting_name)
                        setfield!(preintervention_settings, setting_name, [preintervention_obj])
                    end
                    setfield!(infection_params, setting_name, intervention_obj[1])
                elseif hasfield(WorkplaceGenerationParameters, setting_name)
                    if first_time_change
                        preintervention_obj = getfield(workplace_generation_params, setting_name)
                        setfield!(preintervention_settings, setting_name, [preintervention_obj])
                    end
                    setfield!(workplace_generation_params, setting_name, intervention_obj[1])
                elseif hasfield(WorkplaceClosureParameters, setting_name)
                    if first_time_change
                        preintervention_obj = getfield(workplace_closure_params, setting_name)
                        setfield!(preintervention_settings, setting_name, [preintervention_obj])
                    end
                    setfield!(workplace_closure_params, setting_name, intervention_obj[1])
                else
                    error("Intervention attempted to change setting that does not exist: $setting_name")
                end
            end
        end
    end
    return nothing
end

"""
   redefine_covid_secure!(args)

Apply COVID-secure flag to every workplace (true/false) and, if necessary, construct COVID-secure workplace contact network.

Only called if 'CS_active_flag' is changed via intervention.
CS workplace contacts only generated if flag is switched to true and contacts have not been
previously generated (same contacts are kept between replicates, but not configurations).

Inputs:  `rng` - random number generator,
         `... parameter structures ...` \n
Outputs: None \n
Location: intervention_fns.jl
"""
function redefine_covid_secure!(rng::MersenneTwister,
                                 configuration_options::ConfigurationOptions,
                                 network_params::NetworkParameters,
                                 contacts::ContactStructure)

   # Reconstruct work contacts to conform to capped team sizes

   @unpack cmax, workertypes, CS_active_flag = configuration_options
   @unpack CS_contacts_previously_generated = network_params

   if (CS_active_flag == true) && (CS_contacts_previously_generated == false)
      # Initialise CS related variables
      contacts.work_contacts_same_workplace_CS = Array{Array{Int64,1},1}(undef,cmax)
      contacts.work_contacts_same_workplace_per_node_CS = zeros(Int64,cmax)
      for node_itr = 1:cmax
         contacts.work_contacts_same_workplace_CS[node_itr] = Int64[]
      end

      # Construct COVID-secure work contacts
      CS_workplace_generation!(network_params.worker_nodes,
                                 network_params.nodes_by_workplace,
                                 contacts.work_contacts_same_workplace_CS,
                                 contacts.work_contacts_same_workplace_per_node_CS,
                                 network_params,
                                 rng)

   end

   # Set all workplaces to be COVID-secure or not
   for worker_grp_idx = 1:workertypes
      workplace_count = length(network_params.workplace_info[worker_grp_idx])
      for workplace_itr = 1:workplace_count
         network_params.workplace_info[worker_grp_idx][workplace_itr].covid_secure = CS_active_flag
      end
   end

   return nothing
end

"""
   redefine_CS_team_size!(args)

Reconstruct COVID-secure workplace contacts according to new team size.

Only called if 'CS_team_size' is changed via intervention. CS workplaces must be active to implement change.
If this is the first change made to the CS workplace contact structure, store pre-intervention structure to enable later reversion.

Inputs:  `rng` - random number generator,
         `... parameter structures ...` \n
Outputs: None \n
Location: intervention_fns.jl
"""
function redefine_CS_team_size!(rng::MersenneTwister,
                                 configuration_options::ConfigurationOptions,
                                 network_params::NetworkParameters,
                                 contacts::ContactStructure)

   # Reconstruct work contacts to conform to capped team sizes

   @unpack cmax, workertypes, CS_active_flag = configuration_options
   @unpack CS_contacts_previously_generated, CS_contacts_changed_by_intervention = network_params

   if CS_active_flag

      #If previously generated and not stored, store now
      if (CS_contacts_previously_generated == true) &&
         (CS_contacts_changed_by_intervention == false)
         contacts.work_contacts_same_workplace_CS_store = deepcopy(contacts.work_contacts_same_workplace_CS)
         contacts.work_contacts_same_workplace_per_node_CS_store = deepcopy(contacts.work_contacts_same_workplace_per_node_CS)
         network_params.CS_contacts_changed_by_intervention = true
      end

      # Initialise CS related variables
      contacts.work_contacts_same_workplace_CS = Array{Array{Int64,1},1}(undef,cmax)
      contacts.work_contacts_same_workplace_per_node_CS = zeros(Int64,cmax)
      for node_itr = 1:cmax
         contacts.work_contacts_same_workplace_CS[node_itr] = Int64[]
      end

      # Construct COVID-secure work contacts
      CS_workplace_generation!(network_params.worker_nodes,
                                 network_params.nodes_by_workplace,
                                 contacts.work_contacts_same_workplace_CS,
                                 contacts.work_contacts_same_workplace_per_node_CS,
                                 network_params,
                                 rng)

   else
      error("Attempted to change CS team size without CS active")
   end

   return nothing
end

"""
   redefine_returned_to_work!(args)

Redefine each node's returned to work status, according to new 'workpercent'.

Only called if 'workpercent' is changed via intervention. Who returns and who does not relies on
the 'returned_to_work_hierarchy', which remains constant during a simulation, but is randomly changed between.

Inputs:  `... parameter structures ...` \n
Outputs: None \n
Location: intervention_fns.jl
"""
function redefine_returned_to_work!(configuration_options::ConfigurationOptions,
                                    worker_nodes::Array{WorkerParameters,1},
                                    workplace_generation_params::WorkplaceGenerationParameters,
                                    network_params::NetworkParameters)

   @unpack cmax, workertypes = configuration_options
   @unpack workpercent = workplace_generation_params
   @unpack nodes_by_sector = network_params

   if length(workpercent) != workertypes
      error("Number of sectors specified in intervention is inconsistent.")
   end

   for sector_itr = 1:workertypes
        returned_to_work_hierarchy = workplace_generation_params.returned_to_work_hierarchy[sector_itr]
        n_workers = length(returned_to_work_hierarchy)
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
   redefine_CT_engagement!(args)

Redefine each node's contact-tracing engagement status, according to new 'CT_engagement' or 'adherence' parameter.

Only called if 'CT_engagement' of 'adherence' are changed via intervention. Who engages and who does not
relies on the 'adherence_hierarchy', which remains constant during a simulation, but is randomly changed between.

Inputs:  `... parameter structures ...` \n
Outputs: None \n
Location: intervention_fns.jl
"""
function redefine_CT_engagement!(configuration_options::ConfigurationOptions,
                                 contact_tracing_params::ContactTracingParameters,
                                 CT_vars::ContactTracingVariables,
                                 states::NodeStates,
                                 infection_params::InfectionParameters)

   @unpack cmax, workertypes = configuration_options
   @unpack adherence_hierarchy = states
   @unpack adherence = infection_params
   @unpack CT_engagement = contact_tracing_params

   if configuration_options.contact_tracing_active == false
      error("Attempted to change CT_engagement whilst contact tracing inactive")
   end

   num_engaging = CT_engagement*adherence*cmax

   for node_itr = 1:cmax
      # For each worker, check if they engage with contact tracing
      if adherence_hierarchy[node_itr] <= num_engaging # engage with contact tracing
          CT_vars.Engage_with_CT[node_itr] = true
      else # do not engage with contact tracing
          CT_vars.Engage_with_CT[node_itr] = false
      end
   end

   return nothing
end

"""
   redefine_adherence!(args)

Redefine each node's adherence status, according to new 'adherence' parameter.

Only called if 'adherence' is changed via intervention. Who adheres and who does not relies on the
'adherence_hierarchy', which remains constant during a simulation, but is randomly changed between.

Inputs:  `cmax` - number of nodes,
         `... parameter structures ...` \n
Outputs: None \n
Location: intervention_fns.jl
"""
function redefine_adherence!(cmax::Int64,
                              states::NodeStates,
                              infection_params::InfectionParameters)

   @unpack adherence_hierarchy = states
   @unpack adherence = infection_params

   num_adhering =  adherence*cmax

   for node_itr = 1:cmax
      if adherence_hierarchy[node_itr] <= num_adhering # those who adhere will isolate when they get symptoms
         states.hh_isolation[node_itr] = 1 # adherence to household isolation = 1 if adherent, 0 if not.
      else
         states.hh_isolation[node_itr] = 0
      end
   end

   return nothing
end

"""
   redefine_daily_social_contacts!(args)

Regenerate daily dynamic social contacts from starttime to endtime, according to new settings, and store old settings.

Only called if one or more of 'network_generation_method_dynamic_social', 'social_workday_dd', 'social_nonworkday_dd',
'group_limit', 'dynamic_time_frame', 'n_groups_per_day_distribution', 'n_households_per_group' are changed via intervention.

Inputs:  `rng` - random number generator,
         `starttime` - time new contact structure is implemented,
         `intervention_set_itr` - number of current intervention set,
         `... parameter structures ...`,
         `daily_social_contact_setting_names` - array of setting names related to social contact structure,
         `old_daily_social_contact_settings` - array of pre-intervention setting values \n
Outputs: None \n
Location: intervention_fns.jl
"""
function redefine_daily_social_contacts!(rng::MersenneTwister,
                                          starttime::Int64,
                                          intervention_set_itr::Int64,
                                          configuration_options::ConfigurationOptions,
                                          network_params::NetworkParameters,
                                          contacts::ContactStructure,
                                          contact_component_store::Array{PreinterventionComponentStore,1},
                                          daily_social_contact_setting_names::Array{Symbol,1},
                                          old_daily_social_contact_settings::Array{Any,1})

   @unpack endtime = configuration_options

   store_length = length(contact_component_store)
   stored_times = zeros(Int64, store_length)
   n_relevant_settings = length(daily_social_contact_setting_names)

   # Check if contact structure from new settings has already been stored
   new_settings_already_stored = false
   store_id = -1
   for store_itr = 1:store_length
      current_store = contact_component_store[store_itr]
      same_setting_count = 0
      for setting_itr = 1:n_relevant_settings
         setting_name = daily_social_contact_setting_names[setting_itr]
         current_setting = getfield(network_params, setting_name)
         store_setting = getfield(current_store, setting_name)
         if current_setting == store_setting
            same_setting_count += 1
         end
      end
      # If all settings are the same, they have been stored before
      if same_setting_count == n_relevant_settings
         new_settings_already_stored = true
         store_id = store_itr
         break
      end
   end

   # Check if contact structure from current (old) settings has been stored
   old_settings_already_stored = false
   for store_itr = 1:store_length
      current_store = contact_component_store[store_itr]
      same_setting_count = 0
      for setting_itr = 1:n_relevant_settings
         setting_name = daily_social_contact_setting_names[setting_itr]
         old_setting = old_daily_social_contact_settings[setting_itr]
         store_setting = getfield(current_store, setting_name)
         if old_setting == store_setting
            same_setting_count += 1
         end
      end
      # If all settings are the same, they have been stored before
      if same_setting_count == n_relevant_settings
         old_settings_already_stored = true
         break
      end
   end

   # If current (old) contact structure has not been stored before, store it now
   if old_settings_already_stored == false
      preintervention_component = PreinterventionComponentStore(workday_social_contacts_by_day = contacts.workday_social_contacts_by_day[:, :],
                                                               nonworkday_social_contacts_by_day = contacts.nonworkday_social_contacts_by_day[:, :],
                                                               contacts_changed_by_intervention_time = starttime)
      for setting_itr = 1:n_relevant_settings
         setting_name = daily_social_contact_setting_names[setting_itr]
         old_setting = old_daily_social_contact_settings[setting_itr]
         setfield!(preintervention_component, setting_name, old_setting)
      end

      push!(contact_component_store, preintervention_component)
   end

   # If new contact structure has not been stored before, generate and store it now
   if new_settings_already_stored == false
      generate_social_contacts_each_day!(rng,
                                       starttime,
                                        network_params,
                                        contacts,
                                        configuration_options,
                                        intervention_set_itr)

       new_component = PreinterventionComponentStore(workday_social_contacts_by_day = contacts.workday_social_contacts_by_day[:, :],
                                                    nonworkday_social_contacts_by_day = contacts.nonworkday_social_contacts_by_day[:, :],
                                                    contacts_changed_by_intervention_time = starttime)
       for setting_itr = 1:n_relevant_settings
          setting_name = daily_social_contact_setting_names[setting_itr]
          new_setting = getfield(network_params, setting_name)
          setfield!(new_component, setting_name, new_setting)
       end

       push!(contact_component_store, new_component)
   # Otherwise load it from store
    else
      stored_component = contact_component_store[store_id]
      contacts.workday_social_contacts_by_day[starttime:endtime, :] = stored_component.workday_social_contacts_by_day[starttime:endtime, :]
      contacts.nonworkday_social_contacts_by_day[starttime:endtime, :] = stored_component.nonworkday_social_contacts_by_day[starttime:endtime, :]
   end

   return nothing
end

"""
   redefine_workplace_open!(args)

Redefine each workplace's open status, according to new sector-dependent 'sector_open' status.

Only called if 'sector_open' is changed via intervention.

Inputs: `network_params` - NetworkParameters structure \n
Outputs: None \n
Location: intervention_fns.jl
"""
function redefine_workplace_open!(network_params::NetworkParameters)

   @unpack sector_open = network_params

   n_sectors = length(sector_open)

   for sector_itr = 1:n_sectors
        n_workplaces = length(network_params.workplace_info[sector_itr])
        for workplace_itr = 1:n_workplaces
            network_params.workplace_info[sector_itr][workplace_itr].workplace_open = sector_open[sector_itr]
        end
    end

   return nothing
end

"""
   redefine_symp_isolation!(args)

Adjust NodeStates object tracking symptomatic isolation ('symp_isolation_array') for those currently isolating, according to new isolation time.

Only called if 'symp_isoltime' is changed via intervention. Only affects nodes currently in isolation.

Inputs: `time` - current time in simulation,
         `old_symp_isoltime` - value of symp_isoltime before intervention,
         `... parameter structures ...` \n
Outputs: None \n
Location: intervention_fns.jl
"""
function redefine_symp_isolation!(time::Int64,
                                    old_symp_isoltime::Int64,
                                    configuration_options::ConfigurationOptions,
                                    states::NodeStates,
                                    infection_params::InfectionParameters,
                                    CT_vars::ContactTracingVariables)

   @unpack cmax, endtime = configuration_options
   @unpack symp_isoltime = infection_params

   if configuration_options.isolation < 1
      error("Attempted to change symp_isoltime with isolation inactive")
   end

   array_time_now = time + 1

   change_symp_isoltime = symp_isoltime - old_symp_isoltime

   # If come to the end of latent time, move to infectious time etc
    for node_itr = 1:cmax

      # If in symptomatic isolation
      if states.symp_isolation_array[node_itr, array_time_now] == 1
         # If symp_isoltime has decreased
         if change_symp_isoltime < 0
            # Look back to check number of consecutive days in symp isolation so far (+ any delay to adherence)
            n_symp_isol_days = 1 + states.delay_adherence[node_itr]
            for day_itr = 1:symp_isoltime
               if states.symp_isolation_array[node_itr, array_time_now - day_itr] == 1
                  n_symp_isol_days += 1
               else
                  break
               end
            end
            # If not affected by false negative
            if CT_vars.Test_result_false_negative[node_itr] == false
               # Look forward and remove any excess symp isolation days
               for array_time_idx = array_time_now:(endtime+1)
                  if n_symp_isol_days >= symp_isoltime
                     states.symp_isolation_array[node_itr, array_time_idx] = 0
                  else
                     n_symp_isol_days += 1
                  end
               end
            # Otherwise affected by false negative
            else
               CT_delay = CT_vars.CT_delay_until_test_result[node_itr]
               # Look forward and remove any excess symp isolation days
               for array_time_idx = array_time_now:(endtime+1)
                  if n_symp_isol_days >= min(symp_isoltime, CT_delay)
                     states.symp_isolation_array[node_itr, array_time_idx] = 0
                  else
                     n_symp_isol_days += 1
                  end
               end
            end

         # Otherwise, add extra days on
         else
            # If not affected by false negative
            if CT_vars.Test_result_false_negative[node_itr] == false
               n_extra_days = change_symp_isoltime
               for array_time_idx = array_time_now:(endtime+1)
                  if n_extra_days > 0
                     if states.symp_isolation_array[node_itr, array_time_idx] == 0
                        states.symp_isolation_array[node_itr, array_time_idx] = 1
                        n_extra_days -= 1
                     end
                  end
               end
            # Otherwise is affected by false negative
            else
               # Only change if previous symp_isoltime was less than time to test result
               if (symp_isoltime - change_symp_isoltime) < CT_vars.CT_delay_until_test_result[node_itr]
                  n_extra_days = CT_vars.CT_delay_until_test_result[node_itr] - (symp_isoltime - change_symp_isoltime)
                  for array_time_idx = array_time_now:(endtime+1)
                     if n_extra_days > 0
                        if states.symp_isolation_array[node_itr, array_time_idx] == 0
                           states.symp_isolation_array[node_itr, array_time_idx] = 1
                           n_extra_days -= 1
                        end
                     end
                  end
               end
            end
         end
      end
   end

   return nothing
end

"""
   redefine_household_isolation!(args)

Adjust NodeStates object tracking household isolation ('hh_in_isolation_array') for those currently isolating, according to new isolation time.

Only called if 'household_isoltime' is changed via intervention. Only affects nodes currently in isolation.

Inputs: `time` - current time in simulation,
         `old_household_isoltime` - value of household_isoltime before intervention,
         `... parameter structures ...` \n
Outputs: None \n
Location: intervention_fns.jl
"""
function redefine_household_isolation!(time::Int64,
                                    old_household_isoltime::Int64,
                                    configuration_options::ConfigurationOptions,
                                    states::NodeStates,
                                    infection_params::InfectionParameters,
                                    CT_vars::ContactTracingVariables)

   @unpack cmax, endtime = configuration_options
   @unpack household_isoltime = infection_params

   array_time_now = time + 1

   change_household_isoltime = household_isoltime - old_household_isoltime

   # If come to the end of latent time, move to infectious time etc
    for node_itr = 1:cmax

      # If in household isolation
      if states.hh_in_isolation_array[node_itr, array_time_now] == 1
         # If household_isoltime has decreased
         if change_household_isoltime < 0
            # Look back to check number of consecutive days in hh isolation so far (+ any delay to adherence)
            n_hh_isol_days = 1 + states.delay_adherence[node_itr]
            for day_itr = 1:household_isoltime
               if states.hh_in_isolation_array[node_itr, array_time_now - day_itr] == 1
                  n_hh_isol_days += 1
               else
                  break
               end
            end
            # If not affected by false negative
            if CT_vars.Test_result_false_negative[node_itr] == false
               # Look forward and remove any excess hh isolation days
               for array_time_idx = array_time_now:(endtime+1)
                  if n_hh_isol_days >= household_isoltime
                     states.hh_in_isolation_array[node_itr, array_time_idx] = 0
                  else
                     n_hh_isol_days += 1
                  end
               end
            # Otherwise affected by false negative
            else
               CT_delay = CT_vars.CT_delay_until_test_result[node_itr]
               # Look forward and remove any excess hh isolation days
               for array_time_idx = array_time_now:(endtime+1)
                  if n_hh_isol_days >= min(household_isoltime, CT_delay)
                     states.hh_in_isolation_array[node_itr, array_time_idx] = 0
                  else
                     n_hh_isol_days += 1
                  end
               end
            end

         # Otherwise, add extra days on
         else
            # If not affected by false negative
            if CT_vars.Test_result_false_negative[node_itr] == false
               n_extra_days = change_household_isoltime
               for array_time_idx = array_time_now:(endtime+1)
                  if n_extra_days > 0
                     if states.hh_in_isolation_array[node_itr, array_time_idx] == 0
                        states.hh_in_isolation_array[node_itr, array_time_idx] = 1
                        n_extra_days -= 1
                     end
                  end
               end
            # Otherwise is affected by false negative
            else
               # Only change if previous household_isoltime was less than time to test result
               if (household_isoltime - change_household_isoltime) < CT_vars.CT_delay_until_test_result[node_itr]
                  n_extra_days = CT_vars.CT_delay_until_test_result[node_itr] - (household_isoltime - change_household_isoltime)
                  for array_time_idx = array_time_now:(endtime+1)
                     if n_extra_days > 0
                        if states.hh_in_isolation_array[node_itr, array_time_idx] == 0
                           states.hh_in_isolation_array[node_itr, array_time_idx] = 1
                           n_extra_days -= 1
                        end
                     end
                  end
               end
            end
         end
      end
   end

   return nothing
end

"""
   redefine_delay_adherence!(args)

Randomly redraw each node's delay to adherence, according to new distribution, and carry through changes to contract-tracing variables if appropriate.

Only called if 'delay_adherence_pmf' is changed via intervention.

Inputs: `rng` - random number generator,
         `... parameter structures ...` \n
Outputs: None \n
Location: intervention_fns.jl
"""
function redefine_delay_adherence!(rng::MersenneTwister,
                           configuration_options::ConfigurationOptions,
                           states::NodeStates,
                           infection_params::InfectionParameters,
                           contact_tracing_params::ContactTracingParameters,
                           CT_vars::ContactTracingVariables)

   @unpack cmax = configuration_options
   @unpack CT_days_before_symptom_included = contact_tracing_params

   csum_delay_adherence = cumsum(infection_params.delay_adherence_pmf)

   for node_itr = 1:cmax
      states.delay_adherence[node_itr] = draw_sample_from_pmf(csum_delay_adherence,
                                                              rng;
                                                              idx_offset = 1)

      if CT_vars.cmax > 0
         CT_vars.relevant_prev_days_for_CT[node_itr] = min(CT_days_before_symptom_included + states.delay_adherence[node_itr],
                                                      CT_days_before_symptom_included + 7)
      end
   end

   return nothing
end

"""
   redefine_CT_delay_until_test_result!(args)

Randomly redraw each node's delay to test result, according to new distribution.

Only called if 'CT_delay_until_test_result_pmf' is changed via intervention and contact tracing is active.

Inputs: `rng` - random number generator,
         `... parameter structures ...` \n
Outputs: None \n
Location: intervention_fns.jl
"""
function redefine_CT_delay_until_test_result!(rng::MersenneTwister,
                                             contact_tracing_params::ContactTracingParameters,
                                             CT_vars::ContactTracingVariables)

   if configuration_options.contact_tracing_active == false
      error("Attempted to change CT_delay_until_test_result whilst contact tracing inactive")
   end

   csum_test_result_delay = cumsum(contact_tracing_params.CT_delay_until_test_result_pmf)

   for node_itr = 1:cmax
       CT_vars.CT_delay_until_test_result[node_itr] = draw_sample_from_pmf(csum_test_result_delay,
                                                                 rng;
                                                                 idx_offset = 1)
    end

   return nothing
end

"""
   redefine_relevant_prev_days_for_CT!(args)

Redefine number of days to recall contacts for each node.

Only called if 'CT_days_before_symptom_included' is changed via intervention and contact tracing is active.

Inputs: `... parameter structures ...` \n
Outputs: None \n
Location: intervention_fns.jl
"""
function redefine_relevant_prev_days_for_CT!(configuration_options::ConfigurationOptions,
                                             CT_vars::ContactTracingVariables,
                                             contact_tracing_params::ContactTracingParameters,
                                             states::NodeStates)

   @unpack cmax = configuration_options
   @unpack CT_days_before_symptom_included = contact_tracing_params
   @unpack delay_adherence = states

   if configuration_options.contact_tracing_active == false
      error("Attempted to change CT_days_before_symptom_included whilst contact tracing inactive")
   end

   for node_itr = 1:cmax
      CT_vars.relevant_prev_days_for_CT[node_itr] = min(CT_days_before_symptom_included + delay_adherence[node_itr],
                                        CT_days_before_symptom_included + 7)
   end

   return nothing
end

"""
   redefine_CT_isolation!(args)

Adjust NodeStates object tracking isolation due to contact tracing ('CT_isolation_array') for those currently isolating, according to new isolation time.

Only called if 'CT_caused_isol_limit' is changed via intervention and contact tracing is active.
Only affects nodes currently in isolation.

Inputs: `... parameter structures ...` \n
Outputs: None \n
Location: intervention_fns.jl
"""
function redefine_CT_isolation!(configuration_options::ConfigurationOptions,
                                 contact_tracing_params::ContactTracingParameters,
                                 CT_vars::ContactTracingVariables,
                                 states::NodeStates)

   @unpack CT_caused_isol_limit = contact_tracing_params
   @unpack cmax = configuration_options
   @unpack CT_delay_until_test_result = CT_vars

   if configuration_options.contact_tracing_active == false
      error("Attempted to change CT_caused_isol_limit whilst contact tracing inactive")
   end

   array_time_now = time + 1

   change_CT_isoltime = CT_caused_isol_limit - old_CT_caused_isol_limit

   for node_itr = 1:cmax
      # Change isol time if currently isolating due to CT
      if states.CT_isolation_array[node_itr, array_time_now] == 1
         # If isol time has decreased
         if change_CT_isoltime < 0
            # Look back to check number of consecutive days in isolation so far (+ any delay to test result)
            n_CT_isol_days = 1 + CT_delay_until_test_result[node_itr]
            for day_itr = 1:CT_caused_isol_limit
               if states.CT_isolation_array[node_itr, array_time_now - day_itr] == 1
                  n_CT_isol_days += 1
               else
                  break
               end
            end
            # Look forward and remove any excess isolation days
            for array_time_idx = array_time_now:(endtime+1)
               if n_CT_isol_days >= CT_caused_isol_limit
                  states.CT_isolation_array[node_itr, array_time_idx] = 0
               else
                  n_CT_isol_days += 1
               end
            end
         # Otherwise if isol time has increased
         elseif change_CT_isoltime > 0
            n_extra_days = change_CT_isoltime
            for array_time_idx = array_time_now:(endtime+1)
               if n_extra_days > 0
                  if states.CT_isolation_array[node_itr, array_time_idx] == 0
                     states.CT_isolation_array[node_itr, array_time_idx] = 1
                     n_extra_days -= 1
                  end
               end
            end
         end
      end
   end

   return nothing
end

"""
   redefine_workplace_CT_threshold!(args)

Redefine number of infections required to cause workplace closure, for each workplace, according to new threshold proportion.

Only called if 'workplace_CT_threshold' is changed via intervention and workplace closures are active.

Inputs: `... parameter structures ...` \n
Outputs: None \n
Location: intervention_fns.jl
"""
function redefine_workplace_CT_threshold!(configuration_options::ConfigurationOptions,
                                       workplace_closure_params::WorkplaceClosureParameters,
                                       network_params::NetworkParameters)

   @unpack workertypes = configuration_options
   @unpack num_workplaces, workplace_CT_threshold = workplace_closure_params
   @unpack workplace_sizes = network_params

   if isassigned(workplace_closure_params.workplace_thresholds) == false
      error("Attempted to change workplace_CT_threshold before initialising workplace closures")
   elseif configuration_options.workplace_closure_active == false
      error("Attempted to change workplace_CT_threshold whilst workplace closures inactive")
   end

   for worktypeID = 1:workertypes
         n_workplaces = num_workplaces[worktypeID]
         for work_ID = 1:n_workplaces
             workplace_closure_params.workplace_thresholds[worktypeID][work_ID] = ceil(Int64,workplace_sizes[worktypeID][work_ID]*workplace_CT_threshold)
         end
     end

   return nothing
end

"""
   redefine_workplace_CT_memory!(args)

Redefine WorkplaceClosureParameters object tracking the number of infections in each workplace over a specified timeframe.

Only called if 'workplace_CT_memory' is changed via intervention and workplace closures are active.

Inputs: `... parameter structures ...` \n
Outputs: None \n
Location: intervention_fns.jl
"""
function redefine_workplace_CT_memory!(configuration_options::ConfigurationOptions,
                                       workplace_closure_params::WorkplaceClosureParameters,
                                       old_workplace_CT_memory::Int64)

   @unpack workertypes = configuration_options
   @unpack num_workplaces, workplace_CT_memory, workplace_memory = workplace_closure_params

   if isassigned(workplace_closure_params.workplace_thresholds) == false
      error("Attempted to change workplace_CT_threshold before initialising workplace closures")
   elseif configuration_options.workplace_closure_active == false
      error("Attempted to change workplace_CT_threshold whilst workplace closures inactive")
   end

   for worktypeID = 1:workertypes
      n_workplaces = num_workplaces[worktypeID]
      new_workplace_memory = zeros(n_workplaces, workplace_CT_memory)
      # using previous memory length, start at mod(time - 1) and track back for new memory length
      for time_itr = 1:workplace_CT_memory
         old_WP_time_slot = mod1(time - time_itr, old_workplace_CT_memory)
         new_WP_time_slot = mod1(time - time_itr, workplace_CT_memory)
         if time_itr <= old_workplace_CT_memory
            for workplace_itr = 1:n_workplaces
               new_workplace_memory[workplace_itr, new_WP_time_slot] = workplace_memory[worktypeID][workplace_itr, old_WP_time_slot]
            end
         else
            new_workplace_memory[:, new_WP_time_slot] .= 0
         end
      end
      workplace_closure_params.workplace_memory[worktypeID] = new_workplace_memory
   end

   return nothing
end

"""
   redefine_isolation_arrays!(args)

Redefine NodeStates objects tracking symptomatic and household isolation ('symp_isolation_array','hh_in_isolation_array'), due to isolation becoming active / inactive.

Only called if 'isolation' is changed via intervention.
If isolation is switched on, relevant individuals enter isolation.
If isolation is switched off, anyone isolating is released.

Inputs: `time` - current time in simulation,
         `... parameter structures ...` \n
Outputs: None \n
Location: intervention_fns.jl
"""
function redefine_isolation_arrays!(time::Int64,
                                    configuration_options::ConfigurationOptions,
                                    states::NodeStates,
                                    infection_params::InfectionParameters,
                                    contacts::ContactStructure)

   @unpack cmax, isolation, endtime = configuration_options
   @unpack symp_isoltime, household_isoltime = infection_params
   @unpack household_contacts_per_node, household_contacts = contacts

    # Check symptomatic status of all individuals
    # Iterate over each symptomatic.
    # - If they adhere, they isolate
    # - Check household members. If they adhere (and not symptomatic themselves),
    #       they enter household isolation tracker.
    # If isolation switched on, isolate currently symptomatic individuals
    if isolation == true
         for node_itr = 1:cmax
            # Check symptomatic status of all individuals
            # Two conditions: Not asymptomatic AND is in symptomatic phase of infection (if symptoms displayed)
            if (states.asymp[node_itr] == 0) && (states.timesymp[node_itr] > 0)
               # For each symptomatic.
               # - If the person adheres to new guidance, they now enter symptomatic isolation
               # - Check household members. If they adhere (and not symptomatic themselves),
               #       they enter household isolation tracker.

               # Index case isolation check
               if states.hh_isolation[node_itr] == 1
                   # Isolation due to presence of symptoms
                   # Set tracker to match time elapsed since symptom onset
                   time_left_in_symp_isol = (symp_isoltime - states.timesymp[node_itr])
                        # No +1, as increment on states.timesymp[node_itr] occurs after this loop
                   if time_left_in_symp_isol > 0
                       for time_itr = 1:time_left_in_symp_isol
                           array_time_idx = time + time_itr
                           if array_time_idx <= (endtime + 1)
                              states.symp_isolation_array[node_itr, array_time_idx] = 1
                              states.hh_in_isolation_array[node_itr, array_time_idx] = 0
                              states.CT_isolation_array[node_itr, array_time_idx] = 0
                           end
                       end
                   end
               end

                # Get time since index case first displayed symptoms
                time_left_in_hh_isol = household_isoltime - states.timesymp[node_itr]

                # Irrespective of whether index case self-isolates,
                # adherent members of their household may also isolate.
                if time_left_in_hh_isol > 0
                    for hh = 1:household_contacts_per_node[node_itr]
                        contact_ID = household_contacts[node_itr][hh]
                        if ((states.hh_isolation[contact_ID]==1) &&
                             (states.symp_isolation_array[contact_ID, time+1]==0)) # Household member not already symptomatic themselves

                            # Household member enters household isolation
                            # Populate household isolation time tracker array
                            for time_itr = 1:time_left_in_hh_isol
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
         end
    # Otherwise, release currently isolating individuals
    else
      for node_itr = 1:cmax
         if states.symp_isolation_array[node_itr, time+1] == 1
            for time_itr = 1:symp_isoltime
               array_time_idx = time + time_itr
               states.symp_isolation_array[node_itr, array_time_idx] = 0
            end
         end
         if states.hh_in_isolation_array[node_itr, time+1] == 1
            for time_itr = 1:household_isoltime
               array_time_idx = time + time_itr
               states.hh_in_isolation_array[node_itr, array_time_idx] = 0
            end
         end
      end
    end

   return nothing
end

"""
   load_interventions(intervention_name::String)

Load intervention variables according to input intervention name.

Returns an array{array} of InterventionVariables objects, defining sets of interventions (multiple interventions per replicate) to be implemented for each configuration.

Inputs: `intervention_name` -  name of intervention sets to load in \n
Outputs: `all_intervention_sets` - array{array} of InterventionVariables, defining sets of interventions to test \n
Location: intervention_fns.jl
"""
function load_interventions(intervention_name::String)

    if intervention_name=="none"
        all_intervention_sets = Array{Array{InterventionVariables,1},1}(undef,0)
    elseif intervention_name=="test_intervention"
        intervention1 = InterventionVariables(start_time = [30],
                                            workpercent = [zeros(WORKERTYPES)])
        intervention2 = InterventionVariables(start_time = [130],
                                            workpercent = [ones(WORKERTYPES)])
        all_intervention_sets = [[InterventionVariables(start_time = [30])],[intervention1, intervention2]]

    elseif intervention_name=="lockdown_duration_intervention"
        duration_intervention_options = [20:20:80;]
        n_int_sets = length(duration_intervention_options)
        all_intervention_sets = [[InterventionVariables(start_time = [30],
                                            workpercent = [zeros(WORKERTYPES)]),
                                InterventionVariables(start_time = [30+duration_intervention_options[i]],
                                            workpercent = [ones(WORKERTYPES)])] for i = 1:n_int_sets]

    elseif intervention_name=="synchronised_changedays_intervention"
        toff_intervention_options = [-1:4;]
        n_int_sets = length(toff_intervention_options)
        all_intervention_sets = [[InterventionVariables(start_time = [15],
                                                sameday = [0],
                                                toff = [toff_intervention_options[i]],
                                                ton = [0],
                                                contact_tracing_active = [true],
                                                isolation = [1])] for i = 1:n_int_sets]

    elseif intervention_name=="variable_changedays_intervention"
        toff_intervention_options = [-1:4;]
        n_int_sets = length(toff_intervention_options)
        all_intervention_sets = [[InterventionVariables(start_time = [15],
                                                sameday = [2],
                                                toff = [toff_intervention_options[i]],
                                                ton = [0],
                                                contact_tracing_active = [true],
                                                isolation = [1])] for i = 1:n_int_sets]

    elseif intervention_name=="weeks_on_intervention"
        toff_intervention_options = [1 1 2 2 2 2]
        ton_intervention_options = [1 1 2 2 1 1]
        sameday_intervention_options = [3 4 3 4 3 4]
        n_int_sets = length(toff_intervention_options)
        all_intervention_sets = [[InterventionVariables(start_time = [30],
                                                sameday = [sameday_intervention_options[i]],
                                                toff = [toff_intervention_options[i]],
                                                ton = [ton_intervention_options[i]])] for i = 1:n_int_sets]

    elseif intervention_name=="workpercent_intervention"
        workpercent_intervention_options = [1:-0.1:0;]
        n_int_sets = length(workpercent_intervention_options)
        all_intervention_sets = [[InterventionVariables(start_time = [50],
                                                workpercent = [workpercent_intervention_options[i]*ones(WORKERTYPES)],
                                                contact_tracing_active = [true],
                                                isolation = [1])] for i = 1:n_int_sets]
        # Add scenario where proportion of each type of worker returning to work is not constant across all sectors
        workpercent_config_segment2 = [0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.5, 0.8, 0.8, 0.5, 0.5, 0.8, 0.8, 0.8, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.8, 0.3, 0.3, 0.3, 0.7, 0.5, 0.3, 0.3, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.8, 0.7, 0.3, 0.5, 0.8, 0.8, 0.3]
        intervention_config2 = InterventionVariables(start_time = [15],
                                                workpercent = [workpercent_config_segment2],
                                                contact_tracing_active = [true],
                                                isolation = [1])
        # Concatenate the two sets of configs
        push!(all_intervention_sets,[intervention_config2])

    elseif intervention_name=="amount_backwards_CT_intervention"
        prob_backwards_CT_intervention_options = [0:0.1:1.0;]
        n_int_sets = length(prob_backwards_CT_intervention_options)
        all_intervention_sets = [[InterventionVariables(start_time = [15],
                                                        isolation = [1],
                                                        contact_tracing_active = [true],
                                                        prob_backwards_CT = [prob_backwards_CT_intervention_options[i]],
                                                        perform_CT_from_infector = [true])] for i = 1:n_int_sets]

    elseif intervention_name=="adherence_intervention"
        adherence_intervention_options = [0:0.1:1;]
        n_int_sets = length(adherence_intervention_options)
        all_intervention_sets = [[InterventionVariables(start_time = [15],
                                                adherence = [adherence_intervention_options[i]],
                                                contact_tracing_active = [true],
                                                isolation = [1])]
                                                for i = 1:n_int_sets]

    elseif intervention_name=="CS_intervention"
        CS_team_size_options = [2, 5, 10]
        CS_scale_transrisk_options = [0.25,0.5,0.75,1]
        n_CS_team_size_ops = length(CS_team_size_options)
        n_CS_scale_transrisk_ops = length(CS_scale_transrisk_options)
        n_int_sets = n_CS_team_size_ops*n_CS_scale_transrisk_ops
        all_intervention_sets = Array{Array{InterventionVariables,1},1}(undef, n_int_sets)
        set_idx = 1

        for CS_team_size_itr = 1:n_CS_team_size_ops
            CS_team_size_interv = CS_team_size_options[CS_team_size_itr]
            for CS_scale_transrisk_itr = 1:n_CS_scale_transrisk_ops
                CS_scale = CS_scale_transrisk_options[CS_scale_transrisk_itr]

                # Set up intervention list
                all_intervention_sets[set_idx] = [InterventionVariables(start_time = [15],
                                                        CS_scale_transrisk = [CS_scale*ones(Float64, WORKERTYPES)],
                                                        CS_team_size = [CS_team_size_interv],
                                                        contact_tracing_active = [true],
                                                        isolation = [1],
                                                        CS_active_flag = [true])]

                # Increment set index
                set_idx += 1
            end
        end

    elseif intervention_name=="CS_intervention_no_isol"
        CS_team_size_options = [2, 5, 10]
        CS_scale_transrisk_options = [0.25,0.5,0.75,1]
        n_CS_team_size_ops = length(CS_team_size_options)
        n_CS_scale_transrisk_ops = length(CS_scale_transrisk_options)
        n_int_sets = n_CS_team_size_ops*n_CS_scale_transrisk_ops
        all_intervention_sets = Array{Array{InterventionVariables,1},1}(undef, n_int_sets)
        set_idx = 1

        for CS_team_size_itr = 1:n_CS_team_size_ops
            CS_team_size_interv = CS_team_size_options[CS_team_size_itr]
            for CS_scale_transrisk_itr = 1:n_CS_scale_transrisk_ops
                CS_scale = CS_scale_transrisk_options[CS_scale_transrisk_itr]

                # Set up intervention list
                all_intervention_sets[set_idx] = [InterventionVariables(start_time = [15],
                                                      CS_active_flag = [true],
                                                    CS_scale_transrisk = [CS_scale*ones(Float64, WORKERTYPES)],
                                                    CS_team_size = [CS_team_size_interv])]

                # Increment set index
                set_idx += 1
            end
        end

    elseif intervention_name=="workplace_CT_threshold_intervention"
        workplace_CT_threshold_intervention_options = [0:0.05:0.4;]
        n_int_sets = length(workplace_CT_threshold_intervention_options)
        all_intervention_sets = [[InterventionVariables(start_time = [30],
                                                        contact_tracing_active = [true],
                                                        isolation = [1],
                                                        workplace_CT_threshold = [workplace_CT_threshold_intervention_options[i]],
                                                        workplace_closure_active = [true])] for i = 1:n_int_sets]

    elseif intervention_name=="hardfast_vs_weakslow_intervention"
        workpercent_intervention_options = [0.2:0.2:0.8;]
        intervention_length_options = [40:40:160;]
        n_int_sets = length(workpercent_intervention_options)
        all_intervention_sets = [[InterventionVariables(start_time = [30],
                                                workpercent = [workpercent_intervention_options[i]*ones(WORKERTYPES)]),
                                    InterventionVariables(start_time = [30+intervention_length_options[i]],
                                                        workpercent = [ones(WORKERTYPES)])] for i = 1:n_int_sets]

    elseif intervention_name=="adaptive_workplace_closures"
        max_lockdown_weeks = 4
        workpercent_intervention_options = [0, 0.5, 0.75]
        n_blocks_intervention_options = [1, 2, 4]
        block_gap_weeks_options = [1, 2, 4]          # weeks between lockdowns
        n_int_sets = length(workpercent_intervention_options)*length(n_blocks_intervention_options)*length(block_gap_weeks_options)
        all_intervention_sets = Array{Array{InterventionVariables,1},1}(undef, n_int_sets)
        workpercent_op_count = 1
        n_block_count = 1
        block_gap_count = 1
        for int_set_itr = 1:n_int_sets
            workpercent = workpercent_intervention_options[workpercent_op_count]
            n_blocks = n_blocks_intervention_options[n_block_count]
            block_gap = block_gap_weeks_options[block_gap_count]*7

            block_length = Int64((max_lockdown_weeks/(1-workpercent))/n_blocks)*7

            intervention_set = [InterventionVariables(start_time = [30],
                                                workpercent = [workpercent*ones(WORKERTYPES)]),
                                InterventionVariables(start_time = [30+block_length],
                                                workpercent = [ones(WORKERTYPES)])]

            for block_itr = 1:(n_blocks-1)
                push!(intervention_set, InterventionVariables(start_time = [30+block_itr*(block_length+block_gap)],
                                                    workpercent = [workpercent*ones(WORKERTYPES)]))
                push!(intervention_set, InterventionVariables(start_time = [30+block_itr*(block_length+block_gap)+block_length],
                                                    workpercent = [ones(WORKERTYPES)]))
            end
            all_intervention_sets[int_set_itr] = intervention_set

            block_gap_count += 1
            if block_gap_count > length(block_gap_weeks_options)
                block_gap_count = 1
                n_block_count += 1
            end
            if n_block_count > length(n_blocks_intervention_options)
                n_block_count = 1
                workpercent_op_count += 1
            end
        end

    elseif intervention_name=="multiphase_lockdown"
        workpercent_intervention_options = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
        start_intervention_options = [30, 60, 90, 120]
        length_intervention_options = [1, 2, 4, 8, 12, 16]
        n_blocks_intervention_options = [1, 2, 4, 6]
        block_gap_weeks_options = [1, 2, 4, 8, 12]          # weeks between lockdowns
        n_int_sets = length(workpercent_intervention_options)*length(start_intervention_options)*length(length_intervention_options)*length(n_blocks_intervention_options)*length(block_gap_weeks_options)
        all_intervention_sets = Array{Array{InterventionVariables,1},1}(undef, n_int_sets)
        workpercent_op_count = 1
        start_op_count = 1
        length_op_count = 1
        n_block_count = 1
        block_gap_count = 1
        for int_set_itr = 1:n_int_sets
            workpercent = workpercent_intervention_options[workpercent_op_count]
            start_time = start_intervention_options[start_op_count]
            total_length = length_intervention_options[length_op_count]
            n_blocks = n_blocks_intervention_options[n_block_count]
            block_gap = block_gap_weeks_options[block_gap_count]*7

            block_length = floor(Int64, (total_length*7)/n_blocks)

            intervention_set = [InterventionVariables(start_time = [start_time],
                                                workpercent = [workpercent*ones(WORKERTYPES)],
                                                scaling_social = [workpercent],
                                                scaling_random = [workpercent]),
                                InterventionVariables(start_time = [start_time+block_length],
                                                workpercent = [ones(WORKERTYPES)])]

            for block_itr = 1:(n_blocks-1)
                push!(intervention_set, InterventionVariables(start_time = [start_time+block_itr*(block_length+block_gap)],
                                                    workpercent = [workpercent*ones(WORKERTYPES)],
                                                    scaling_social = [workpercent],
                                                    scaling_random = [workpercent]))
                push!(intervention_set, InterventionVariables(start_time = [start_time+block_itr*(block_length+block_gap)+block_length],
                                                    workpercent = [ones(WORKERTYPES)]))
            end
            all_intervention_sets[int_set_itr] = intervention_set

            block_gap_count += 1
            if block_gap_count > length(block_gap_weeks_options)
                block_gap_count = 1
                n_block_count += 1
            end
            if n_block_count > length(n_blocks_intervention_options)
                n_block_count = 1
                length_op_count += 1
            end
            if length_op_count > length(length_intervention_options)
                length_op_count = 1
                start_op_count += 1
            end
            if start_op_count > length(start_intervention_options)
                start_op_count = 1
                workpercent_op_count += 1
            end
        end

    elseif intervention_name=="multiphase_lockdown_contact_struct"
        workpercent_intervention_options = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
        start_intervention_options = [30, 60, 90, 120]
        length_intervention_options = [1, 2, 4, 8, 12, 16]
        n_blocks_intervention_options = [1, 2, 4, 6]
        block_gap_weeks_options = [1, 2, 4, 8, 12]          # weeks between lockdowns
        n_int_sets = length(workpercent_intervention_options)*length(start_intervention_options)*length(length_intervention_options)*length(n_blocks_intervention_options)*length(block_gap_weeks_options)
        all_intervention_sets = Array{Array{InterventionVariables,1},1}(undef, n_int_sets)
        workpercent_op_count = 1
        start_op_count = 1
        length_op_count = 1
        n_block_count = 1
        block_gap_count = 1
        for int_set_itr = 1:n_int_sets
            workpercent = workpercent_intervention_options[workpercent_op_count]
            start_time = start_intervention_options[start_op_count]
            total_length = length_intervention_options[length_op_count]
            n_blocks = n_blocks_intervention_options[n_block_count]
            block_gap = block_gap_weeks_options[block_gap_count]*7

            block_length = floor(Int64, (total_length*7)/n_blocks)

            intervention_set = [InterventionVariables(start_time = [start_time],
                                                workpercent = [workpercent*ones(WORKERTYPES)],
                                                social_contact_scaling = [workpercent],
                                                random_contact_scaling = [workpercent]),
                                InterventionVariables(start_time = [start_time+block_length],
                                                workpercent = [ones(WORKERTYPES)])]

            for block_itr = 1:(n_blocks-1)
                push!(intervention_set, InterventionVariables(start_time = [start_time+block_itr*(block_length+block_gap)],
                                                    workpercent = [workpercent*ones(WORKERTYPES)],
                                                    social_contact_scaling = [workpercent],
                                                    random_contact_scaling = [workpercent]))
                push!(intervention_set, InterventionVariables(start_time = [start_time+block_itr*(block_length+block_gap)+block_length],
                                                    workpercent = [ones(WORKERTYPES)]))
            end
            all_intervention_sets[int_set_itr] = intervention_set

            block_gap_count += 1
            if block_gap_count > length(block_gap_weeks_options)
                block_gap_count = 1
                n_block_count += 1
            end
            if n_block_count > length(n_blocks_intervention_options)
                n_block_count = 1
                length_op_count += 1
            end
            if length_op_count > length(length_intervention_options)
                length_op_count = 1
                start_op_count += 1
            end
            if start_op_count > length(start_intervention_options)
                start_op_count = 1
                workpercent_op_count += 1
            end
        end

    elseif intervention_name=="transrisk_intervention"
        scaling_intervention_options = [0:0.2:1;]
        n_int_sets = length(scaling_intervention_options)
        all_intervention_sets = [[InterventionVariables(start_time = [30],
                                                scaling_work_static = [scaling_intervention_options[i]],
                                                scaling_work_dynamic = [scaling_intervention_options[i]],
                                                scaling_social = [scaling_intervention_options[i]],
                                                scaling_random = [scaling_intervention_options[i]],
                                                reset_time = [80])] for i = 1:n_int_sets]

    elseif intervention_name=="CT_engagement_intervention"
        CT_engagement_intervention_options = [0:0.2:1;]
        n_int_sets = length(CT_engagement_intervention_options)
        all_intervention_sets = [[InterventionVariables(start_time = [30],
                                                        contact_tracing_active = [true],
                                                        isolation = [1],
                                                        CT_engagement = [CT_engagement_intervention_options[i]])]
                                                                                for i = 1:n_int_sets]

    elseif intervention_name=="contact_structure_intervention"
        contact_scaling_intervention_options = [0:0.2:1;]
        n_int_sets = length(contact_scaling_intervention_options)
        all_intervention_sets = [[InterventionVariables(start_time = [30],
                                                social_contact_scaling = [contact_scaling_intervention_options[i]],
                                                random_contact_scaling = [contact_scaling_intervention_options[i]],
                                                reset_time = [80])] for i = 1:n_int_sets]

    elseif intervention_name=="adaptive_lockdown_adherence"
        workpercent_intervention_options = [0, 0.5, 1.0]
        start_intervention_options = [30, 60, 90]
        length_intervention_options = [1, 2, 4, 8]
        n_blocks_intervention_options = [1, 2]
        block_gap_weeks_options = [1, 2, 4]          # weeks between lockdowns
        n_int_sets = length(workpercent_intervention_options)*length(start_intervention_options)*length(length_intervention_options)*length(n_blocks_intervention_options)*length(block_gap_weeks_options)
        all_intervention_sets = Array{Array{InterventionVariables,1},1}(undef, n_int_sets)

        adherence = 0.7
        adherence_period = 28   # time scale (days) at which adherence naturally decreases
        n_adherence_periods = ceil(Int64, ENDTIME/adherence_period)
        adherence_period_decrease = 0.05

        workpercent_op_count = 1
        start_op_count = 1
        length_op_count = 1
        n_block_count = 1
        block_gap_count = 1
        for int_set_itr = 1:n_int_sets
            workpercent = workpercent_intervention_options[workpercent_op_count]
            start_time = start_intervention_options[start_op_count]
            total_length = length_intervention_options[length_op_count]
            n_blocks = n_blocks_intervention_options[n_block_count]
            block_gap = block_gap_weeks_options[block_gap_count]*7

            block_length = floor(Int64, (total_length*7)/n_blocks)

            policy_change_times = zeros(Int64, n_blocks*2)
            lockdown_active_array = zeros(Int64, ENDTIME)
            for block_itr = 1:n_blocks
                policy_change_times[(2*(block_itr-1)+1)] = start_time+(block_itr-1)*(block_length+block_gap)
                policy_change_times[(2*block_itr)] = start_time+(block_itr-1)*(block_length+block_gap)+block_length
                lockdown_active_array[(start_time+(block_itr-1)*(block_length+block_gap)):(start_time+(block_itr-1)*(block_length+block_gap)+block_length)] .= 1
            end

            cumulative_lockdown_days = cumsum(lockdown_active_array)
            cumulative_lockdown_periods = floor.(Int64, cumulative_lockdown_days./adherence_period)

            adherence_over_time = ones(Float64, ENDTIME).*adherence
            for time_itr = 1:ENDTIME
                adherence_over_time[time_itr] = adherence*(1-adherence_period_decrease)^cumulative_lockdown_periods[time_itr]
            end

            intervention_set = InterventionVariables[]
            push!(intervention_set, InterventionVariables(start_time = [15],
                                                         isolation = [1],
                                                         contact_tracing_active = [true]))

            time_itr = adherence_period + 1
            while time_itr <= ENDTIME

                if (time_itr+adherence_period-1) <= ENDTIME
                    current_time_period = [time_itr:(time_itr+adherence_period-1);]
                else
                    current_time_period = [time_itr:ENDTIME;]
                end

                # Find policy changes in current time period
                policy_change_ids = findall((policy_change_times.>=current_time_period[1]).*(policy_change_times.<=current_time_period[end]).==1)
                n_policy_change = length(policy_change_ids)

                if adherence_over_time[time_itr] != adherence_over_time[time_itr-adherence_period]
                    push!(intervention_set, InterventionVariables(start_time = [time_itr],
                                                    adherence = [adherence_over_time[time_itr]]))
                end
                if n_policy_change > 0
                    for policy_change_itr = 1:n_policy_change
                        push!(intervention_set, InterventionVariables(start_time = [policy_change_times[policy_change_ids[policy_change_itr]]],
                                                        workpercent = [workpercent*ones(WORKERTYPES)],
                                                        social_contact_scaling = [workpercent],
                                                        random_contact_scaling = [workpercent]))
                    end

                end
                time_itr += adherence_period
            end

            all_intervention_sets[int_set_itr] = intervention_set

            block_gap_count += 1
            if block_gap_count > length(block_gap_weeks_options)
                block_gap_count = 1
                n_block_count += 1
            end
            if n_block_count > length(n_blocks_intervention_options)
                n_block_count = 1
                length_op_count += 1
            end
            if length_op_count > length(length_intervention_options)
                length_op_count = 1
                start_op_count += 1
            end
            if start_op_count > length(start_intervention_options)
                start_op_count = 1
                workpercent_op_count += 1
            end
        end
    elseif intervention_name == "social_easing_intervention"
        group_limit = [3:3:15;]
        n_households_per_group = [1:15;]
        all_group_limits = repeat(group_limit, inner = length(n_households_per_group))
        all_n_households_per_group = repeat(n_households_per_group, outer = length(group_limit))
        all_intervention_sets = [[InterventionVariables(start_time = [2],
                                                        group_limit = [all_group_limits[ii]],
                                                        n_households_per_group = [all_n_households_per_group[ii]],
                                                        network_generation_method_dynamic_social = ["household_groups"])]
                                    for ii = 1:length(all_group_limits)]
    elseif intervention_name == "social_easing_test"
        all_intervention_sets = [[InterventionVariables(start_time = [20],
                                                        group_limit = [100],
                                                        n_households_per_group = [50],
                                                        network_generation_method_dynamic_social = ["household_groups"])],
                                [InterventionVariables(start_time = [50],
                                                    group_limit = [100],
                                                    n_households_per_group = [50],
                                                    network_generation_method_dynamic_social = ["household_groups"])]]
    elseif intervention_name == "isolation_test"
        all_intervention_sets = [[InterventionVariables(start_time = [20],
                                                        isolation = [1],
                                                        workpercent = [0.2*ones(WORKERTYPES)])],
                                [InterventionVariables(start_time = [50],
                                                    isolation = [0],
                                                    workpercent = [0.2*ones(WORKERTYPES)])]]
    end

    n_intervention_sets = length(all_intervention_sets)
    for intervention_set_itr = 1:n_intervention_sets
        intervention_set = all_intervention_sets[intervention_set_itr]
        n_interventions = length(intervention_set)
        for intervention_itr = 1:n_interventions
            intervention = intervention_set[intervention_itr]
            if isassigned(intervention.start_time) == false
                error("No start time given for intervention $intervention_itr in set $intervention_set_itr")
            elseif length(intervention.start_time) > 1
                error("Multiple start times given for intervention $intervention_itr in set $intervention_set_itr")
            end
            if isassigned(intervention.reset_time) == true
                if length(intervention.reset_time) > 1
                    error("Multiple reset times given for intervention $intervention_itr in set $intervention_set_itr")
                end
            end
        end
    end

    return all_intervention_sets::Array{Array{InterventionVariables,1},1}
end

"""
   condition_close_example(args)

Example of triggered intervention condition function. Intervention is triggered when daily incidence
exceeds 1% of the population.

Inputs: `intervention_trigger_input_data` - InterventionDataFeeds structure storing infection data,
         `time` - current time in simulation,
         `configuration_options` - ConfigurationOptions structure (for cmax) \n
Outputs: `output_bool` - boolean indicating if condition has been met or not \n
Location: intervention_fns.jl
"""
function condition_close_example(intervention_trigger_input_data::InterventionDataFeeds,
                                        time::Int64,
                                        configuration_options::ConfigurationOptions)

   @unpack cmax = configuration_options
   @unpack rep_inf_this_timestep = intervention_trigger_input_data

   # Check if daily incidence condition surpassed
   daily_incidence = sum(rep_inf_this_timestep)/n_nodes
   if (daily_incidence > 0.01)
      output_bool = true
   else
      output_bool = false
   end

   return output_bool::Bool
end

"""
   condition_open_example(args)

Example of triggered intervention condition function. Intervention is triggered when daily incidence
falls below 0.1% of the population.

Inputs: `intervention_trigger_input_data` - InterventionDataFeeds structure storing infection data,
         `time` - current time in simulation,
         `configuration_options` - ConfigurationOptions structure (for cmax) \n
Outputs: `output_bool` - boolean indicating if condition has been met or not \n
Location: intervention_fns.jl
"""
function condition_open_example(intervention_trigger_input_data::InterventionDataFeeds,
                                        time::Int64,
                                        configuration_options::ConfigurationOptions)

   @unpack cmax = configuration_options
   @unpack rep_inf_this_timestep = intervention_trigger_input_data

   # Check if daily incidence condition declined to low enough level
   daily_incidence = sum(rep_inf_this_timestep)/n_nodes
   if daily_incidence < 0.001
      output_bool = true
   else
      output_bool = false
   end

   return output_bool::Bool
end

"""
   affect_close_example(network_params::NetworkParameters)

Example of triggered intervention affect function. When triggered, sectors 9 and 39 are closed.

Inputs:  `network_params` - NetworkParameters structure \n
Outputs: None \n
Location: intervention_fns.jl
"""
function affect_close_example!(network_params::NetworkParameters)

   @unpack workplace_info, sector_open = network_params

   # Get number of sectors in use
   n_sectors = length(workplace_info)

   # Throw an error if number of sectors is not compatible with this function
   if n_sectors != 41
      error("Number of sectors in use is $(n_sectors). Selected affect function requires use of 41 sectors.")
   end

   # The sectors to be closed
   close_sector_IDs = [9,39]
   n_sectors_closed = length(close_sector_IDs)

   # Iterate over each sector to be closed
   for close_sector_itr = 1:n_sectors_closed
      # Get sector of interest
      close_sector_ID = close_sector_IDs[close_sector_itr]

      if (sector_open[close_sector_ID] == true)
         # Get number of workplaces in that sector
         n_workplaces = length(workplace_info[close_sector_ID])

         # Set status of each workplace in that sector to be closed
         for workplace_itr = 1:n_workplaces
            workplace_info[close_sector_ID][workplace_itr].workplace_open = false
         end

         # Update sector_open field
         sector_open[close_sector_ID] = false
      end
   end
end

"""
   affect_open_example(network_params::NetworkParameters)

Example of triggered intervention affect function. When triggered, sectors 9 and 39 are opened.

Inputs:  `network_params` - NetworkParameters structure \n
Outputs: None \n
Location: intervention_fns.jl
"""
function affect_open_example!(network_params::NetworkParameters)

   @unpack workplace_info, sector_open = network_params

   # Get number of sectors in use
   n_sectors = length(workplace_info)

   # Throw an error if number of sectors is not compatible with this function
   if n_sectors != 41
      error("Number of sectors in use is $(n_sectors). Selected affect function requires use of 41 sectors.")
   end

   # The sectors to be opened
   open_sector_IDs = [9,39]
   n_sectors_open = length(open_sector_IDs)

   # Iterate over each sector to be closed
   for open_sector_itr = 1:n_sectors_open
      # Get sector of interest and number of workplaces in that sector
      open_sector_ID = open_sector_IDs[open_sector_itr]

      if (sector_open[open_sector_ID] == false)
         # Get number of workplaces in that sector
         n_workplaces = length(workplace_info[open_sector_ID])

         # Set status of each workplace in that sector to be closed
         for workplace_itr = 1:n_workplaces
            workplace_info[open_sector_ID][workplace_itr].workplace_open = true
         end

         # Update sector_open field
         sector_open[open_sector_ID] = true
      end
   end
end

"""
   check_triggered_interventions!(args)

Check if the state of the outbreak is sufficient to trigger the implementation of an intervention, and implement if so.

Inputs: `time` - current time in simulation,
         `... parameter structures ...` \n
Outputs: None \n
Location: intervention_fns.jl
"""
function check_triggered_interventions!(time::Int64,
                                        network_params::NetworkParameters,
                                        states::NodeStates,
                                        output::SimulationOutputs)

    # Package health outcome measures that may be used in decision
    # to enact an intervention
    intervention_trigger_input_data = InterventionDataFeeds(rep_inf_this_timestep = states.rep_inf_this_timestep,
                                                                numlat = output.numlat[output_time_idx,count,intervention_set_itr],
                                                                numinf = output.numinf[output_time_idx,count,intervention_set_itr],
                                                                numrep = output.numrep[output_time_idx,count,intervention_set_itr],
                                                                newinf = output.newinf[output_time_idx,count,intervention_set_itr])

    for interv_trig_itr = 1:n_intervention_fns
        # Check if condition is met
        condition_fn = intervention_fns[interv_trig_itr,1]
        eval_condition_fn = condition_fn(intervention_trigger_input_data,
                                          time,
                                          network_params)
        if eval_condition_fn == true
            # If condition satisfied, apply the affect
            chosen_affect_fn = intervention_fns[interv_trig_itr,2]
            chosen_affect_fn(network_params)
        end
    end
end
