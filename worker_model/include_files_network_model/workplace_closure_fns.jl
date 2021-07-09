"""
File containing functions related to workplace closures, triggered by a specified proportion
of infections within a workplace
"""

"""
    initialise_workplace_closure_params!(args)

Initialises objects within the WorkplaceClosureParameters structure, previously defaulted to empty arrays.

Called during initialisation phase if workplace closures are active in the configuration,
or within intervention functions if activated during the outbreak.

Inputs: `... parameter structures ...` \n
Outputs: None \n
Location: workplace_closure_fns.jl
"""
function initialise_workplace_closure_params!(workplace_closure_params::WorkplaceClosureParameters,
                                                configuration_options::ConfigurationOptions,
                                                network_params::NetworkParameters)

   @unpack workertypes = configuration_options
   @unpack workplace_sizes = network_params

   workplace_closure_params.workplace_memory = Array{Array{Int64,2},1}(undef, workertypes) # initialise the memory for each workplace
   workplace_closure_params.work_closed_time = Array{Array{Int64,1},1}(undef, workertypes) # initialise the timer for work closure
   workplace_closure_params.workplace_thresholds = Array{Array{Int64,1},1}(undef, workertypes) # initialise the threshold for work closure
   workplace_closure_params.num_workplaces = zeros(Int64, workertypes)

   for worktypeID = 1:workertypes
      n_workplaces = length(workplace_sizes[worktypeID])
      workplace_closure_params.num_workplaces[worktypeID] = n_workplaces
      workplace_closure_params.work_closed_time[worktypeID] = zeros(Int64, n_workplaces)
      workplace_closure_params.workplace_memory[worktypeID] = zeros(Int64, n_workplaces,workplace_closure_params.workplace_CT_memory)
      workplace_closure_params.workplace_thresholds[worktypeID] = zeros(Int64, n_workplaces)
      for work_ID = 1:n_workplaces
          workplace_closure_params.workplace_thresholds[worktypeID][work_ID] = ceil(Int64, workplace_sizes[worktypeID][work_ID]*workplace_closure_params.workplace_CT_threshold)
      end
   end

   return nothing
end

"""
    reinitialise_workplace_closure_params!(args)

Reinitialises objects within the WorkplaceClosureParameters structure, previously initialised.

Called during reinitialisation phase of each replicate if workplace closures parameters have been changed via intervention.

Inputs: `... parameter structures ...` \n
Outputs: None \n
Location: workplace_closure_fns.jl
"""
function reinitialise_workplace_closure_params!(workplace_closure_params::WorkplaceClosureParameters,
                                                configuration_options::ConfigurationOptions,
                                                network_params::NetworkParameters)

   @unpack workertypes = configuration_options
   @unpack workplace_sizes = network_params

   for worktypeID = 1:workertypes
      n_workplaces = length(workplace_sizes[worktypeID])
      workplace_closure_params.num_workplaces[worktypeID] = n_workplaces
      workplace_closure_params.work_closed_time[worktypeID] = zeros(Int64, n_workplaces)
      workplace_closure_params.workplace_memory[worktypeID] = zeros(Int64, n_workplaces, workplace_closure_params.workplace_CT_memory)
      workplace_closure_params.workplace_thresholds[worktypeID] = zeros(Int64, n_workplaces)
      for work_ID = 1:n_workplaces
          workplace_closure_params.workplace_thresholds[worktypeID][work_ID] = ceil(Int64, workplace_sizes[worktypeID][work_ID]*workplace_closure_params.workplace_CT_threshold)
      end
   end

   return nothing
end

"""
    reactive_workplace_closures!(args)

If workplace closures are active, check if each workplace has enough infections to trigger closure,
or has been closed long enough to be re-opened.

Inputs: `... parameter structures ...` \n
Outputs: None \n
Location: workplace_closure_fns.jl
"""
function reactive_workplace_closures!(configuration_options::ConfigurationOptions,
                                    network_params::NetworkParameters,
                                    workplace_closure_params::WorkplaceClosureParameters)

    @unpack workertypes = configuration_options
    @unpack workplace_CT_memory, time_WC, workplace_memory,
            num_workplaces, workplace_thresholds = workplace_closure_params

    for worktypeID = 1:workertypes
        for work_ID = 1:num_workplaces[worktypeID]
            if workplace_closure_params.work_closed_time[worktypeID][work_ID]>0
                # if the workplace is closed, move on the time counter
                workplace_closure_params.work_closed_time[worktypeID][work_ID] += 1
            else
                # otherwise decide if the workplace should be closed
                total_infections = 0
                for ii=1:workplace_CT_memory
                    total_infections += workplace_memory[worktypeID][work_ID,ii]
                end

                # Number of infections in time window exceeds threshold
                # Set workplace to be closed
                if total_infections>workplace_thresholds[worktypeID][work_ID]
                    workplace_closure_params.work_closed_time[worktypeID][work_ID] = 1
                    network_params.workplace_info[worktypeID][work_ID].workplace_open = false
                end
            end

            if workplace_closure_params.work_closed_time[worktypeID][work_ID]>time_WC
                # if the workplace has been closed long enough, open it
                workplace_closure_params.work_closed_time[worktypeID][work_ID] = 0

                # Check sector is not closed. If sector is open, workplace set to be open
                if (network_params.sector_open[worktypeID] == true)
                    network_params.workplace_info[worktypeID][work_ID].workplace_open = true
                end
            end
        end
    end
end
