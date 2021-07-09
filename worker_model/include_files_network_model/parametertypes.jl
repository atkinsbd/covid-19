"""
File containing the parameter structures to be used throughout the model
"""

"""
   ConfigurationOptions()

Structure containing broad configuration options, such as which interventions are active,
default worker pattern, number of nodes and length of simulation.

Location: parametertypes.jl
"""
@with_kw mutable struct ConfigurationOptions

   # Set if contact tracing is active or not (Bool type variable)
    contact_tracing_active::Bool = false
    # Set if workplace closures is active or not (Bool type variable)
    workplace_closure_active::Bool = false
    # set if backwards contact tracing is active or not (Bool type variable)
    perform_CT_from_infector::Bool = false
    # Set if isolation due to symptoms and household isolation are active (0 = inactive)
    isolation::Int64 = 0
    # Set if Covid-Secure workplaces are active (Bool type variable)
    CS_active_flag::Bool = false
    # Set proportion of recovered individuals at start of simulation
    recov_propn::Float64 = 0.

    # Set default working pattern (see populate_atwork! function)
    sameday::Int64 = 3
    ton::Int64 = 1
    toff::Int64 = 0

    # Variables to be set from input arguments
    rng_seed::Int64 = 0
    cmax::Int64 = 0
    endtime::Int64 = 0
    countfinal::Int64 = 0
    workertypes::Int64 = 0

    # Set method of seeding initial infection states
    seed_initial_states_fn::Function = seed_states_with_uncertainty
end

"""
   WorkerParameters()

Structure allocated to each node, containing information on workplace and household membership,
working-from-home vs returned-to-work status and transmission risk in each setting.

Location: parametertypes.jl
"""
@with_kw mutable struct WorkerParameters
   returned_to_work::Int64 = 0            # If 1, worker returns to workplace as designated by atwork schedule.
   sector_ID::Int64                       # The sector the worker has been allocated to
   workplace_ID::Int64                    # The ID of the workplace (within that sector)
   household_ID::Int64 = 0                # Household the individual has been assigned to
   transrisk_household::Float64 = 0.      # Secondary attack rate within household setting
   transrisk_static_work::Float64 = 0.    # Secondary attack rate within static workplace setting
   transrisk_dynamic_work::Float64 = 0.   # Secondary attack rate within dynamic workplace setting
   transrisk_social::Float64 = 0.         # Secondary attack rate within social setting
   transrisk_random::Float64 = 0.         # Secondary attack rate within random setting
end

"""
   WorkplaceParameters()

Structure assigned to each workplace, containing information on whether the workplace is open
and if COVID-secure measures are in place or not.

Location: parametertypes.jl
"""
@with_kw mutable struct WorkplaceParameters
   covid_secure::Bool = false    # If true, workplace is covid secure. Can alter transmission risk & contacts
   workplace_open::Bool = true   # For imposing workplace closures, can set flag to false
end

"""
   NetworkParameters()

Structure containing parameters required to generate the network.

Location: parametertypes.jl
"""
@with_kw mutable struct NetworkParameters

   # Method used to generate contacts (can be "configuration" or "ER")
   network_generation_method::String = "configuration"

   # Method used to generate dynamic social contacts (see generate_social_contacts_each_day! function)
   network_generation_method_dynamic_social::String = "cluster"

   # Parameters used within "ER" network generation method
   prob_anyworker_contact::Float64 = 0.0                                 # Contact probability for workers in different sectors
   prob_social_contact::Float64 = 0.0                                    # Contact probability for social contacts with any other node (used to form static friend group)
   prob_random_contact::Float64 = 0.0                                    # Contact probability for random contacts with any other node
   prob_workertype_contact::Array{Float64,1} = Array{Float64,1}(undef,0) # Sector specific contact probability for workers in same sector but different workplace
   dd_within_workplace::Array{Float64,1} = Array{Float64,1}(undef,0)     # Mean degree for workplaces in each sector, used to generate contacts within workplaces
   dynamic_conts_mean::Array{Float64,1} = Array{Float64,1}(undef,0)      # Mean number of dynamic contacts for each sector, used to generate normally distributed number of dynamic contacts
   dynamic_conts_sd::Array{Float64,1} = Array{Float64,1}(undef,0)        # Standard deviation of dynamic contacts for each sector, used to generate normally distributed number of dynamic contacts

   # Parameters used within "configuration" network generation method
   workplace_degree_distribution::Array{Distribution,1} = repeat([Distributions.LogNormal(1.896,1.233)], WORKERTYPES)        # Sector specific degree distribution for contacts at work
   between_workplace_contact_probs::Array{Float64,1} = [0.05]                                           # Probability that each contact made is with other workplace compared to within workplace
   max_contacts_work_dynamic::Int64 = 100                                                               # Maximum number of dynamic work contacts for an individual per day
   social_group_size_distribution::Distribution = Distributions.LogNormal(3.14,1.41)                    # Distribution of friendship group sizes
   max_contacts_social::Int64 = 100                                                                     # Maximum contacts allowed for friendship group
   friend_of_friend_prob::Float64 = 0.5                                                                 # Probability of making contacts with friends-of-friends opposed to others
   workplace_dynamic_degree_distribution::Array{Distribution,1} = repeat([Distributions.LogNormal(1.262,1.315)], WORKERTYPES)# Sector-specific degree distribution of dynamic contacts in workplaces

   # Parameters used within "configuration", "cluster" and "repeated" dynamic social network generation method
   social_workday_dd::Distribution = Distributions.LogNormal(1.397,1.27)                          # Distribution of number of social contacts made per workday
   social_nonworkday_dd::Distribution = Distributions.LogNormal(1.536,1.153)                      # Distribution of number of social contacts made per nonworkday
   cluster_social_contacts::Bool = true                                                           # Whether or not to cluster daily social contacts ("configuration" method only)

   # Parameters used within "ER" dynamic social network generation method
   n_social_mean_workday::Int64 = 1                                      # Average number of social contacts on a workday, used to generate Poisson distributed number of social contacts each day
   n_social_mean_nonworkday::Int64 = 5                                   # Average number of social contacts on a nonworkday, used to generate Poisson distributed number of social contacts each day

   # Parameters used in "repeated", "fixed_daily_degree", "groups", "household_groups" dynamic social network generation methods
   group_limit::Int64 = 100                                                  # Maximum group size for single meeting
   dynamic_time_frame::Int64 = 1                                             # Number of days to repeat same social contacts
   n_groups_per_day_distribution::Distribution = Distributions.Poisson(1)    # Distribution of number of social groups each individual meets per day ("groups" only)
   n_households_per_group::Int64 = group_limit                               # Number of households making up each social group meeting ("household_groups" only)

   # Household size distribution (1, 2, 3, 4, 5, 6+)
   household_size_distribution::Array{Float64,1} = [0.314, 0.423, 0.075, 0.025, 0.006, 0.002]/sum([0.314, 0.423, 0.075, 0.025, 0.006, 0.002])

   # Flag to indicate if sector is open (true) or closed (false)
   sector_open::Array{Bool,1} = Array{Bool,1}(undef,0)

   # Scaling on number of contacts for lockdowns - acts as probability of 'activating' an already generated contact
   social_contact_scaling::Float64 = 1.0
   random_contact_scaling::Float64 = 1.0

   # COVID-secure workplace related parameters
   CS_team_size::Int64 = 10                          # Size of fully-connected teams
   CS_contacts_previously_generated::Bool = false    # Flag recording if CS contacts have been generated already
   CS_contacts_changed_by_intervention::Bool = false # Flag recording if CS contacts have been generated and changed via intervention

   # Workplace / worker related storage
   worker_nodes::Array{WorkerParameters,1} = Array{WorkerParameters,1}(undef,0)                     # Worker information array. Entry per worker.
   workplace_sizes::Array{Array{Int64,1},1} = Array{Int64,1}[]                                      # Object storing sizes of workplaces, by sector and workplace ID.
   workplace_info::Array{Array{WorkplaceParameters,1},1} = Array{WorkplaceParameters,1}[]           # Object storing open and covid-secure status of each workplace
   nodes_by_workplace::Array{Array{Array{Int64,1},1},1} = Array{Array{Array{Int64,1},1},1}(undef,0) # Object storing node IDs by sector and workplace membership
   nodes_by_sector::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}(undef,0)                      # Object storing node IDs by sector membership only
end

"""
   InfectionParameters()

Structure containing parameters relating to the infection process.

Location: parametertypes.jl
"""
@with_kw mutable struct InfectionParameters

   # Parameters defining epidemiological process
   d_incub::Distribution = Erlang(6,0.88)                                         # Distribution of length of incubation period
   inftime::Int64 = 4                                                             # time to spend in infectious, pre-symptomatic state
   symptime::Int64 = 10                                                           # time to spend in symptomatic state
   dist_infectivity::Array{Float64,1} =  [0.0369, 0.0491, 0.0835, 0.1190, 0.1439,
                                          0.1497, 0.1354, 0.1076, 0.0757, 0.0476,
                                          0.0269, 0.0138, 0.0064, 0.0044]         # Distribution of infectiousness (Corrected He (4 day pre-symptomatic period) & 10 day symptomatic period)
   delay_household_infection_pmf::Array{Float64,1} = dist_infectivity./sum(dist_infectivity) # Distribution of delay in household infection between household members (first entry = 0 days)
   probasymp_dist::Uniform{Float64} = Uniform(0.5,0.8)   # Distribution of probability of being asymptomatic
   probasymp::Float64 = -1.                              # Probability of being asymptomatic (redrawn every replicate)

   # Parameters defining risk of transmission
   transrisk_household_group_mean::Array{Float64,1} = [0.48,0.40,0.33,0.22]  # Mean transmission risk in households of different sizes
   transrisk_household_group_sd::Array{Float64,1} = [0.06,0.06,0.05,0.05]    # SD transmission risk in households of different sizes
   transrisk_household_mean_average::Float64 = 0.37                          # Mean transmission risk in households averaged across all sizes
   transrisk_household_sd_average::Float64 = 0.03                            # SD transmission risk in households averaged across all sizes

   transrisk_static_work_mean::Array{Float64,1} = (transrisk_household_mean_average/0.5).*
                                                   [0.13, 0.13, 0.13, 0.13, 0.13, 0.13,
                                                   0.21, 0.21, 0.21, 0.30, 0.30, 0.17,
                                                   0.27, 0.55, 0.17, 0.17, 0.17, 0.59,
                                                   0.13, 0.15, 0.18, 0.13, 0.17, 0.13,
                                                   0.19, 0.12, 0.17, 0.19, 0.31, 0.29,
                                                   0.29, 0.29, 0.18, 0.18, 0.05, 0.18,
                                                   0.19, 0.15, 0.22, 0.13, 0.17]  # Baseline risk in workplace, static and dynamic - sector dependent (currently only works for 41 sectors)

   transrisk_dynamic_work_mean::Array{Float64,1} = (transrisk_household_mean_average/0.5).*[0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.06, 0.06, 0.06, 0.33, 0.33, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.14, 0.21, 0.20, 0.14, 0.21, 0.25, 0.04, 0.35, 0.25, 0.25, 0.35, 0.33, 0.39, 0.39, 0.39, 0.36, 0.36, 0.25, 0.36, 0.35, 0.49, 0.56, 0.50, 0.25]
   transrisk_static_work_sd::Array{Float64,1} = [transrisk_static_work_mean[i]*(transrisk_household_sd_average/transrisk_household_mean_average) for i=1:length(transrisk_static_work_mean)]
   transrisk_dynamic_work_sd::Array{Float64,1} = [transrisk_dynamic_work_mean[i]*(transrisk_household_sd_average/transrisk_household_mean_average) for i=1:length(transrisk_dynamic_work_mean)]

   transrisk_social_mean::Float64 = (transrisk_household_mean_average/0.5)*0.355   # Baseline risk for social transmission, set relative to household transmission
   transrisk_social_sd::Float64 = transrisk_social_mean*(transrisk_household_sd_average/transrisk_household_mean_average)

   transrisk_random_mean::Float64 = (transrisk_household_mean_average/0.5)*0.355   # Baseline risk for random transmission, set relative to household transmission
   transrisk_random_sd::Float64 = transrisk_random_mean*(transrisk_household_sd_average/transrisk_household_mean_average)

   # Parameters that scale the risk of transmission
   iso_trans_scaling::Float64 = 1.                                # Transmission scaling for those symptomatic (self-isolation)
   asymp_trans_scaling_dist::Uniform{Float64} = Uniform(0.3,0.7)  # Distribution of transmission scaling for asymptomatics
   asymp_trans_scaling::Float64 = -1.                             # Transmission scaling for asymptomatics (redrawn each replicate)
   scaling_social::Float64 = 0.8                                  # Scale transmission risk in social setting
   scaling_work_static::Float64 = 0.8                             # Scale transmission risk in static work setting
   scaling_work_dynamic::Float64 = 0.8                            # Scale transmission risk in dynamic work setting
   scaling_household::Float64 = 0.8                               # Scale transmission risk in household setting
   scaling_random::Float64 = 0.8                                  # Scale transmission risk in random setting
   CS_scale_transrisk::Array{Float64,1} = ones(Float64, WORKERTYPES) # Transmission scaling within a COVID-secure workplace (set by configuration or intervention)

   # Parameters relating to enforced isolation
   symp_isoltime::Int64 = 10                          # Number of days spent in isolation caused by developing symptoms
   household_isoltime::Int64 = 10                     # Number of days spent in isolation caused by symptomatic household member
   adherence::Float64 = 0.7                           # Proportion of population who will adhere
   delay_adherence_pmf::Array{Float64,1} = [1.,0.,0.] # Distribution of delay in reporting symptoms (first entry = 0 days)

   # Functions used to assign transmission risks in different settings
   assign_household_transrisk_fn::Function = assign_household_transmit_household_size!
   assign_workplace_static_transrisk_fn::Function = assign_workplace_static_transmit!
   assign_workplace_dynamic_transrisk_fn::Function = assign_workplace_dynamic_transmit!
   assign_social_transrisk_fn::Function = assign_social_transmit!
   assign_random_transrisk_fn::Function = assign_random_transmit!
end

"""
   ContactTracingParameters()

Structure containing parameters relating to contact tracing.

Location: parametertypes.jl
"""
@with_kw mutable struct ContactTracingParameters

   CT_engagement::Float64 = 1.  # Proportion of adhering individuals that will engage with CT (total engagement proportion = adherence*CT_engagement)

   CT_delay_until_test_result_pmf::Array{Float64,1} = [0., 0., 1.,] # Distribution of delay between testing and result (first entry = day of test)
   csum_test_result_delay = cumsum(CT_delay_until_test_result_pmf)  # Cumulative distribution of delay between test and result, used to draw random delay for each node each replicate

   CT_days_before_symptom_included::Int64 = 2  # Number of days before symptoms develop that CT will cover

   # Propotion of tests returning a positive if infected depending on days since infection (first entry = day after infection)
   test_detection_prob_vec::Array{Float64,1} = [0.,0.11,0.59,0.785,0.83,0.845,0.84,0.82,0.79,0.76,  # Days 1-10
                                                0.72,0.68,0.64,0.59,0.54,0.485,0.445,0.405,0.37,0.335, # Days 11-20
                                                0.30,0.27,0.24,0.22,0.20,0.18,0.16,0.15,0.14,0.13] # Days 21-30

   CT_caused_isol_limit::Int64 = 10  # Amount of time spent in isolation if contact traced

   dynamic_contacts_recalled_propn::Array{Float64,1} = [0.5,0.4,0.3,0.2,0.1]  # Proportion of dynamic work contacts remembered x days ago (first entry = 1 day before engagement with CT)
   social_contacts_recalled_propn::Array{Float64,1} = [1.,1.,1.,1.,1.]        # Proportion of social contacts remembered x days ago (first entry = 1 day before engagement with CT)

   prob_backwards_CT::Float64 = 0.  # Proportion of people that can identify their infector (set to 0 for no backwards CT)

   infector_engage_with_CT_prob::Float64 = 1.0  # For those not calready adhering, probability they engage with CT if idenfitied as possible infector

end

"""
   ContactTracingVariables()

Structure containing storage objects required for performing contact tracing.
Objects have zero length unless contact tracing is activated.

Location: parametertypes.jl
"""
@with_kw mutable struct ContactTracingVariables

   cmax::Int64 = 0 # sizes to initialise arrays (number of nodes)

   relevant_prev_days_for_CT::Array{Int64,1} = zeros(Int64,cmax) # The number of days prior to symptoms that each node remembers

   Inds_to_be_contacted::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}(undef,cmax) # array{array} of IDs storing those to be contacted in CT per node

   Time_to_test_result::Array{Int64,1} = -1*ones(Int64,cmax) # Variables for waiting for test results, set at -1 until activated

   Test_result_false_negative::Array{Bool,1} = Array{Bool,1}(undef,cmax) # Boolean vector to store whether a false negative test result would be returned, redrawn each replicate

   Engage_with_CT::Array{Bool,1} = Array{Bool,1}(undef,cmax) # Boolean vector to store if individual will engage with CT or not, redrawn each replicate

   CT_engagement_hierarchy::Array{Int64,1} = Array{Int64,1}(undef,cmax) # Ordered array of node IDs, to decide who engages with CT, reshuffled each replicate

   Recall_infector::Array{Int64,1} = zeros(Int64,cmax) # Array to keep track of whether an infected recalls their infector, redrawn each replicate

   CT_delay_until_test_result::Array{Int64,1} = zeros(Int64,cmax) # Delay before test result is returned for each node, redrawn each replicate
end

"""
   ContactStructure()

Structure containing storage objects for contact layers.

Location: parametertypes.jl
"""
@with_kw mutable struct ContactStructure

   # Work contact layer
   work_contacts_same_workplace::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}(undef,0)  # array{array} of work contacts at the same workplace
   work_contacts_same_workplace_per_node::Array{Int64,1} = zeros(Int64,0)                    # number of work contacts at same workplace per node

   work_contacts_other_workplace::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}(undef,0) # array{array} of work contacts at different workplace
   work_contacts_other_workplace_per_node::Array{Int64,1} = zeros(Int64,0)                   # number of work contacts at different workplace per node

   work_contacts_same_workplace_CS::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}(undef,0) # array{array} of work contacts at the same workplace under COVID-secure restrictions
   work_contacts_same_workplace_per_node_CS::Array{Int64,1} = zeros(Int64,0)                   # number of work contacts at same workplace per node under COVID-secure restrictions
   work_contacts_same_workplace_CS_store::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}(undef, 0) # storage for CS contact layer if changed by intervention
   work_contacts_same_workplace_per_node_CS_store::Array{Int64,1} = zeros(Int64, 0)                   # storage for CS contact layer if changed by intervention

   dynamic_worker_contacts::Array{Array{Int64,1},2} = Array{Array{Int64,1},2}(undef,0,0)  # 2D array{array} of dynamic work contacts made each day for each node

   # Social contact layer
   social_contacts::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}(undef,0) # array{array} of static social contacts, representing friendship groups of each node (sampled from each day to obtain actual social contacts)
   social_contacts_per_node::Array{Int64,1} = zeros(Int64,0)                   # number of static social contacts per node

   workday_social_contacts_by_day::Array{Array{Int64,1},2} = Array{Array{Int64,1},2}(undef,0,0)     # 2D array{array} of social contacts made each work-day for each node
   nonworkday_social_contacts_by_day::Array{Array{Int64,1},2} = Array{Array{Int64,1},2}(undef,0,0)  # 2D array{array} of social contacts made each nonwork-day for each node

   social_meeting_group_sizes::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}(undef,0)  # array{array} of meeting group sizes for each intervention set (used in generate_social_contacts_each_day_household_groups functions)

   # Household contact layer
   household_contacts::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}(undef,0)  # array{array} of contacts within same household
   household_contacts_per_node::Array{Int64,1} = zeros(Int64,0)                    # number of contacts within same household per node
   nodes_by_household::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}(undef,0)  # array{array} of node IDs arranged into household groups (used in generate_social_contacts_each_day_household_groups functions)

   # Random contact layer
   dynamic_random_contacts::Array{Array{Int64,1},2} = Array{Array{Int64,1},2}(undef,0,0) # 2D array{array} of dynamic random contacts made each day by each node

end

"""
   SimulationOutputs()

Structure containing storage objects for outputs to be saved from simulations.

Location: parametertypes.jl
"""
@with_kw mutable struct SimulationOutputs

   # Parameters required for dimensions of storage arrays
   n_intervention_sets::Int64 = 0
   cmax::Int64 = 0
   endtime::Int64 = 0
   countfinal::Int64 = 0

   # 2D outputs
   numlat::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # number latently infected
   numinf::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # number infectious
   numrep::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # number reporting infection
   prevlat::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # Prevalence for latently infected
   prevsymp::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # Prevalence for symptomatic infectious
   prevasymp::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # Prevalence for asymptomatic infectious
   prevpresymp::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # Prevalence of pre-symptomatic symptomatic infectious
   prevrec::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # Prevalence for recovereds
   newinf::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # number of new infections
   newasymp::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # number of newly infected asymptomatics
   workersinf::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # number of workers that are newly infected
   workersasymp::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # number of workers that are newly infected asymptomatics
   num_isolating::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # number of individuals isolating at each timepoint
   num_household_isolating::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # number of individuals isolating due to household members having symptoms at each timepoint
   num_symp_isolating::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # number of individuals isolating due to symptoms at each timepoint
   num_isolating_CTcause::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # number of individuals isolating as a contact of a positive test
   num_CT::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # total number of recallable contacts
   num_infected::Array{Int64,3} = zeros(Int64,cmax,countfinal,n_intervention_sets) # number of infections caused by that node
   dynamic_infection_count::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets) # number infections over dynamic links
   Rt::Array{Float64,3} = zeros(Float64,endtime+1,countfinal,n_intervention_sets) # real time R value (number of secondary infections caused by nodes that were newly infected on that day)
   transmission_setting::Array{Int64,4} = zeros(Int64, endtime+1, countfinal, n_intervention_sets, 5) # records proportion of infections ocurring in each setting: social, work-static, work-dynamic, household, other
   num_init_infected::Array{Array{Int64,2}} = Array{Array{Int64,2}}(undef,countfinal) # number of infections caused by the initially infected nodes
   social_meeting_groups::Array{Float64,2} = zeros(Float64, 8, n_intervention_sets+1) # summary of sizes of social meeting groups (used in generate_social_contacts_each_day_household_groups style contact generation only)

   # 2D outputs. Testing related
   tests_performed::Array{Int64,3} = zeros(Int64,endtime+1,countfinal,n_intervention_sets)

   # 3D outputs. Testing related
    # Slices: 1 - True positive; 2 - false negative; 3 - True negative; 4 - false positive.
   test_outcomes::Array{Int64,4} = zeros(Int64,endtime+1,countfinal,n_intervention_sets,4)

   # 1D outputs
   infected_by::Array{Int64,2} = zeros(Int64,cmax,n_intervention_sets)  # who each node was infected by
   var_num_infected::Array{Int64,2} = zeros(Int64,countfinal,n_intervention_sets) # variance in the number of infections caused by each node
   mean_init_generation_time::Array{Float64,2} = zeros(Float64,countfinal,n_intervention_sets) # mean initial generation time

end

"""
   InterventionDataFeeds()

Structure keeping track of variables used to trigger an intervention (triggered interventions not yet complete)

Location: parametertypes.jl
"""
@with_kw mutable struct InterventionDataFeeds
   rep_inf_this_timestep::Array{Int64,1} = Array{Float64,1}[] # Indicator of whether a node reports symptomatic infection during current timestep
   numlat::Int64 = 0                                          # Number of nodes currently in latent phase
   numinf::Int64 = 0                                          # Number of nodes currently infectious
   numrep::Int64 = 0                                          # Number of nodes currently reporting infection
   newinf::Int64 = 0                                          # Number of new infections
end

"""
   NodeStates()

Structure containing storage objects to track the status of each node, including working, adherence, isolation and infection status.
These objects are reinitialised every replicate (or timestep where stated)

Location: parametertypes.jl
"""
@with_kw mutable struct NodeStates

   # Parameters required for dimensionality of storage objects
   cmax::Int64 = 0
   endtime::Int64 = 0

   # Working status
   atwork::Array{Int64,2} = zeros(Int64,cmax,endtime) # whether the node is at work on each day

   # Isolation status
   hh_in_isolation_array::Array{Int64,2} = zeros(Int64,cmax,endtime+1) # Household isolation. Daily record for each node
   symp_isolation_array::Array{Int64,2} = zeros(Int64,cmax,endtime+1) # Symptomatic isolation. Daily record for each node
   CT_isolation_array::Array{Int64,2} = zeros(Int64,cmax,endtime+1) # Isolation by contact tracing. Daily record for each node
   daily_record_atworkplace::Array{Int64,2} = zeros(Int64, endtime, cmax) # Flag array to keep a log of whether each node is at workplace each timestep
   daily_record_inisol::Array{Int64,2} = zeros(Int64, endtime, cmax) # Flag array to keep a log of whether each node is in isolation each timestep

   # Adherence status
   hh_isolation::Array{Int64,1} = zeros(Int64,cmax) # Whether each node adheres to isolation guidance. (1) Yes. (0) No.
   delay_adherence::Array{Int64,1} = zeros(Int64,cmax) # Number of days between symptom onset and reporting for each node
   adherence_hierarchy::Array{Int64,1} = Array{Int64,1}(undef,0) # # Ordered array of node IDs, to decide who adheres, reshuffled each replicate

   # Infection status
   asymp::Array{Int64,1} = zeros(Int64,cmax) # whether the node will be asymptomatic
   timelat::Array{Int64,1} = zeros(Int64,cmax) # time currently spent in latent state
   timeinf::Array{Int64,1} = zeros(Int64,cmax) # time currently spent in infectious state
   timesymp::Array{Int64,1} = zeros(Int64,cmax) # time currently spent in symptomatic state
   lattime::Array{Int64,1} = zeros(Int64,cmax) # time to spend in latent state
   acquired_infection::Array{Int64,1} = zeros(Int64, cmax) # time node acquired infection
   time_to_symps::Array{Int64,1} = zeros(Int64, cmax) # Number of days between infection and symptoms, randomly drawn for each node
   infected_by::Array{Int64,1} = zeros(Int64, cmax) # ID of node that caused each node's infection
   rep_inf_this_timestep::Array{Int64,1} = zeros(Int64,cmax) # whether each node reports symptoms this timestep (reset each timestep)

end

"""
   WorkplaceGenerationParameters()

Structure containing parameters required for workplace generation.

Location: parametertypes.jl
"""
@with_kw mutable struct WorkplaceGenerationParameters

   workpercent::Array{Float64,1} = ones(Float64, WORKERTYPES)                        # percentage of each sector that is physically present at the workplace
   workforce_proportion::Array{Float64,1} = ones(Float64, WORKERTYPES)./WORKERTYPES  # proportion of nodes belonging to each sector
   workplace_size_mean::Array{Float64,1} = ones(Int64, WORKERTYPES).*10              # mean size of workplaces (used in workplace_size_gen_using_normal_dist only)
   workplace_size_sd::Array{Float64,1} = ones(Int64, WORKERTYPES).*10                # sd size of workplaces (used in workplace_size_gen_using_normal_dist only)

   workplace_size_gen_fn::Function = workplace_size_sampled_from_empirical_pmf       # function used to generate workplace sizes (workplace_size_gen_using_normal_dist or workplace_size_sampled_from_empirical_pmf)

   returned_to_work_hierarchy::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}(undef, WORKERTYPES) # Ordered array of node IDs, per sector, used to decide who is physically present at workplace, reshuffled each replicate
end

"""
   WorkplaceClosureParameters()

Structure containing parameters and storage objects required to perform workplace closures.
Storage objects have 0 length unless initialised by configuration or intervention.

Location: parametertypes.jl
"""
@with_kw mutable struct WorkplaceClosureParameters

   # Parameters
   workplace_CT_memory::Int64 = 7        # Number of days in the past that workplace infections are 'remembered'
   workplace_CT_threshold::Float64 = 0.5 # Proportion of infected employees required to trigger workplace closure
   time_WC::Int64 = 14                   # Number of days to close workplaces for

   # Storage objects
   num_workplaces::Array{Int64,1} = Array{Int64,1}(undef, 0)                         # Number of workplaces in each sector
   work_closed_time::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}(undef, 0)     # Number of days each workplace has been closed for (sector x workplace)
   workplace_thresholds::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}(undef, 0) # Number of infected employees required to trigger closure in each workplace (sector x workplace)
   workplace_memory::Array{Array{Int64,2},1} = Array{Array{Int64,1},1}(undef, 0)     # Number of infections on each remembered day for each workplace (sector x workplace x day)
   WP_memory_slot::Int64 = 0                                                         # Modulo - moves through memory object to update memory storage
end

"""
   ConfigurationVariables()

Structure containing all variables that can be changed within the configuration.
These changes are applied from the start of a simulation and remain the same for each replicate,
unless changed during a replicate by intervention. Each variable name corresponds to a
counterpart within other structures, with an extra leading dimension: e.g. an integer variable
in the NetworkParameters structure will be represented here as an array{Int64,1}. This allows
every variable to remain undefined. Only those variables that are defined here will be altered
during the configuration set up. For explanations of each variable, see the corresponding variable
in relevant parameter structure.

Location: parametertypes.jl
"""
@with_kw mutable struct ConfigurationVariables

   n_configs::Int64 = 0

   ton::Array{Int64,1} = Array{Int64,1}(undef,0)
   toff::Array{Int64,1} = Array{Int64,1}(undef,0)
   sameday::Array{Int64,1} = Array{Int64,1}(undef,0)
   prob_backwards_CT::Array{Float64,1} = Array{Float64,1}(undef,0)
   infector_engage_with_CT_prob::Array{Float64,1} = Array{Float64,1}(undef,0)
   CT_engagement::Array{Float64,1} = Array{Float64,1}(undef,0)
   adherence::Array{Float64,1} = Array{Float64,1}(undef,0)
   CS_team_size::Array{Int64,1} = Array{Int64,1}(undef,0)
   dynamic_time_frame::Array{Int64,1} = Array{Int64,1}(undef,0)
   group_limit::Array{Int64,1} = Array{Int64,1}(undef,0)
   max_contacts_social::Array{Int64,1} = Array{Int64,1}(undef,0)
   max_contacts_work_dynamic::Array{Int64,1} = Array{Int64,1}(undef,0)
   n_households_per_group::Array{Int64,1} = Array{Int64,1}(undef,0)
   scaling_social::Array{Float64,1} = Array{Float64,1}(undef,0)
   scaling_work_static::Array{Float64,1} = Array{Float64,1}(undef,0)
   scaling_work_dynamic::Array{Float64,1} = Array{Float64,1}(undef,0)
   scaling_household::Array{Float64,1} = Array{Float64,1}(undef,0)
   scaling_random::Array{Float64,1} = Array{Float64,1}(undef,0)
   friend_of_friend_prob::Array{Float64,1} = Array{Float64,1}(undef,0)
   isolation::Array{Int64,1} = Array{Int64,1}(undef,0)
   recov_propn::Array{Float64,1} = Array{Float64,1}(undef,0)
   prob_anyworker_contact::Array{Float64,1} = Array{Float64,1}(undef,0)
   prob_social_contact::Array{Float64,1} = Array{Float64,1}(undef,0)
   prob_random_contact::Array{Float64,1} = Array{Float64,1}(undef,0)
   prob_workertype_contact::Array{Float64,1} = Array{Float64,1}(undef,0)
   n_social_mean_workday::Array{Int64,1} = Array{Int64,1}(undef,0)
   n_social_mean_nonworkday::Array{Int64,1} = Array{Int64,1}(undef,0)
   symp_isoltime::Array{Int64,1} = Array{Int64,1}(undef,0)
   household_isoltime::Array{Int64,1} = Array{Int64,1}(undef,0)
   CT_days_before_symptom_included::Array{Int64,1} = Array{Int64,1}(undef,0)
   CT_caused_isol_limit::Array{Int64,1} = Array{Int64,1}(undef,0)
   inftime::Array{Int64,1} = Array{Int64,1}(undef,0)
   symptime::Array{Int64,1} = Array{Int64,1}(undef,0)
   workplace_CT_threshold::Array{Float64,1} = Array{Float64,1}(undef,0)
   workplace_CT_memory::Array{Int64,1} = Array{Int64,1}(undef,0)
   social_contact_scaling::Array{Float64,1} = Array{Float64,1}(undef,0)
   random_contact_scaling::Array{Float64,1} = Array{Float64,1}(undef,0)
   iso_trans_scaling::Array{Float64,1} = Array{Float64,1}(undef,0)
   time_WC::Array{Int64,1} = Array{Int64,1}(undef,0)

   transrisk_social_mean::Array{Float64,1} = Array{Float64,1}(undef,0)
   transrisk_social_sd::Array{Float64,1} = Array{Float64,1}(undef,0)
   transrisk_random_mean::Array{Float64,1} = Array{Float64,1}(undef,0)
   transrisk_random_sd::Array{Float64,1} = Array{Float64,1}(undef,0)

   transrisk_dynamic_work_mean::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
   transrisk_dynamic_work_sd::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
   transrisk_static_work_mean::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
   transrisk_static_work_sd::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
   transrisk_household_group_mean::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
   transrisk_household_group_sd::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)

   perform_CT_from_infector::Array{Bool,1} = Array{Bool,1}(undef,0)
   contact_tracing_active::Array{Bool,1} = Array{Bool,1}(undef,0)
   workplace_closure_active::Array{Bool,1} = Array{Bool,1}(undef,0)
   CS_active_flag::Array{Bool,1} = Array{Bool,1}(undef,0)
   cluster_social_contacts::Array{Bool,1} = Array{Bool,1}(undef,0)
   sector_open::Array{Array{Bool,1},1} = Array{Array{Bool,1},1}(undef,0)

   network_generation_method::Array{String,1} = Array{String,1}(undef,0)
   network_generation_method_dynamic_social::Array{String,1} = Array{String,1}(undef,0)

   seed_initial_states_fn::Array{Function,1} = Array{Function,1}(undef,0)
   workplace_size_gen_fn::Array{Function,1} = Array{Function,1}(undef,0)
   assign_household_transrisk_fn::Array{Function,1} = Array{Function,1}(undef,0)
   assign_workplace_static_transrisk_fn::Array{Function,1} = Array{Function,1}(undef,0)
   assign_workplace_dynamic_transrisk_fn::Array{Function,1} = Array{Function,1}(undef,0)
   assign_social_transrisk_fn::Array{Function,1} = Array{Function,1}(undef,0)
   assign_random_transrisk_fn::Array{Function,1} = Array{Function,1}(undef,0)

   social_workday_dd::Array{Distribution,1} = Array{Distribution,1}(undef,0)
   social_nonworkday_dd::Array{Distribution,1} = Array{Distribution,1}(undef,0)
   n_groups_per_day_distribution::Array{Distribution,1} = Array{Distribution,1}(undef,0)
   social_group_size_distribution::Array{Distribution,1} = Array{Distribution,1}(undef,0)
   asymp_trans_scaling_dist::Array{Distribution,1} = Array{Distribution,1}(undef,0)
   d_incub::Array{Distribution,1} = Array{Distribution,1}(undef,0)
   probasymp_dist::Array{Uniform{Float64},1} = Array{Uniform{Float64},1}(undef,0)

   workplace_degree_distribution::Array{Array{Distribution,1},1} = Array{Array{Distribution,1},1}(undef,0)
   workplace_dynamic_degree_distribution::Array{Array{Distribution,1},1} = Array{Array{Distribution,1},1}(undef,0)

   workpercent::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
   between_workplace_contact_probs::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
   household_size_distribution::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
   dd_within_workplace::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
   dynamic_conts_mean::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
   dynamic_conts_sd::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
   dist_infectivity::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
   delay_adherence_pmf::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
   CT_delay_until_test_result_pmf::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
   workforce_proportion::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
   workplace_size_mean::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
   workplace_size_sd::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
   CS_scale_transrisk::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
   test_detection_prob_vec::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
   dynamic_contacts_recalled_propn::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
   social_contacts_recalled_propn::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
end

"""
   InterventionVariables()

Structure containing all variables that can be changed via intervention during a replicate.
These changes are applied at 'start_time' for every replicate, then reverted at the end of each replicate,
or at 'reset_time' if defined.
Each variable name corresponds to a counterpart within other structures, with an extra leading dimension:
e.g. an integer variable in the NetworkParameters structure will be represented here as an array{Int64,1}.
This allows every variable to remain undefined. Only those variables that are defined here will be altered
via intervention. For explanations of each variable, see the corresponding variable in relevant parameter structure.

Location: parametertypes.jl
"""
@with_kw mutable struct InterventionVariables

   start_time::Array{Int64,1} = Array{Int64,1}(undef,0)  # Start time of intervention
   reset_time::Array{Int64,1} = Array{Int64,1}(undef,0)  # End time of each intervention (not required)

   sameday::Array{Int64,1} = Array{Int64,1}(undef,0)
   ton::Array{Int64,1} = Array{Int64,1}(undef,0)
   toff::Array{Int64,1} = Array{Int64,1}(undef,0)
   group_limit::Array{Int64,1} = Array{Int64,1}(undef,0)
   dynamic_time_frame::Array{Int64,1} = Array{Int64,1}(undef,0)
   symp_isoltime::Array{Int64,1} = Array{Int64,1}(undef,0)
   household_isoltime::Array{Int64,1} = Array{Int64,1}(undef,0)
   CT_days_before_symptom_included::Array{Int64,1} = Array{Int64,1}(undef,0)
   CT_caused_isol_limit::Array{Int64,1} = Array{Int64,1}(undef,0)
   workplace_CT_memory::Array{Int64,1} = Array{Int64,1}(undef,0)
   isolation::Array{Int64,1} = Array{Int64,1}(undef,0)
   time_WC::Array{Int64,1} = Array{Int64,1}(undef,0)
   CS_team_size::Array{Int64,1} = Array{Int64,1}(undef,0)
   n_households_per_group::Array{Int64,1} = Array{Int64,1}(undef,0)

   CT_engagement::Array{Float64,1} = Array{Float64,1}(undef,0)
   workplace_CT_threshold::Array{Float64,1} = Array{Float64,1}(undef,0)
   adherence::Array{Float64,1} = Array{Float64,1}(undef,0)
   social_contact_scaling::Array{Float64,1} = Array{Float64,1}(undef,0)
   random_contact_scaling::Array{Float64,1} = Array{Float64,1}(undef,0)
   iso_trans_scaling::Array{Float64,1} = Array{Float64,1}(undef,0)
   scaling_work_static::Array{Float64,1} = Array{Float64,1}(undef,0)
   scaling_work_dynamic::Array{Float64,1} = Array{Float64,1}(undef,0)
   scaling_social::Array{Float64,1} = Array{Float64,1}(undef,0)
   scaling_random::Array{Float64,1} = Array{Float64,1}(undef,0)
   infector_engage_with_CT_prob::Array{Float64,1} = Array{Float64,1}(undef,0)

   workpercent::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
   delay_adherence_pmf::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
   CT_delay_until_test_result_pmf::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
   test_detection_prob_vec::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)
   CS_scale_transrisk::Array{Array{Float64,1},1} = Array{Array{Float64,1},1}(undef,0)

   sector_open::Array{Array{Bool,1},1} = Array{Array{Bool,1},1}(undef,0)

   CS_active_flag::Array{Bool,1} = Array{Bool,1}(undef,0)
   contact_tracing_active::Array{Bool,1} = Array{Bool,1}(undef,0)
   workplace_closure_active::Array{Bool,1} = Array{Bool,1}(undef,0)
   perform_CT_from_infector::Array{Bool,1} = Array{Bool,1}(undef,0)
   prob_backwards_CT::Array{Float64,1} = Array{Float64,1}(undef,0)

   network_generation_method_dynamic_social::Array{String,1} = Array{String,1}(undef,0)

   social_workday_dd::Array{Distribution,1} = Array{Distribution,1}(undef,0)
   social_nonworkday_dd::Array{Distribution,1} = Array{Distribution,1}(undef,0)
   n_groups_per_day_distribution::Array{Distribution,1} = Array{Distribution,1}(undef,0)
end

"""
   PreinterventionComponentStore()

Structure containing storage objects used to store the pre-intervention state of the network.
This is only used if the social contact structure is changed via intervention.

Location: parametertypes.jl
"""
@with_kw mutable struct PreinterventionComponentStore

   contacts_changed_by_intervention_time::Int64 = -1   # Time that intervention took place

   workday_social_contacts_by_day::Array{Array{Int64,1},2} = Array{Array{Int64,1},2}(undef, 0, 0)
   nonworkday_social_contacts_by_day::Array{Array{Int64,1},2} = Array{Array{Int64,1},2}(undef, 0, 0)

   network_generation_method_dynamic_social::String = ""
   social_workday_dd::Distribution = Distributions.Normal(0,0)
   social_nonworkday_dd::Distribution = Distributions.Normal(0,0)
   group_limit::Int64 = -1
   dynamic_time_frame::Int64 = -1
   n_groups_per_day_distribution::Distribution = Distributions.Normal(0,0)
   n_households_per_group::Int64 = -1

end
