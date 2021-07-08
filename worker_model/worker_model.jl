"""
Purpose:
Network model to explore impact of working patterns and other work-related interventions
on transmission in a worker population
"""

"""
Set paths & load environment
"""

#Set paths
cd(dirname(@__FILE__))

#Load environment
using Pkg
Pkg.activate("../../")

"""
Load packages
"""
#Required packages
using MAT, Distributions, CSV
using LinearAlgebra, Random, DelimitedFiles, Parameters

# using Profile  # For use with memory allocation checks

"""
Load supporting files
"""
# Files containing other functions needed to run the model
include("include_files_network_model/parametertypes.jl")
include("include_files_network_model/contact_tracing_fns.jl")
include("include_files_network_model/network_generation_fns.jl")
include("include_files_network_model/workplace_size_generation_fns.jl")
include("include_files_network_model/intervention_fns.jl")
include("include_files_network_model/seed_initial_states_fns.jl")
include("include_files_network_model/main_function.jl")
include("include_files_network_model/infection_process_fns.jl")
include("include_files_network_model/configuration_fns.jl")
include("include_files_network_model/initialisation_fns.jl")
include("include_files_network_model/workplace_closure_fns.jl")
include("include_files_network_model/assign_transmission_risk_fns.jl")

"""
Set variables from ARGS
"""
args = ARGS

# If running locally from REPL, will not have any input ARGS
# Can set values for input parameters to main run function here
# args/ARGS list
# args[1] JOB_ID
# args[2] RNGseed: To be used to initialise the random number generator
# args[3] COUNTFINAL: Number of replicates requested
# args[4] ENDTIME: Timesteps for each individual simulation
# args[5] n_nodes: Overall size of network
# args[6] n_sectors: Defining total number of work sectors
# args[7] runsets: Scenarios to run: [configuration, intervention] (see load_configurations() and load_interventions() for options)
if length(ARGS)==0
    args = ["1", "100", "100", "365", "1000", "41", "[\"default\" \"none\"]"]
end

# To run from command line, example:
# julia worker_pattern_model.jl 1 100 10 365 10000 41 set_ten_initial_infected '["run_one_run"]'

# Set identifier for job
const JOB_ID = parse(Int64, args[1])

# Set RNG seed
const RNGSEED_BASE = parse(Int64, args[2])
const RNGSEED = JOB_ID*RNGSEED_BASE

# Set simulation run params
const COUNTFINAL = parse(Int64, args[3])  # Number of simulations to be performed per scenario
const ENDTIME = parse(Int64, args[4]) # Timesteps for each individual simulation
const CMAX = parse(Int64, args[5]) # Number of individuals in the network
const WORKERTYPES = parse(Int64, args[6]) # Defining total number of work sectors

# Specify scenario to be run
const RUNSETS = eval(Meta.parse(args[7]))

for run_it = 1:length(RUNSETS[:,1])

    # Get network configuration / intervention pair
    configuration_name = RUNSETS[run_it, 1]
    intervention_name = RUNSETS[run_it, 2]

    # Load configuration variables
    configuration_variables = load_configurations(configuration_name)

    # Load intervention sets (one set per simulation - a set can contain multiple interventions)
    all_intervention_sets = load_interventions(intervention_name)

    n_configs = configuration_variables.n_configs
    n_intervention_sets = max(length(all_intervention_sets), 1)

    # Initialise output arrays. Store counts for each network configuration
    # put COUNTFINAL in the 2nd place so that we can concatenate variables later
    numlat_save = zeros(ENDTIME+1,COUNTFINAL,n_intervention_sets,n_configs)
    numinf_save = zeros(ENDTIME+1,COUNTFINAL,n_intervention_sets,n_configs)
    numrep_save = zeros(ENDTIME+1,COUNTFINAL,n_intervention_sets,n_configs)
    prevlat_save = zeros(ENDTIME+1,COUNTFINAL,n_intervention_sets,n_configs)
    prevsymp_save = zeros(ENDTIME+1,COUNTFINAL,n_intervention_sets,n_configs)
    prevasymp_save = zeros(ENDTIME+1,COUNTFINAL,n_intervention_sets,n_configs)
    prevpresymp_save = zeros(ENDTIME+1,COUNTFINAL,n_intervention_sets,n_configs)
    prevrec_save = zeros(ENDTIME+1,COUNTFINAL,n_intervention_sets,n_configs)
    newinf_save = zeros(ENDTIME+1,COUNTFINAL,n_intervention_sets,n_configs)
    workersinf_save = zeros(ENDTIME+1,COUNTFINAL,n_intervention_sets,n_configs)
    workersasymp_save = zeros(ENDTIME+1,COUNTFINAL,n_intervention_sets,n_configs)
    newasymp_save = zeros(ENDTIME+1,COUNTFINAL,n_intervention_sets,n_configs)
    num_isolating_save = zeros(ENDTIME+1,COUNTFINAL,n_intervention_sets,n_configs)      # Number isolating on given day
    num_symp_isolating_save = zeros(ENDTIME+1,COUNTFINAL,n_intervention_sets,n_configs) # Number isolating due to having symptoms
    num_household_isolating_save = zeros(ENDTIME+1,COUNTFINAL,n_intervention_sets,n_configs) # Number isolating due to having symptoms
    num_CT_save = zeros(ENDTIME+1,COUNTFINAL,n_intervention_sets,n_configs)             # Number contact traced on given day
    num_infected_save = Array{Array{Int64,3}}(undef,n_configs) #zeros(Int64,cmax,COUNTFINAL,n_configs)
    dynamic_infection_count_save = zeros(Int64,ENDTIME+1,COUNTFINAL,n_intervention_sets,n_configs)
    var_num_infected_save = zeros(Float64,1,COUNTFINAL,n_intervention_sets,n_configs) # so that the COUNTFINAL is in the 2nd place
    mean_init_generation_time_save = zeros(Float64,1,COUNTFINAL,n_intervention_sets,n_configs) # so that the COUNTFINAL is in the 2nd place
    num_isolating_CTcause_save = zeros(ENDTIME+1,COUNTFINAL,n_intervention_sets,n_configs)
    Rt_save = zeros(ENDTIME+1,COUNTFINAL,n_intervention_sets,n_configs)
    transmission_setting_save = zeros(Int64, ENDTIME+1, COUNTFINAL, n_intervention_sets, 5, n_configs)
    tests_performed_save = zeros(ENDTIME+1,COUNTFINAL,n_intervention_sets,n_configs)
    test_outcomes_save = zeros(ENDTIME+1,COUNTFINAL,n_intervention_sets,4,n_configs)
    num_init_infected_save = Array{Array{Array{Int64,2}}}(undef,n_configs)
    social_meeting_groups_save = zeros(Float64, 8, n_intervention_sets+1, n_configs)

    for config_idx = 1:n_configs

        # Load default parameters
        configuration_options, network_params,
        contact_tracing_params, infection_params,
        workplace_generation_params, workplace_closure_params = load_defaults()

        apply_config_changes!(config_idx, configuration_variables, configuration_options,
                            network_params, contact_tracing_params, infection_params,
                            workplace_generation_params, workplace_closure_params)

        # Set up array to store number infected by each node
        num_infected_save[config_idx] = zeros(Int64, configuration_options.cmax,
                                            COUNTFINAL,
                                            n_intervention_sets)

        # Call function (located in include_files_network_model\main_function.jl)
        @time  output = worker_pattern_network_run(config_idx, configuration_variables,
                                                configuration_options, network_params,
                                                contact_tracing_params, infection_params,
                                                workplace_generation_params, workplace_closure_params,
                                                all_intervention_sets)
        @unpack numlat, numinf, numrep,
                prevlat, prevsymp, prevasymp, prevpresymp, prevrec,
                newinf, newasymp, workersinf, workersasymp, infected_by,
                num_isolating, num_household_isolating, num_symp_isolating, num_isolating_CTcause,
                num_CT, num_infected, dynamic_infection_count,
                var_num_infected, num_init_infected, mean_init_generation_time,
                Rt, transmission_setting, tests_performed, test_outcomes, social_meeting_groups = output

        # For this iteration's network config, write results to output storage arrays
        numlat_save[:,:,:,config_idx] = numlat
        numinf_save[:,:,:,config_idx] = numinf
        numrep_save[:,:,:,config_idx] = numrep
        prevlat_save[:,:,:,config_idx] = prevlat
        prevsymp_save[:,:,:,config_idx] = prevsymp
        prevasymp_save[:,:,:,config_idx] = prevasymp
        prevpresymp_save[:,:,:,config_idx] = prevpresymp
        prevrec_save[:,:,:,config_idx] = prevrec
        newinf_save[:,:,:,config_idx] = newinf
        workersinf_save[:,:,:,config_idx] = workersinf
        workersasymp_save[:,:,:,config_idx] = workersasymp
        newasymp_save[:,:,:,config_idx] = newasymp
        num_isolating_save[:,:,:,config_idx] = num_isolating
        num_symp_isolating_save[:,:,:,config_idx] = num_symp_isolating
        num_household_isolating_save[:,:,:,config_idx] = num_household_isolating
        num_CT_save[:,:,:,config_idx] = num_CT
        num_infected_save[config_idx][:,:,:] = num_infected
        dynamic_infection_count_save[:,:,:,config_idx] = dynamic_infection_count
        var_num_infected_save[1,:,:,config_idx] = var_num_infected
        num_init_infected_save[config_idx] = num_init_infected
        mean_init_generation_time_save[1,:,:,config_idx] = mean_init_generation_time
        Rt_save[:,:,:,config_idx] = Rt
        transmission_setting_save[:,:,:,:,config_idx] = transmission_setting
        tests_performed_save[:,:,:,config_idx] = tests_performed
        test_outcomes_save[:,:,:,:,config_idx] = test_outcomes

        num_isolating_CTcause_save[:,:,:,config_idx] = num_isolating_CTcause

        social_meeting_groups_save[:,:,config_idx] = social_meeting_groups
    end

    output_file = "/Users/benjaminatkins/Documents/worker_model_results/worker_model_output_$(configuration_name)_$(intervention_name)_cmax$(CMAX)_sectors$(WORKERTYPES)_#$(JOB_ID).mat"
    file = matopen(output_file, "w")

    write(file,"numlat",numlat_save)
    write(file,"numinf",numinf_save)
    write(file,"numrep",numrep_save)
    write(file,"prevlat",prevlat_save)
    write(file,"prevsymp",prevsymp_save)
    write(file,"prevasymp",prevasymp_save)
    write(file,"prevpresymp",prevpresymp_save)
    write(file,"prevrec",prevrec_save)
    write(file,"newinf",newinf_save)
    write(file,"workersinf",workersinf_save)
    write(file,"workersasymp",workersasymp_save)
    write(file,"newasymp",newasymp_save)
    write(file, "num_isolating", num_isolating_save)
    write(file, "num_symp_isolating", num_symp_isolating_save)
    write(file, "num_household_isolating", num_household_isolating_save)
    write(file, "num_infected", num_infected_save)
    write(file, "dynamic_infection_count", dynamic_infection_count_save)
    write(file, "var_num_infected_save", var_num_infected_save)
    write(file, "num_init_infected_save", num_init_infected_save)
    write(file, "mean_init_generation_time_save", mean_init_generation_time_save)
    write(file, "Rt_save", Rt_save)
    write(file, "transmission_setting", transmission_setting_save)
    write(file, "tests_performed", tests_performed_save)
    write(file, "test_outcomes", test_outcomes_save)
    write(file, "num_CT", num_CT_save)
    write(file, "num_isolating_CTcause", num_isolating_CTcause_save)
    write(file, "social_meeting_groups", social_meeting_groups_save)
    close(file)
end
