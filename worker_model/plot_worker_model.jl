"""
Purpose:
Plot results from network model - requires network model to have been run and results saved
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
using MAT, Distributions, CSV, PyPlot
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
    args = ["1", "100", "100", "365", "1000", "41", "[\"default\" \"test_intervention\"]"]
end

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

pygui(true)

for run_it = 1:length(RUNSETS[:,1])

    # Get network configuration / intervention pair
    configuration_name = RUNSETS[run_it, 1]
    intervention_name = RUNSETS[run_it, 2]

    # Load configuration variables
    configuration_variables = load_configurations(configuration_name)
    n_configs = configuration_variables.n_configs

    # Load intervention sets (one set per simulation - a set can contain multiple interventions)
    all_intervention_sets = load_interventions(intervention_name)
    n_intervention_sets = max(length(all_intervention_sets), 1)

    file = matopen("/Users/benjaminatkins/Documents/worker_model_results/worker_model_output_$(configuration_name)_$(intervention_name)_cmax$(CMAX)_sectors$(WORKERTYPES)_#$(JOB_ID).mat", "r")
    newinf = read(file, "newinf")
    close(file)

    newinf_weekly = Array{Float64,4}(undef, floor(Int64,ENDTIME/7), COUNTFINAL, n_intervention_sets, n_configs)
    for i = 1:floor(Int64,ENDTIME/7)
        newinf_weekly[i,:,:,:] = sum(newinf[((i-1)*7+1):(i*7),:,:,:],dims=(1))[1,:,:,:]
    end

    fig, ax = PyPlot.subplots(n_intervention_sets, n_configs, figsize = (10*n_configs,2*n_intervention_sets), sharex=true, sharey=true)

    for row_itr = 1:n_intervention_sets
        for col_itr = 1:n_configs

            if (n_intervention_sets > 1) && (n_configs > 1)
                ax[row_itr, col_itr].plot(newinf_weekly[:,:,row_itr,col_itr])
            elseif n_configs == 1
                ax[row_itr].plot(newinf_weekly[:,:,row_itr,col_itr])
            elseif n_intervention_sets == 1
                ax[col_itr].plot(newinf_weekly[:,:,row_itr,col_itr])
            end
        end
    end

    fig.savefig("/Users/benjaminatkins/Documents/worker_model_results/weekly_newinf_$(configuration_name)_$(intervention_name)_cmax$(CMAX)_sectors$(WORKERTYPES)_#$(JOB_ID).png")
end
