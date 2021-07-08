"""
File containing functions used to generate workplace sizes
"""

"""
    workplace_size_sampled_from_empirical_pmf(args)

Generate single workplace size according to a discrete CDF that is loaded in.

Inputs: `rng` - random number generator,
        `worker_grp_idx` - current sector number,
        `workplace_generation_params` - WorkplaceGenerationParameters structure \n
Outputs: `workplace_size` - size of one workplace \n
Location: workplace_size_generation_fns.jl
"""
function workplace_size_sampled_from_empirical_pmf(rng::MersenneTwister,
                                                        worker_grp_idx::Int64,
                                                        workplace_generation_params::WorkplaceGenerationParameters)

    # Get data from spreadsheet/load from parameter type
    work_size_cdf_array = readdlm("../../Data/worker_model/work_sector_proportions_by_size.csv")
        # Loads a work_sector x n_bins array. CDF per row.

    # Sample based on worker_grp_idx
    # Bins: 0-4, 5-9, 10-19, 20-49,	50-99, 100-249, 250+
    selected_cdf = work_size_cdf_array[:,worker_grp_idx]
    chosen_bin = draw_sample_from_pmf(selected_cdf,rng)

    # If not in last bin, sample uniformly
    if chosen_bin == 1
        workplace_size = rand(rng,1:4)
    elseif chosen_bin == 2
        workplace_size = rand(rng,5:9)
    elseif chosen_bin == 3
        workplace_size = rand(rng,10:19)
    elseif chosen_bin == 4
        workplace_size = rand(rng,20:49)
    elseif chosen_bin == 5
        workplace_size = rand(rng,50:99)
    elseif chosen_bin == 6
        workplace_size = rand(rng,100:249)
    else
        # From last bin. Draw from a given distribution
        workplace_size = round(Int64,250 + rand(rng,Gamma(1, 100)))
    end

    return workplace_size::Int64
end

"""
    workplace_size_sampled_from_empirical_pmf(args)

Generate single workplace size according to a normal distribution with specified mean and SD.

Inputs: `rng` - random number generator,
        `worker_grp_idx` - current sector number,
        `workplace_generation_params` - WorkplaceGenerationParameters structure \n
Outputs: `workplace_size` - size of one workplace \n
Location: workplace_size_generation_fns.jl
"""
function workplace_size_gen_using_normal_dist(rng::MersenneTwister,
                                                    worker_grp_idx::Int64,
                                                    workplace_generation_params::WorkplaceGenerationParameters)

    @unpack workplace_size_mean, workplace_size_sd = workplace_generation_params
    workplace_size = ceil(Int64,abs(randn(rng)*workplace_size_sd[worker_grp_idx]+workplace_size_mean[worker_grp_idx]))

    return workplace_size::Int64
end
