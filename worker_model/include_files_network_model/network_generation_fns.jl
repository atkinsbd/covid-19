"""
File containing functions to generate each layer of the network
"""

"""
    generate_workplaces_and_allocate_workers!(args)

Randomly generate workplace sizes and assign each node to a workplace.

The following information is stored in NetworkParameters structure during this process:
- `worker_nodes`: array of WorkerParameters structures containing info on each node
- `workplace_sizes`: array{array} of workplace sizes, by sector and workplace number
- `workplace_info`: array{array} of WorkplaceParameters structures, by sector and workplace number
- `nodes_by_workplace`: array{array{array}} of node IDs, sorted into workplaces by sector and workplace number
- `nodes_by_sector`: array{array} of node ids, collected across workplaces into sectors

Inputs: `... parameter structures ...`,
        `rng` - random number generator \n
Outputs: None \n
Location: network_generation_fns.jl
"""
function generate_workplaces_and_allocate_workers!(configuration_options::ConfigurationOptions,
                                            workplace_generation_params::WorkplaceGenerationParameters,
                                            network_params::NetworkParameters,
                                            rng::MersenneTwister)

    # Get workplace params based on number of work sectors in use
    @unpack workpercent, workforce_proportion, workplace_size_mean, workplace_size_sd, workplace_size_gen_fn = workplace_generation_params
    @unpack cmax, workertypes, CS_active_flag = configuration_options

    # total number of workers in each worker type
    worker_numbers = [round(Int64,cmax*workforce_proportion[i]) for i=1:length(workforce_proportion)]
    # correct for rounding
    if sum(worker_numbers) != cmax
        diff = sum(worker_numbers) - cmax
        worker_numbers[end] -= diff
    end

    # Initialise places to store workplace sizes per workplace per worker group
    # & workplace information (covid secure, workplace open)
    workplace_sizes = Array{Array{Int64,1},1}(undef, workertypes)
    workplace_info = Array{Array{WorkplaceParameters,1},1}(undef, workertypes)

    # Initialise ID vars
    workplace_ids = Array{Array{Int64,1},1}(undef, workertypes)  # representing workplace id per person per group
    workertype_ids = Int64[]  # representing workertype id per person

    # Initialise array to store ids of workers within each workplace
    nodes_by_workplace = Array{Array{Array{Int64,1},1},1}(undef,workertypes)

    if isassigned(network_params.sector_open) == false
        # Initialise sector open/closure status
        sector_open = Array{Bool,1}(undef,workertypes)
        for sector_itr = 1:workertypes
            sector_open[sector_itr] = true
        end
        network_params.sector_open = sector_open
    else
        sector_open = network_params.sector_open
    end

    ## Now generate workplaces within each worker type
    for worker_grp_idx = 1:workertypes
        # initialise arrays for this worker type
        workplace_info[worker_grp_idx] =  WorkplaceParameters[]
        workplace_sizes[worker_grp_idx] = Int64[]
        workplace_ids[worker_grp_idx] = Int64[]

        # Initialise array to store ids of workers within each workplace
        nodes_by_workplace[worker_grp_idx] = Array[]

        # Set up tracking counters
        workplace_count = 0   # track number of workplaces in this worker type
        worker_count = 0        # track number of workers assigned to workplaces

        # Add workplaces until worker_count exceeds target number of workers in that sector
        while worker_count < worker_numbers[worker_grp_idx]

            # Update workplace number counter
            workplace_count += 1

            # generate workplace sizes and add to array
            # For workplace_size_gen_fn options, see "include_files_network_model/workplace_size_generation_fns.jl"
            push!(workplace_sizes[worker_grp_idx], workplace_size_gen_fn(rng,worker_grp_idx,workplace_generation_params))
            append!(workplace_ids[worker_grp_idx], repeat([workplace_count], workplace_sizes[worker_grp_idx][end]))

            # Add empty workplace to nodes_by_workplace array
            push!(nodes_by_workplace[worker_grp_idx], Int64[])

            # Update tracking counter values
            worker_count += workplace_sizes[worker_grp_idx][end]
        end

        # remove extra workers
        workplace_sizes[worker_grp_idx][end] = worker_numbers[worker_grp_idx] - sum(workplace_sizes[worker_grp_idx][1:end-1])
        workplace_ids[worker_grp_idx] = workplace_ids[worker_grp_idx][1:sum(workplace_sizes[worker_grp_idx])]

        # randomly shuffle workers between workplaces
        shuffle!(rng,workplace_ids[worker_grp_idx])

        append!(workertype_ids, repeat([worker_grp_idx],sum(workplace_sizes[worker_grp_idx])))

        # Create workplace parameter type for this worker type
        workplace_info[worker_grp_idx] =  Array{WorkplaceParameters,1}(undef,workplace_count)
        for workplace_itr = 1:workplace_count
            workplace_info[worker_grp_idx][workplace_itr] = WorkplaceParameters(covid_secure = CS_active_flag,
                                                                                workplace_open = sector_open[worker_grp_idx])
        end

    end
    # randomly shuffle workertypes
    shuffle!(rng,workertype_ids)

    worker_nodes = Array{WorkerParameters,1}(undef,cmax)    ## (atwork, workertype, workplaceid)

    for ii = 1:cmax
        worker_nodes[ii] = WorkerParameters(sector_ID = workertype_ids[ii],
                                            workplace_ID = workplace_ids[workertype_ids[ii]][1])

        # All workers added to workplace, regardless of returned_to_work
        push!(nodes_by_workplace[workertype_ids[ii]][workplace_ids[workertype_ids[ii]][1]], ii)
        deleteat!(workplace_ids[workertype_ids[ii]],1) # remove worker from top of list
    end

    network_params.nodes_by_sector = [collect(Iterators.flatten(nodes_by_workplace[sector_itr][:])) for sector_itr = 1:workertypes]

    network_params.worker_nodes::Array{WorkerParameters,1} = worker_nodes
    network_params.workplace_sizes::Array{Array{Int64,1},1} = workplace_sizes
    network_params.workplace_info::Array{Array{WorkplaceParameters,1},1} = workplace_info
    network_params.nodes_by_workplace::Array{Array{Array{Int64,1},1},1} = nodes_by_workplace

    return nothing
end

"""
    generate_contacts(args)

Generate static contact networks for households, workplaces and social groups, according to configuration options and network parameters.

Returns the number of households and a ContactStructure structure containing the following information:
- `work_contacts_same_workplace`: array{array} of node IDs, for static contacts made by each node in same workplace
- `work_contacts_same_workplace_per_node`: array containing number of static contacts made in same workplace per node
- `work_contacts_other_workplace`: array{array} of node IDs, for static contacts made by each node in other workplace
- `work_contacts_other_workplace_per_node`: array containing number of static contacts made in other workplace per node
- `social_contacts`: array{array} of node IDs, for static social groups of each node
- `social_contacts_per_node`: array containing size of social group per node (not including self)
- `household_contacts`: array{array} of node IDs, for static household contacts made by each node
- `household_contacts_per_node`: array containing size of households per node (not including self)
- `nodes_by_household`: array{array} of node IDs by household

If COVID-secure workplace measures are in place, the ContactStructure will also contain:
- `work_contacts_same_workplace_CS`: array{array} of node IDs, for static contacts made by each node in COVID-secure workplace
- `work_contacts_same_workplace_per_node_CS`: array containing number of static contacts made in COVID-secure workplace per node

Inputs: `... parameter structures ...`,
        `rng` - random number generator \n
Outputs: `contacts` - ContactStructure as described above,
         `n_households` - number of households \n
Location: network_generation_fns.jl
"""
function generate_contacts(configuration_options::ConfigurationOptions,
                            network_params::NetworkParameters,
                            rng::MersenneTwister)

    @unpack worker_nodes, workplace_sizes,prob_workertype_contact,prob_anyworker_contact,
        prob_social_contact,dd_within_workplace, household_size_distribution,
        workplace_info, workplace_degree_distribution,
        between_workplace_contact_probs, social_group_size_distribution,
        friend_of_friend_prob, max_contacts_social, network_generation_method,
        CS_team_size, nodes_by_workplace = network_params

    @unpack cmax, endtime, CS_active_flag, workertypes = configuration_options

    # Initialise vector of vectors storing IDs of contacts for each node
    # at workplace (non-covid secure setting)
    work_contacts_same_workplace = Array{Array{Int64,1},1}(undef,cmax)
    work_contacts_other_workplace = Array{Array{Int64,1},1}(undef,cmax)

    # Initialise vector of vectors storing IDs of contacts for each node in social
    # and household settings
    social_contacts = Array{Array{Int64,1},1}(undef,cmax)
    household_contacts = Array{Array{Int64,1},1}(undef,cmax)

    # If CS settings are active,
    # initialise vector of vectors storing IDs of contacts for each node
    # at workplace (covid-secure setting) & vector giving total contacts made by node
    if CS_active_flag == true
        work_contacts_same_workplace_CS = Array{Array{Int64,1},1}(undef,cmax)
        work_contacts_same_workplace_per_node_CS = zeros(Int64,cmax)
    end

    for ii = 1:cmax
        work_contacts_same_workplace[ii] = Int64[]
        work_contacts_other_workplace[ii] = Int64[]
        social_contacts[ii] = Int64[]
        household_contacts[ii] = Int64[]

        if CS_active_flag == true
            work_contacts_same_workplace_CS[ii] = Int64[]
        end
    end

    # Initialise vectors giving total contacts made by node
    work_contacts_same_workplace_per_node = zeros(Int64,cmax)
    work_contacts_other_workplace_per_node = zeros(Int64,cmax)
    social_contacts_per_node = zeros(Int64,cmax)
    household_contacts_per_node = zeros(Int64,cmax)

    # Initialise array to store social groups
    nodes_by_social_group = Array[]

    if network_generation_method == "configuration"

        if CS_active_flag == false

            # Construct links for workplace

            if length(workplace_degree_distribution)==1
                workplace_degree_distribution = repeat(workplace_degree_distribution, workertypes)
            end
            if length(between_workplace_contact_probs)==1
                between_workplace_contact_probs = repeat(between_workplace_contact_probs, workertypes)
            end

            worker_count = 0

            # Cycle through sectors and workplaces
            for worker_grp_idx = 1:workertypes

                # nodes_outside_cluster = collect(Iterators.flatten(nodes_by_workplace[worker_grp_idx]))

                # find nodes in the given sector (these are nodes outside the cluster of a given workplace)
                length_array = 0
                for ii=1:length(nodes_by_workplace[worker_grp_idx])
                    length_array += length(nodes_by_workplace[worker_grp_idx][ii])
                end

                nodes_outside_cluster = zeros(Int64,length_array)
                init_pos = 1
                for ii=1:length(nodes_by_workplace[worker_grp_idx])
                    end_pos = init_pos+length(nodes_by_workplace[worker_grp_idx][ii])-1
                    nodes_outside_cluster[init_pos:end_pos] = nodes_by_workplace[worker_grp_idx][ii]
                    init_pos = end_pos+1
                end

                # set the degree distribution and external contact probablity for that sector
                degree_distribution = workplace_degree_distribution[worker_grp_idx]
                external_contact_prob = between_workplace_contact_probs[worker_grp_idx]

                # iterate through the workplaces in that sector
                for workplace_idx = 1:length(workplace_sizes[worker_grp_idx])

                    # nodes in that workplace are in a cluster
                    nodes_within_cluster = nodes_by_workplace[worker_grp_idx][workplace_idx]

                    if length(nodes_within_cluster) > 0
                        worker_count += length(nodes_within_cluster)
                        configuration_model!(external_contact_prob,
                                            degree_distribution,
                                            nodes_within_cluster,
                                            nodes_outside_cluster,
                                            work_contacts_same_workplace,
                                            work_contacts_other_workplace,
                                            work_contacts_same_workplace_per_node,
                                            work_contacts_other_workplace_per_node,
                                            network_params,
                                            rng)
                    end
                end
            end
            total_ext = sum(work_contacts_other_workplace_per_node)/2
            total_int = sum(work_contacts_same_workplace_per_node)/2
            println("External workplace contacts: $total_ext")
            println("Internal workplace contacts: $total_int")

        ## CS workplace generation
        else
            CS_workplace_generation!(nodes_by_workplace,
                                    work_contacts_same_workplace_CS,
                                    work_contacts_same_workplace_per_node_CS,
                                    network_params,
                                    rng)

            total_int = sum(work_contacts_same_workplace_per_node_CS)/2
            println("Internal CS workplace contacts: $total_int")
        end

        degree_distribution = social_group_size_distribution

        configuration_model!(cmax,
                                friend_of_friend_prob,
                                degree_distribution,
                                social_contacts,
                                social_contacts_per_node,
                                max_contacts_social,
                                network_params,
                                rng)


    elseif network_generation_method == "ER"

        ER_model!(worker_nodes, cmax, workplace_sizes,
                            dd_within_workplace, prob_workertype_contact, prob_anyworker_contact,
                            work_contacts_same_workplace::Array{Array{Int64,1},1},
                            work_contacts_other_workplace::Array{Array{Int64,1},1},
                            work_contacts_same_workplace_per_node::Array{Int64,1},
                            work_contacts_other_workplace_per_node::Array{Int64,1},
                            rng::MersenneTwister)

        # Construct network links for social
        for ii = 1:(cmax-1)
            for jj = (ii+1):cmax
                ### Social contacts with anyone ###
                if rand(rng) < prob_social_contact
                    # Assign IDs of contact to tuple vectors
                    push!(social_contacts[ii],jj)
                    push!(social_contacts[jj],ii)
                    # Increment number of contacts for nodes ii & jj
                    social_contacts_per_node[ii] += 1
                    social_contacts_per_node[jj] += 1
                end
            end
        end
    end

    # Create household contacts & assign a household ID
    nodes_by_household::Array{Array{Int64,1},1} = Array[]
    csum_household_size_distribution = cumsum(household_size_distribution) # Get cumulative sum, used in each repeat of while loop
    worker_id = 1     # Initialise variables to be incremented in while loop
    household_id = 1
    while worker_id <= cmax
        household_size = findfirst(x->x>rand(rng), csum_household_size_distribution)  # generate random household size

        # check we don't exceed remaining workers left
        household_size = min(household_size, cmax-worker_id+1)

        # create fully-connected household
        for ii = 0:(household_size-1)
            for jj = 0:(household_size-1)
                # Get household contacts
                if ii != jj
                    push!(household_contacts[worker_id+ii], worker_id+jj)
                    household_contacts_per_node[worker_id+ii] += 1
                end
            end

            # Assign household ID to node ii
            worker_nodes[worker_id+ii].household_ID = household_id
        end

        push!(nodes_by_household, shuffle(rng, [worker_id:(worker_id+household_size-1);]))

        worker_id += household_size
        household_id += 1
    end

    # Assign number of households in total to output variable
    n_households = household_id

    # Define what is returned from the function
    if CS_active_flag == true
        contacts = ContactStructure(
        work_contacts_same_workplace=work_contacts_same_workplace,
        work_contacts_same_workplace_per_node = work_contacts_same_workplace_per_node,
        work_contacts_other_workplace=work_contacts_other_workplace,
        work_contacts_other_workplace_per_node = work_contacts_other_workplace_per_node,
        social_contacts = social_contacts,
        social_contacts_per_node = social_contacts_per_node,
        household_contacts = household_contacts,
        household_contacts_per_node = household_contacts_per_node,
        nodes_by_household = nodes_by_household,
        work_contacts_same_workplace_CS = work_contacts_same_workplace_CS,
        work_contacts_same_workplace_per_node_CS = work_contacts_same_workplace_per_node_CS)
    else
        contacts = ContactStructure(
        work_contacts_same_workplace=work_contacts_same_workplace,
        work_contacts_same_workplace_per_node = work_contacts_same_workplace_per_node,
        work_contacts_other_workplace=work_contacts_other_workplace,
        work_contacts_other_workplace_per_node = work_contacts_other_workplace_per_node,
        social_contacts = social_contacts,
        social_contacts_per_node = social_contacts_per_node,
        household_contacts = household_contacts,
        household_contacts_per_node = household_contacts_per_node,
        nodes_by_household = nodes_by_household)
    end

    return contacts::ContactStructure,
        n_households::Int64
end

"""
    CS_workplace_generation!(args)

Generate workplace contact layer under COVID-secure restrictions.

Workplaces are split into fully connected, mutually exclusive groups of specified size.

Inputs: `nodes_by_workplace` - array{array{array}} of node IDs grouped by sector and workplace,
        `work_contacts_same_workplace_CS` - array{array} of contact IDs per node,
        `work_contacts_same_workplace_per_node_CS` - array of number of contacts per node,
        `network_params` - NetworkParameters structure,
        `rng` - random number generator \n
Outputs: None \n
Location: network_generation_fns.jl
"""
function CS_workplace_generation!(nodes_by_workplace::Array{Array{Array{Int64,1},1},1},
                                    work_contacts_same_workplace_CS::Array{Array{Int64,1},1},
                                    work_contacts_same_workplace_per_node_CS::Array{Int64,1},
                                    network_params::NetworkParameters,
                                    rng::MersenneTwister)

    @unpack workplace_sizes, CS_team_size = network_params

    workertypes = length(nodes_by_workplace)

    # Cycle through sectors and workplaces
    for worker_grp_idx = 1:workertypes

        n_workplaces = length(nodes_by_workplace[worker_grp_idx])

        # iterate through the workplaces in that sector
        for workplace_idx = 1:n_workplaces

            # Find nodes in this workplace
            nodes_within_workplace = nodes_by_workplace[worker_grp_idx][workplace_idx]
            # Number of nodes in workplace
            n_nodes_within_workplace = workplace_sizes[worker_grp_idx][workplace_idx]
            # Number of isolated, fully connected groups of size CS_team_size
            n_CS_groups = floor(Int64, n_nodes_within_workplace/CS_team_size)

            # Create fully connected groups
            for CS_group_idx = 1:n_CS_groups
                current_group = nodes_within_workplace[((CS_group_idx-1)*CS_team_size+1):(CS_group_idx*CS_team_size)]
                for worker_idx = 1:CS_team_size
                    node_id = current_group[worker_idx]
                    work_contacts_same_workplace_CS[node_id] = current_group[1:end .!= worker_idx]
                    work_contacts_same_workplace_per_node_CS[node_id] = CS_team_size - 1
                end
            end

            # Group together any remaining workers in workplace
            workers_remaining = Int64(n_nodes_within_workplace%CS_team_size)

            if workers_remaining > 1
                current_group = nodes_within_workplace[(n_CS_groups*CS_team_size+1):end]
                for worker_idx = 1:workers_remaining
                    node_id = current_group[worker_idx]
                    work_contacts_same_workplace_CS[node_id] = current_group[1:end .!= worker_idx]
                    work_contacts_same_workplace_per_node_CS[node_id] = workers_remaining - 1
                end
            end
        end
    end

    network_params.CS_contacts_previously_generated = true

    return nothing
end


"""
    configuration_model!(args)

Generate STATIC WORKPLACE contact layer using configuration model style algorithm.

Contacts are generated according to a specified degree distribution, occurring within
the same workplace or with another workplace (in the same sector) according to a specified probability.
Repeat contacts and self-contacts are not allowed.

Inputs: `external_contact_prob` - probability that contact is made outside of workplace,
        `degree_distribution` - desired degree distribution of contacts,
        `nodes_within_cluster` - array of node IDs belonging to cluster (e.g. workplace),
        `nodes_outside_cluster` - array of node IDs not belonging to cluster, but able to make contacts, (e.g. other workplaces in sector),
        `work_contacts_same_workplace` - array{array} of contact IDs within same workplace per node,
        `work_contacts_other_workplace` - array{array} of contact IDs with other workplace per node,
        `work_contacts_same_workplace_per_node` - array of number of contacts within same workplace per node,
        `work_contacts_other_workplace_per_node` - array of number of contacts with other workplace per node,
        `network_params` - NetworkParameters structure,
        `rng` - random number generator \n
Outputs: None \n
Location: network_generation_fns.jl
"""
function configuration_model!(external_contact_prob::Float64,
                                degree_distribution::Distribution,
                                nodes_within_cluster::Array{Int64,1},
                                nodes_outside_cluster::Array{Int64,1},
                                work_contacts_same_workplace::Array{Array{Int64,1},1},
                                work_contacts_other_workplace::Array{Array{Int64,1},1},
                                work_contacts_same_workplace_per_node::Array{Int64,1},
                                work_contacts_other_workplace_per_node::Array{Int64,1},
                                network_params::NetworkParameters,
                                rng::MersenneTwister)


    n_nodes = length(nodes_within_cluster)
    n_nodes_external = length(nodes_outside_cluster)

    edges_per_node = Distributions.rand(rng, degree_distribution, n_nodes)

    # Round degree to nearest whole number
    # Decrease by proportion of contacts made external to workplace
    # Limit maximum number of contacts to workplace_size - 1
    for node_id = 1:n_nodes
        edges_per_node[node_id] = round(Int64, (edges_per_node[node_id])/(1+external_contact_prob))
        edges_per_node[node_id] = min((n_nodes - 1), edges_per_node[node_id])
    end

    half_edges = cumsum(edges_per_node)

    n_stubs = half_edges[end]

    edges_within_group = zeros(Int64, n_nodes, n_nodes)
    [edges_within_group[i,i] = 1 for i=1:n_nodes]

    n_nodes_external_remaining = n_nodes_external

    while half_edges[end] > 1

        stub_id1 = round(Int64, rand(rng)*(half_edges[end]-1) + 1)
        node_id1 = findfirst(x -> x >= stub_id1, half_edges)

        # Contact made external to cluster
        if (rand(rng) < external_contact_prob) & (n_nodes_external_remaining > n_nodes)

            node_id2 = round(Int64, rand(rng)*(n_nodes_external-1) + 1)

            # Don't allow links within cluster or repeated links
            while (nodes_outside_cluster[node_id2] in work_contacts_other_workplace[nodes_within_cluster[node_id1]]) || (node_id2 in nodes_within_cluster)
                node_id2 = round(Int64, rand(rng)*(n_nodes_external-1) + 1)
            end

            push!(work_contacts_other_workplace[nodes_within_cluster[node_id1]], nodes_outside_cluster[node_id2])
            push!(work_contacts_other_workplace[nodes_outside_cluster[node_id2]], nodes_within_cluster[node_id1])

            work_contacts_other_workplace_per_node[nodes_within_cluster[node_id1]] += 1
            work_contacts_other_workplace_per_node[nodes_outside_cluster[node_id2]] += 1

            n_nodes_external_remaining -= 1

            for ii=node_id1:n_nodes
                half_edges[ii] -= 1
            end

        # Contact made within cluster
        else
            # This ensures no self-links or repeated edges within cluster
            nodes_remaining = findall(edges_within_group[node_id1,:] .== 0)
            half_edges_remaining = cumsum(edges_per_node[nodes_remaining])
            # removals = findall(half_edges_remaining.==0)
            # deleteat!(nodes_remaining,removals)
            # deleteat!(half_edges_remaining,removals)
            nodes_remaining = nodes_remaining[half_edges_remaining.>0]
            half_edges_remaining = half_edges_remaining[half_edges_remaining.>0]

            # Contact not made (half edge lost)
            if length(nodes_remaining) < 1

                for ii=node_id1:n_nodes
                    half_edges[ii] -= 1
                end
            else
                stub_id2 = round(Int64, rand(rng)*(half_edges_remaining[end]-1) + 1)
                node_id2 = nodes_remaining[findfirst(x -> x >= stub_id2, half_edges_remaining)]

                push!(work_contacts_same_workplace[nodes_within_cluster[node_id1]], nodes_within_cluster[node_id2])
                push!(work_contacts_same_workplace[nodes_within_cluster[node_id2]], nodes_within_cluster[node_id1])

                work_contacts_same_workplace_per_node[nodes_within_cluster[node_id1]] += 1
                work_contacts_same_workplace_per_node[nodes_within_cluster[node_id2]] += 1

                edges_within_group[node_id1,node_id2] += 1

                for ii=node_id1:n_nodes
                    half_edges[ii] -= 1
                end
                for ii=node_id2:n_nodes
                    half_edges[ii] -= 1
                end
            end
        end
    end

end

"""
    configuration_model!(args)

Generate STATIC SOCIAL contact layer using "friends-of-friends" configuration model style algorithm.

Contacts are generated according to a specified degree distribution, occurring with a 'friend-of-friend'
or any other node according to a specified probability. Repeat contacts and self-contacts are not allowed.

Inputs: `n_nodes` - number of nodes in network,
        `friend_of_friend_prob` - probability that contact is made with friend-of-friend,
        `degree_distribution` - desired degree distribution of contacts,
        `social_contacts` - array{array} of social contact IDs per node,
        `social_contacts_per_node` - array of number of social contacts per node,
        `max_contacts_social` - maximum number of social contactgs allowed per node,
        `network_params` - NetworkParameters structure,
        `rng` - random number generator \n
Outputs: None \n
Location: network_generation_fns.jl
"""
function configuration_model!(n_nodes::Int64,
                                friend_of_friend_prob::Float64,
                                degree_distribution::Distribution,
                                social_contacts::Array{Array{Int64,1},1},
                                social_contacts_per_node::Array{Int64,1},
                                max_contacts_social::Int64,
                                network_params::NetworkParameters,
                                rng::MersenneTwister)

    edges_per_node = Distributions.rand(rng, degree_distribution, n_nodes)

    # Round degree to nearest whole number
    # Limit maximum number of contacts to max_contacts_social
    for node_id = 1:n_nodes
        edges_per_node[node_id] = round(Int64, edges_per_node[node_id])
        edges_per_node[node_id] = min(max_contacts_social, edges_per_node[node_id])
    end

    half_edges = cumsum(edges_per_node)

    if half_edges[end] % 2 != 0
        half_edges[end] += 1
    end

    n_stubs = half_edges[end]

    edges_remaining_per_node = copy(edges_per_node)
    fof = Int64[]
    while half_edges[end] > 1

        stub_id1 = round(Int64, rand(rng)*(half_edges[end]-1) + 1)
        node_id1 = findfirst(x -> x >= stub_id1, half_edges)

        # fof = collect(Iterators.flatten(social_contacts[social_contacts[node_id1]]))
        fof = find_fof(social_contacts,social_contacts[node_id1],edges_remaining_per_node,node_id1,fof)

        # Contact made with friend of friend
        if (length(fof) > 0) & (rand(rng) < friend_of_friend_prob)

            node_id2 = fof[ceil(Int64, rand(rng)*length(fof))]

        # Contact made with not friend of friend
        else

            stub_id2 = round(Int64, rand(rng)*(half_edges[end]-1) + 1)
            node_id2 = findfirst(x -> x >= stub_id2, half_edges)

        end

        push!(social_contacts[node_id1], node_id2)
        push!(social_contacts[node_id2], node_id1)

        social_contacts_per_node[node_id1] += 1
        social_contacts_per_node[node_id2] += 1

        edges_remaining_per_node[node_id1] -= 1
        edges_remaining_per_node[node_id2] -= 1


        for ii=node_id1:n_nodes
            half_edges[ii] -= 1
        end
        for ii=node_id2:n_nodes
            half_edges[ii] -= 1
        end

    end

end

"""
    ER_model!(args)

Generate STATIC WORKPLACE contact layer using ER style algorithm.

Contacts are generated according to a specified mean degree. Every unique contact pair is tested,
and accepted with a probability defined by the specified mean degree. This can be different for
contacts within the same workplace and different workplaces.

Inputs: `worker_nodes` - array of WorkerParameters structures containing worker info,
        `cmax` - number of nodes in network,
        `workplace_sizes` - array{array} of workplace sizes, grouped by sector,
        `dd_within_workplace` - desired mean degree within workplaces, per sector,
        `prob_workertype_contact` - probability of contacts with different workplaces, same sector (per sector),
        `prob_anyworker_contact` - probability of contacts with different workplaces, different sector (single value),
        `work_contacts_same_workplace` - array{array} of contact IDs within same workplace per node,
        `work_contacts_other_workplace` - array{array} of contact IDs with other workplace per node,
        `work_contacts_same_workplace_per_node` - array of number of contacts within same workplace per node,
        `work_contacts_other_workplace_per_node` - array of number of contacts with other workplace per node,
        `rng` - random number generator \n
Outputs: None \n
Location: network_generation_fns.jl
"""
function ER_model!(worker_nodes, cmax, workplace_sizes,
                    dd_within_workplace, prob_workertype_contact, prob_anyworker_contact,
                    work_contacts_same_workplace::Array{Array{Int64,1},1},
                    work_contacts_other_workplace::Array{Array{Int64,1},1},
                    work_contacts_same_workplace_per_node::Array{Int64,1},
                    work_contacts_other_workplace_per_node::Array{Int64,1},
                    rng::MersenneTwister)

    for ii = 1:(cmax-1)

        worker_grp_idx::Int64 = worker_nodes[ii].sector_ID
        workplace_idx::Int64 = worker_nodes[ii].workplace_ID

        for jj = (ii+1):cmax

            # If returned to work, WORK contacts consist of contacts within workplace,
            # + lower level contacts with workers of same worker group
            # Else, if WFH, no work contacts

            ### On days when at work, increased contacts with others at work  ###
            ## For workers in the same group and workplace, edges form according to ER graph
            if (worker_nodes[jj].sector_ID == worker_grp_idx) & (worker_nodes[jj].workplace_ID == workplace_idx)
                if rand(rng) < dd_within_workplace[worker_grp_idx]/(workplace_sizes[worker_grp_idx][workplace_idx] - 1) ## ER component
                    # Assign IDs of contact to tuple vectors
                    push!(work_contacts_same_workplace[ii],jj)
                    push!(work_contacts_same_workplace[jj],ii)

                    # Increment number of contacts for nodes ii & jj
                    work_contacts_same_workplace_per_node[ii] += 1
                    work_contacts_same_workplace_per_node[jj] += 1
                end

            ## same worker group, different workplace
            elseif (worker_nodes[jj].sector_ID == worker_grp_idx)

                if rand(rng) < prob_workertype_contact[worker_grp_idx]
                    # Assign IDs of contact to tuple vectors
                    push!(work_contacts_other_workplace[ii],jj)
                    push!(work_contacts_other_workplace[jj],ii)

                    # Increment number of contacts for nodes ii & jj
                    work_contacts_other_workplace_per_node[ii] += 1
                    work_contacts_other_workplace_per_node[jj] += 1
                end

            ## different worker group
            else
                if rand(rng) < prob_anyworker_contact
                    # Assign IDs of contact to tuple vectors
                    push!(work_contacts_other_workplace[ii],jj)
                    push!(work_contacts_other_workplace[jj],ii)

                    # Increment number of contacts for nodes ii & jj
                    work_contacts_other_workplace_per_node[ii] += 1
                    work_contacts_other_workplace_per_node[jj] += 1
                end
            end
        end
    end

end

"""
    generate_dynamic_worker_contacts_configuration(args)

Generates daily dynamic workplace contacts, from day 1 to endtime, according to specified sector-specific degree distributions and upper limit (default = 100).

Contacts are generated by iterating through each node and generating a random number of contacts from the appropriate distribution. Contacts can be made with any other node, except oneself.

Inputs: `rng` - random number generator,
        `cmax` - number of nodes in network,
        `endtime` - length of simulation,
        `worker_nodes` - array of WorkerParameters structures containing worker info,
        `dynamic_work_dd` - desired degree distribution for each sector,
        `max_contacts_work_dynamic` - maximum number of dynamic work contacts allowed \n
Outputs: `dynamic_worker_contacts` - 2D array{array} of contacted node IDs, for each day x each node \n
Location: network_generation_fns.jl
"""
function generate_dynamic_worker_contacts_configuration(rng::MersenneTwister,
                                            cmax::Int64,
                                            endtime::Int64,
                                            worker_nodes::Array{WorkerParameters,1},
                                            dynamic_work_dd::Array{Distribution,1},
                                            max_contacts_work_dynamic::Int64)

    """
    Initialise vector of vectors storing IDs of contacts for each node
    """
    dynamic_worker_contacts = Array{Array{Int64,1},2}(undef,endtime,cmax)

    """
    Iterate over all nodes
    For those returning to work in role with dynamic contacts,
    assign dynamic contacts for each timestep
    """
    for node_itr = 1:cmax
        # Get dynamic worker group type for node_itr
        node_dynamic_grp_ID = worker_nodes[node_itr].sector_ID

        for time_itr = 1:endtime
            # Generate number of dynamic contacts from appropriate distribution
            gg = round(Int64, Distributions.rand(rng, dynamic_work_dd[node_dynamic_grp_ID]))
            # Limit to specified maximum
            gg = min(gg, max_contacts_work_dynamic)
            # If dynamic worker contacts made, assign to output variable
            if gg > 0
                dynamic_worker_contacts[time_itr,node_itr] = zeros(gg)

                # Generate required number of contacts
                for contact_itr = 1:gg
                    gg1 = ceil(Int64,rand(rng)*cmax) # Get IDs of nodes connected by dynamic link on current timestep
                    while gg1 == node_itr # Redraw if returned the index node
                        gg1 = ceil(Int64,rand(rng)*cmax) # Get IDs of nodes connected by dynamic link on current timestep
                    end
                    dynamic_worker_contacts[time_itr,node_itr][contact_itr] = gg1
                end
            else # No dynamic worker contacts made on given timestep
                dynamic_worker_contacts[time_itr,node_itr] = Int64[]
            end
        end
    end

    return dynamic_worker_contacts::Array{Array{Int64,1},2}

end

"""
    find_fof(args)

Find and return all contacts of a given node's contacts (e.g. 'friends-of-friends')

Inputs: `contacts` - all contacts occurring for each node each day,
        `social_contacts_of_node` - contacts of given node to find contacts of,
        `edges_remaining_per_node` - edges remaining per node, according to desired degree distribution,
        `node_id1` - ID of current node,
        `fof` - array of node IDs of 'friends-of-friends' \n
Outputs: `fof` - array of node IDs of 'friends-of-friends' \n
Location: network_generation_fns.jl
"""
function find_fof(contacts::Array{Array{Int64,1},1},social_contacts_of_node::Array{Int64,1},
    edges_remaining_per_node::Array{Float64,1},node_id1::Int64,fof::Array{Int64,1})

    # # Find groups of friends already meeting today with edges remaining
    fof = collect(Iterators.flatten(contacts[social_contacts_of_node]))

    # Remove self and anyone already friends with and those without edges remaining
    # This version doesn't keep duplicates, so is no more likely to make an edge
    # with a node that is the friend of multiple friends
    # makes fewer allocations and is a lot faster
    setdiff!(fof,social_contacts_of_node)
    setdiff!(fof,[node_id1])
    removals = findall(edges_remaining_per_node[fof].==0)
    deleteat!(fof,removals)

    return fof
end

"""
    generate_social_contacts_each_day!(args)

Calls appropriate function to generate daily, dynamic social contacts, from starttime to endtime, for both work and non-work days.

Daily contacts are stored in the ContactStructure structure, under:
- `workday_social_contacts_by_day`: 2D array{array} of node IDs, for daily workday contacts per day x node
- `nonworkday_social_contacts_by_day`: 2D array{array} of node IDs, for daily non-workday contacts per day x node

Function called is specified by network_generation_method_dynamic_social in the NetworkParameters structure. Options are:
- `"ER"`: (BROKEN) for each node each day, a Poisson distributed number of contacts are generated from their friend group, with repetition
- `"configuration"`: (SLOW) for each node each day, generate contacts from social group using configuration model style algorithm with specified degree distribution. Clustering can be switched on (default) or off
- `"cluster"`: (DEFAULT) for each day, randomly cluster friend groups into fully connected subgroups, with sizes drawn from specified degree distribution. Each node is a member of one subgroup per day
- `"repeated"`: similar to "cluster", but contacts are repeated for dynamic_time_frame (default = 1) number of days
- `"groups"`: each day, all nodes (no social group structure) are split into a random number of groups of fixed size. Days can be repeated
- `"fixed_daily_degree"`: each day, all nodes (no social group structure) are split into groups according to a fixed group size and fixed number of daily contacts. Days can be repeated
- `"household_groups"`: each day, groups are made up of a specified number of households. Size of group is uniformly sampled between 0 and specified max. Days can be repeated

Inputs: `rng` - random number generator,
        `starttime` - time from which to generate daily social contacts,
        `... parameter structures ...`,
        `intervention_set_itr` - number of current intervention set \n
Outputs: None \n
Location: network_generation_fns.jl
"""
function generate_social_contacts_each_day!(rng::MersenneTwister,
                                            starttime::Int64,
                                            network_params::NetworkParameters,
                                            contacts::ContactStructure,
                                            configuration_options::ConfigurationOptions,
                                            intervention_set_itr::Int64)

    @unpack cmax, endtime = configuration_options
    @unpack network_generation_method_dynamic_social, social_workday_dd,
            social_nonworkday_dd, cluster_social_contacts, group_limit,
            dynamic_time_frame, n_groups_per_day_distribution,
            n_households_per_group = network_params
    @unpack social_contacts, social_contacts_per_node = contacts

    if length(contacts.workday_social_contacts_by_day) == 0
        contacts.workday_social_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)
        contacts.nonworkday_social_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)
    end

    if network_generation_method_dynamic_social == "ER"
        workday_social_contacts_by_day,
        nonworkday_social_contacts_by_day =  generate_social_contacts_each_day_ER(rng,
                                                                            cmax,
                                                                            starttime,
                                                                            endtime,
                                                                            social_contacts,
                                                                            social_contacts_per_node,
                                                                            n_social_mean_workday,
                                                                            n_social_mean_nonworkday)
    elseif network_generation_method_dynamic_social == "configuration"
        workday_social_contacts_by_day,
        nonworkday_social_contacts_by_day =  generate_social_contacts_each_day_configuration(rng,
                                                                                cmax,
                                                                                starttime,
                                                                                endtime,
                                                                                social_contacts,
                                                                                social_contacts_per_node,
                                                                                social_workday_dd,
                                                                                social_nonworkday_dd,
                                                                                cluster_social_contacts)

    elseif network_generation_method_dynamic_social == "cluster"
        workday_social_contacts_by_day,
        nonworkday_social_contacts_by_day =  generate_social_contacts_each_day_cluster(rng,
                                                                                cmax,
                                                                                starttime,
                                                                                endtime,
                                                                                social_contacts,
                                                                                social_workday_dd,
                                                                                social_nonworkday_dd)

    elseif network_generation_method_dynamic_social == "repeated"
        workday_social_contacts_by_day,
        nonworkday_social_contacts_by_day =  generate_social_contacts_each_day_repeated(rng,
                                                                                cmax,
                                                                                starttime,
                                                                                endtime,
                                                                                social_contacts,
                                                                                social_workday_dd,
                                                                                group_limit,
                                                                                dynamic_time_frame)

    elseif network_generation_method_dynamic_social == "groups"
        workday_social_contacts_by_day,
        nonworkday_social_contacts_by_day =  generate_social_contacts_each_day_groups(rng,
                                                                                cmax,
                                                                                starttime,
                                                                                endtime,
                                                                                n_groups_per_day_distribution,
                                                                                group_limit,
                                                                                dynamic_time_frame)
        contacts.social_contacts_per_node .+= 1

    elseif network_generation_method_dynamic_social == "fixed_daily_degree"
        workday_social_contacts_by_day,
        nonworkday_social_contacts_by_day =  generate_social_contacts_each_day_fixed_daily_degree(rng,
                                                                                cmax,
                                                                                starttime,
                                                                                endtime,
                                                                                12,
                                                                                group_limit,
                                                                                dynamic_time_frame)
        contacts.social_contacts_per_node .+= 1
    elseif network_generation_method_dynamic_social == "household_groups"
        workday_social_contacts_by_day,
        nonworkday_social_contacts_by_day =  generate_social_contacts_each_day_household_groups_stochastic(rng,
                                                                                        cmax,
                                                                                        starttime,
                                                                                        endtime,
                                                                                        group_limit,
                                                                                        1,
                                                                                        dynamic_time_frame,
                                                                                        n_households_per_group,
                                                                                        contacts,
                                                                                        intervention_set_itr)
        contacts.social_contacts_per_node .+= 1
    else
        println("Invalid dynamic social contact generation method.")
    end

    contacts.workday_social_contacts_by_day[starttime:endtime, :] = workday_social_contacts_by_day[starttime:endtime, :]
    contacts.nonworkday_social_contacts_by_day[starttime:endtime, :] = nonworkday_social_contacts_by_day[starttime:endtime, :]

    return nothing
end

"""
    generate_social_contacts_each_day_configuration(args)

(SLOW) Generate daily, dynamic social contacts, from starttime to endtime, for both work and non-work days, according to specified degree distribution.

For each node each day, generate contacts from social group using configuration model style algorithm with specified degree distribution. Clustering can be switched on (default) or off

Inputs: `rng` - random number generator,
        `cmax` - number of nodes in network,
        `starttime` - time from which to generate daily social contacts,
        `endtime` - length of simulation,
        `social_contacts` - array{array} of fixed social contacts per node,
        `social_contacts_per_node` - array of number of fixed social contacts per node,
        `social_workday_dd` - desired degree distribution of social contacts on workdays,
        `social_nonworkday_dd` - desired degree distribution of social contacts on non-workdays,
        `cluster_social_contacts` - boolean flaggin whether or not to cluster social contacts using f-o-f \n
Outputs: `workday_social_contacts_by_day` - 2D array{array} of social contacts made per day per node (workday),
         `nonworkday_social_contacts_by_day` - 2D array{array} of social contacts made per day per node (non-workday) \n
Location: network_generation_fns.jl
"""
function generate_social_contacts_each_day_configuration(rng::MersenneTwister,
                                            cmax::Int64,
                                            starttime::Int64,
                                            endtime::Int64,
                                            social_contacts::Array{Array{Int64,1},1},
                                            social_contacts_per_node::Array{Int64,1},
                                            social_workday_dd::Distribution,
                                            social_nonworkday_dd::Distribution,
                                            cluster_social_contacts::Bool)

    """
    Initialise vector of vectors storing IDs of contacts for each node
    """
    workday_social_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)
    nonworkday_social_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)

    """
    Configuration model
    """
    n_nodes = cmax

    edges_per_node_workday = Distributions.rand(rng, social_workday_dd, n_nodes)
    edges_per_node_nonworkday = Distributions.rand(rng, social_workday_dd, n_nodes)

    # Round degree to nearest whole number
    # Limit maximum number of contacts to size of friendship group
    for node_id = 1:n_nodes
        friendship_group_size = social_contacts_per_node[node_id]
        edges_per_node_workday[node_id] = round(Int64, edges_per_node_workday[node_id])
        edges_per_node_workday[node_id] = min(friendship_group_size, edges_per_node_workday[node_id])
        edges_per_node_nonworkday[node_id] = round(Int64, edges_per_node_nonworkday[node_id])
        edges_per_node_nonworkday[node_id] = min(friendship_group_size, edges_per_node_nonworkday[node_id])
        # Initialise array to store daily contacts
        for time_itr = 1:endtime
            workday_social_contacts_by_day[time_itr,node_id] = Int64[]
            nonworkday_social_contacts_by_day[time_itr,node_id] = Int64[]
        end
    end

    fof = Int64[]
    for time_itr = starttime:endtime

        """
        WORKDAY
        """
        half_edges = cumsum(edges_per_node_workday)

        edges_remaining_per_node = copy(edges_per_node_workday)
        while half_edges[end] > 1

            stub_id1 = round(Int64, rand(rng)*(half_edges[end]-1) + 1)
            node_id1 = findfirst(x -> x >= stub_id1, half_edges)

            if cluster_social_contacts == true
                # Find groups of friends already meeting today with edges remaining
                # fof = collect(Iterators.flatten(workday_social_contacts_by_day[time_itr,social_contacts[node_id1]]))
                fof = find_fof(workday_social_contacts_by_day,social_contacts[node_id1],edges_remaining_per_node,node_id1,fof)
            else
                fof = Int64[]
            end

            # Contact made with group of friends already meeting
            if length(fof) > 0

                node_id2 = fof[ceil(Int64, rand(rng)*length(fof))]

                push!(workday_social_contacts_by_day[time_itr,node_id1], node_id2)
                push!(workday_social_contacts_by_day[time_itr,node_id2], node_id1)

                edges_remaining_per_node[node_id1] -= 1
                edges_remaining_per_node[node_id2] -= 1


                for ii=node_id1:n_nodes
                    half_edges[ii] -= 1
                end

                for ii=node_id2:n_nodes
                    half_edges[ii] -= 1
                end

            # Contact made with anyone from friend group
            else
                friend_group = social_contacts[node_id1]

                # Remove anyone already meeting today
                friend_group = friend_group[[friend_group[i] ∉ workday_social_contacts_by_day[time_itr,node_id1] for i=1:length(friend_group)]]
                # Remove those with no half edges remaining
                friend_group = friend_group[edges_remaining_per_node[friend_group].>0]

                # If possible, connect to random friend
                if length(friend_group) > 0
                    node_id2 = friend_group[ceil(Int64, rand(rng)*length(friend_group))]

                    push!(workday_social_contacts_by_day[time_itr,node_id1], node_id2)
                    push!(workday_social_contacts_by_day[time_itr,node_id2], node_id1)

                    edges_remaining_per_node[node_id1] -= 1
                    edges_remaining_per_node[node_id2] -= 1


                    for ii=node_id1:n_nodes
                        half_edges[ii] -= 1
                    end

                    for ii=node_id2:n_nodes
                        half_edges[ii] -= 1
                    end

                # Otherwise discard half-edge
                else
                    edges_remaining_per_node[node_id1] -= 1

                    for ii=node_id1:n_nodes
                        half_edges[ii] -= 1
                    end
                end
            end
        end

        """
        NON-WORKDAY
        """
        half_edges = cumsum(edges_per_node_nonworkday)

        edges_remaining_per_node = copy(edges_per_node_nonworkday)
        while half_edges[end] > 1

            stub_id1 = round(Int64, rand(rng)*(half_edges[end]-1) + 1)
            node_id1 = findfirst(x -> x >= stub_id1, half_edges)

            # fof = collect(Iterators.flatten(nonworkday_social_contacts_by_day[time_itr,social_contacts[node_id1]]))
            fof = find_fof(view(nonworkday_social_contacts_by_day[time_itr,:]),social_contacts[node_id1],edges_remaining_per_node,node_id1,fof)

            # Contact made with group of friends already meeting
            if (length(fof) > 0) & (cluster_social_contacts == true)

                node_id2 = fof[ceil(Int64, rand(rng)*length(fof))]

                push!(nonworkday_social_contacts_by_day[time_itr,node_id1], node_id2)
                push!(nonworkday_social_contacts_by_day[time_itr,node_id2], node_id1)

                edges_remaining_per_node[node_id1] -= 1
                edges_remaining_per_node[node_id2] -= 1


                for ii=node_id1:n_nodes
                    half_edges[ii] -= 1
                end

                for ii=node_id2:n_nodes
                    half_edges[ii] -= 1
                end

            # Contact made with anyone from friend group
            else
                friend_group = social_contacts[node_id1]

                # Remove anyone already meeting today
                friend_group = friend_group[[friend_group[i] ∉ nonworkday_social_contacts_by_day[time_itr,node_id1] for i=1:length(friend_group)]]
                # Remove those with no half edges remaining
                friend_group = friend_group[edges_remaining_per_node[friend_group].>0]

                # If possible, connect to random friend
                if length(friend_group) > 0
                    node_id2 = friend_group[ceil(Int64, rand(rng)*length(friend_group))]

                    push!(nonworkday_social_contacts_by_day[time_itr,node_id1], node_id2)
                    push!(nonworkday_social_contacts_by_day[time_itr,node_id2], node_id1)

                    edges_remaining_per_node[node_id1] -= 1
                    edges_remaining_per_node[node_id2] -= 1


                    for ii=node_id1:n_nodes
                        half_edges[ii] -= 1
                    end

                    for ii=node_id2:n_nodes
                        half_edges[ii] -= 1
                    end

                # Otherwise discard half-edge
                else
                    edges_remaining_per_node[node_id1] -= 1

                    for ii=node_id1:n_nodes
                        half_edges[ii] -= 1
                    end
                end
            end
        end
    end

    return workday_social_contacts_by_day::Array{Array{Int64,1},2},
            nonworkday_social_contacts_by_day::Array{Array{Int64,1},2}
end

"""
    generate_social_contacts_each_day_cluster(args)

Generate daily, dynamic social contacts, from starttime to endtime, for both work and non-work days, according to specified degree distribution.

For each day, randomly cluster friend groups into fully connected subgroups, with sizes drawn from specified degree distribution. Each node is a member of one subgroup per day.

Inputs: `rng` - random number generator,
        `cmax` - number of nodes in network,
        `starttime` - time from which to generate daily social contacts,
        `endtime` - length of simulation,
        `social_contacts` - array{array} of fixed social contacts per node,
        `social_workday_dd` - desired degree distribution of social contacts on workdays,
        `social_nonworkday_dd` - desired degree distribution of social contacts on non-workdays,
Outputs: `workday_social_contacts_by_day` - 2D array{array} of social contacts made per day per node (workday),
         `nonworkday_social_contacts_by_day` - 2D array{array} of social contacts made per day per node (non-workday) \n
Location: network_generation_fns.jl
"""
function generate_social_contacts_each_day_cluster(rng::MersenneTwister,
                                            cmax::Int64,
                                            starttime::Int64,
                                            endtime::Int64,
                                            social_contacts::Array{Array{Int64,1},1},
                                            social_workday_dd::Distribution,
                                            social_nonworkday_dd::Distribution)

    """
    Initialise vector of vectors storing IDs of contacts for each node
    """
    workday_social_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)
    nonworkday_social_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)

    for time_itr = starttime:endtime
        # println(time_itr)
        """
        WORKDAY
        """
        # Initialise array to track who has already been contacted
        edges_unassigned = ones(Int64, cmax)

        # Iteratively create clusters of social contacts
        while sum(edges_unassigned) > 0

            # how many nodes are unassigned
            num_unassigned = sum(edges_unassigned)

            # Choose unassigned node at random
            unassigned_id = rand(rng,1:num_unassigned)

            # find the associated node
            ii = 0
            node_id = 0
            while ii<unassigned_id
                node_id += 1
                if edges_unassigned[node_id]==1
                    ii +=1
                end
            end

            # Find friends of chosen node
            friend_list = social_contacts[node_id]

            # Remove friends who are already assigned
            friend_list = friend_list[edges_unassigned[friend_list].==1]

            # Draw random cluster size
            cluster_size = Distributions.rand(rng, social_workday_dd)

            # Round and limit to number of unassigned friends
            cluster_size = round(Int64, min(cluster_size, length(friend_list)))

            if cluster_size > 0
                # Choose friends in cluster at random
                shuffle!(rng, friend_list)
                friends_in_cluster = friend_list[1:cluster_size]

                # Record contacts and assign nodes
                workday_social_contacts_by_day[time_itr,node_id] = copy(friends_in_cluster)
                for friend_itr = 1:cluster_size
                    # ID of current friend node
                    friend_id = friends_in_cluster[friend_itr]

                    # for each friend in the cluster
                    workday_social_contacts_by_day[time_itr,friend_id] = zeros(Int64,cluster_size)
                    # find other friends not in the cluster
                    other_friends = social_contacts[friend_id][[social_contacts[friend_id][i] ∉ friends_in_cluster for i=1:length(social_contacts[friend_id])]]
                    shuffle!(rng,other_friends)
                    jj = 1
                    for ii=1:cluster_size
                        # if they are in the social contacts of friend_id, then take them
                        if (ii!=friend_itr) && (friends_in_cluster[ii] ∈ social_contacts[friend_id])
                            workday_social_contacts_by_day[time_itr,friend_id][ii] = friends_in_cluster[ii]
                        elseif jj <= length(other_friends)
                            # otherwise pick another of their social contacts at random
                            workday_social_contacts_by_day[time_itr,friend_id][ii] = other_friends[jj]
                            jj += 1
                        end
                    end
                    if jj>length(other_friends)
                        workday_social_contacts_by_day[time_itr,friend_id] = workday_social_contacts_by_day[time_itr,friend_id][workday_social_contacts_by_day[time_itr,friend_id].!=0]
                    end

                    # # Calculate shortage of contacts caused by non-overlapping friend groups
                    # spots_remaining = cluster_size - (ii-1)
                    #
                    # # friends_of_friends_in_cluster = friends_in_cluster[[friends_in_cluster[i] ∈ social_contacts[friend_id] for i=1:cluster_size]]
                    # # # Assign intersection to today's contacts
                    # # # if length(friends_of_friends_in_cluster) > 0
                    # #     workday_social_contacts_by_day[time_itr,friend_id] = friends_of_friends_in_cluster
                    # # # else
                    # # #     workday_social_contacts_by_day[time_itr,friend_id] = Int64[]
                    # # # end
                    # # # Calculate shortage of contacts caused by non-overlapping friend groups
                    # # spots_remaining = cluster_size - length(friends_of_friends_in_cluster)
                    # # Assign remaining spots to other friends at random
                    # other_friends = social_contacts[friend_id][[social_contacts[friend_id][i] ∉ friends_in_cluster for i=1:length(social_contacts[friend_id])]]
                    # spots_remaining = min(spots_remaining, length(other_friends))
                    # shuffle!(rng,other_friends)
                    # for spot_itr = 1:spots_remaining
                    #     push!(workday_social_contacts_by_day[time_itr,friend_id], other_friends[spot_itr])
                    # end
                end

                edges_unassigned[friends_in_cluster] .= 0

            else
                workday_social_contacts_by_day[time_itr,node_id] = Int64[]
            end
            edges_unassigned[node_id] = 0
        end

        """
        NON-WORKDAY
        """
        # Initialise array to track who has already been contacted
        edges_unassigned = ones(Int64, cmax)

        # Iteratively create clusters of social contacts
        while sum(edges_unassigned) >0

            # how many nodes are unassigned
            num_unassigned = sum(edges_unassigned)

            # Choose unassigned node at random
            unassigned_id = rand(rng,1:num_unassigned)

            # find the associated node
            ii = 0
            node_id = 0
            while ii<unassigned_id
                node_id += 1
                if edges_unassigned[node_id]==1
                    ii +=1
                end
            end

            # Find friends of chosen node
            friend_list = social_contacts[node_id]

            # Remove friends who are already assigned
            friend_list = friend_list[edges_unassigned[friend_list].==1]

            # Draw random cluster size
            cluster_size = Distributions.rand(rng, social_nonworkday_dd)

            # Round and limit to number of unassigned friends
            cluster_size = round(Int64, min(cluster_size, length(friend_list)))

            if cluster_size > 0
                # Choose friends in cluster at random
                shuffle!(rng, friend_list)
                friends_in_cluster = friend_list[1:cluster_size]

                # Record contacts and assign nodes
                nonworkday_social_contacts_by_day[time_itr,node_id] = copy(friends_in_cluster)
                for friend_itr = 1:cluster_size
                    # ID of current friend node
                    friend_id = friends_in_cluster[friend_itr]

                    # for each friend in the cluster
                    nonworkday_social_contacts_by_day[time_itr,friend_id] = zeros(Int64,cluster_size)
                    # find other friends not in the cluster
                    other_friends = social_contacts[friend_id][[social_contacts[friend_id][i] ∉ friends_in_cluster for i=1:length(social_contacts[friend_id])]]
                    shuffle!(rng,other_friends)
                    jj = 1
                    for ii=1:cluster_size
                        # if they are in the social contacts of friend_id, then take them
                        if (ii!=friend_itr) && (friends_in_cluster[ii] ∈ social_contacts[friend_id])
                            nonworkday_social_contacts_by_day[time_itr,friend_id][ii] = friends_in_cluster[ii]
                        elseif jj <= length(other_friends)
                            # otherwise pick another of their social contacts at random
                            nonworkday_social_contacts_by_day[time_itr,friend_id][ii] = other_friends[jj]
                            jj += 1
                        end
                    end
                    # if there weren't enough other friends to fill the contacts, remove extra zeros
                    if jj>length(other_friends)
                        nonworkday_social_contacts_by_day[time_itr,friend_id] = workday_social_contacts_by_day[time_itr,friend_id][workday_social_contacts_by_day[time_itr,friend_id].!=0]
                    end

                    # # Find intersection of friend group with current cluster
                    # friends_of_friends_in_cluster = friends_in_cluster[[friends_in_cluster[i] ∈ social_contacts[friend_id] for i=1:cluster_size]]
                    # # Assign intersection to today's contacts
                    # if length(friends_of_friends_in_cluster) > 0
                    #     nonworkday_social_contacts_by_day[time_itr,friend_id] = friends_of_friends_in_cluster
                    # else
                    #     nonworkday_social_contacts_by_day[time_itr,friend_id] = Int64[]
                    # end
                    # # Calculate shortage of contacts caused by non-overlapping friend groups
                    # spots_remaining = cluster_size - length(friends_of_friends_in_cluster)
                    # # Assign remaining spots to other friends at random
                    # other_friends = social_contacts[friend_id][[social_contacts[friend_id][i] ∉ friends_in_cluster for i=1:length(social_contacts[friend_id])]]
                    # spots_remaining = min(spots_remaining, length(other_friends))
                    # shuffle!(other_friends)
                    # for spot_itr = 1:spots_remaining
                    #     push!(nonworkday_social_contacts_by_day[time_itr,friend_id], other_friends[spot_itr])
                    # end
                end

                edges_unassigned[friends_in_cluster] .= 0

            else
                nonworkday_social_contacts_by_day[time_itr,node_id] = Int64[]
            end
            edges_unassigned[node_id] = 0
        end
    end

    return workday_social_contacts_by_day::Array{Array{Int64,1},2},
            nonworkday_social_contacts_by_day::Array{Array{Int64,1},2}
end

"""
    generate_social_contacts_each_day_repeated(args)

Generate daily, dynamic social contacts, from starttime to endtime, for both work and non-work days, according to specified degree distribution.

Similar to "cluster", but contacts are repeated for dynamic_time_frame (default = 1) number of days.
Used to assess role of static vs dynamic social contacts.

Inputs: `rng` - random number generator,
        `cmax` - number of nodes in network,
        `starttime` - time from which to generate daily social contacts,
        `endtime` - length of simulation,
        `social_contacts` - array{array} of fixed social contacts per node,
        `social_workday_dd` - desired degree distribution of social contacts on workdays,
        `group_limit` - maximum number of people allowed to meet in one group,
        `dynamic_time_frame` - number of days to repeat contacts for \n
Outputs: `workday_social_contacts_by_day` - 2D array{array} of social contacts made per day per node (workday),
         `nonworkday_social_contacts_by_day` - 2D array{array} of social contacts made per day per node (non-workday) \n
Location: network_generation_fns.jl
"""
function generate_social_contacts_each_day_repeated(rng::MersenneTwister,
                                            cmax::Int64,
                                            starttime::Int64,
                                            endtime::Int64,
                                            social_contacts::Array{Array{Int64,1},1},
                                            social_workday_dd::Distribution,
                                            group_limit::Int64,
                                            dynamic_time_frame::Int64)

    """
    Initialise vector of vectors storing IDs of contacts for each node
    """
    workday_social_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)

    time_itr = starttime
    while time_itr <= endtime
        # println(time_itr)
        """
        WORKDAY
        """
        # Initialise array to track who has already been contacted
        edges_unassigned = ones(Int64, cmax)

        # Iteratively create clusters of social contacts
        while sum(edges_unassigned) > 0

            # how many nodes are unassigned
            num_unassigned = sum(edges_unassigned)

            # Choose unassigned node at random
            unassigned_id = rand(rng,1:num_unassigned)

            # find the associated node
            ii = 0
            node_id = 0
            while ii<unassigned_id
                node_id += 1
                if edges_unassigned[node_id]==1
                    ii +=1
                end
            end

            # Find friends of chosen node
            friend_list = social_contacts[node_id]

            # Remove friends who are already assigned
            friend_list = friend_list[edges_unassigned[friend_list].==1]

            # Draw random cluster size
            cluster_size = Distributions.rand(rng, social_workday_dd)

            # Round and limit to number of unassigned friends or group limit
            cluster_size = round(Int64, min(cluster_size, length(friend_list), group_limit))

            if cluster_size > 0
                # Choose friends in cluster at random
                shuffle!(rng, friend_list)
                friends_in_cluster = friend_list[1:cluster_size]

                # Record contacts and assign nodes
                workday_social_contacts_by_day[time_itr,node_id] = copy(friends_in_cluster)
                for friend_itr = 1:cluster_size
                    # ID of current friend node
                    friend_id = friends_in_cluster[friend_itr]

                    # for each friend in the cluster
                    workday_social_contacts_by_day[time_itr,friend_id] = zeros(Int64,cluster_size)
                    # find other friends not in the cluster
                    other_friends = social_contacts[friend_id][[social_contacts[friend_id][i] ∉ friends_in_cluster for i=1:length(social_contacts[friend_id])]]
                    shuffle!(rng,other_friends)
                    jj = 1
                    for ii=1:cluster_size
                        # if they are in the social contacts of friend_id, then take them
                        if (ii!=friend_itr) && (friends_in_cluster[ii] ∈ social_contacts[friend_id])
                            workday_social_contacts_by_day[time_itr,friend_id][ii] = friends_in_cluster[ii]
                        elseif jj <= length(other_friends)
                            # otherwise pick another of their social contacts at random
                            workday_social_contacts_by_day[time_itr,friend_id][ii] = other_friends[jj]
                            jj += 1
                        end
                    end
                    if jj>length(other_friends)
                        workday_social_contacts_by_day[time_itr,friend_id] = workday_social_contacts_by_day[time_itr,friend_id][workday_social_contacts_by_day[time_itr,friend_id].!=0]
                    end
                end

                edges_unassigned[friends_in_cluster] .= 0

            else
                workday_social_contacts_by_day[time_itr,node_id] = Int64[]
            end
            edges_unassigned[node_id] = 0
        end

        if (dynamic_time_frame > 1) && (time_itr < endtime)
            days_left = endtime - time_itr
            if days_left >= dynamic_time_frame
                for time_rep = (time_itr+1):(time_itr+dynamic_time_frame-1)
                    workday_social_contacts_by_day[time_rep,:] = workday_social_contacts_by_day[time_itr,:]
                end
            else
                for time_rep = (time_itr+1):endtime
                    workday_social_contacts_by_day[time_rep,:] = workday_social_contacts_by_day[time_itr,:]
                end
            end
        end

        time_itr += dynamic_time_frame
    end

    nonworkday_social_contacts_by_day = workday_social_contacts_by_day

    return workday_social_contacts_by_day::Array{Array{Int64,1},2},
            nonworkday_social_contacts_by_day::Array{Array{Int64,1},2}
end


"""
    generate_social_contacts_each_day_groups(args)

Generate daily, dynamic social contacts, from starttime to endtime, for both work and non-work days, according to specified group size.

Each day, all nodes (no social group structure) are split into a random number of groups of fixed size. Days can be repeated.
Used for rule of 6 analysis

Inputs: `rng` - random number generator,
        `cmax` - number of nodes in network,
        `starttime` - time from which to generate daily social contacts,
        `endtime` - length of simulation,
        `n_groups_per_day_distribution` - distribution of number of meeting groups per day,
        `group_limit` - number of people in each group,
        `dynamic_time_frame` - number of days to repeat contacts for \n
Outputs: `workday_social_contacts_by_day` - 2D array{array} of social contacts made per day per node (workday),
         `nonworkday_social_contacts_by_day` - 2D array{array} of social contacts made per day per node (non-workday) \n
Location: network_generation_fns.jl
"""
function generate_social_contacts_each_day_groups(rng::MersenneTwister,
                                            cmax::Int64,
                                            starttime::Int64,
                                            endtime::Int64,
                                            n_groups_per_day_distribution::Distribution,
                                            group_limit::Int64,
                                            dynamic_time_frame::Int64)

    """
    Initialise vector of vectors storing IDs of contacts for each node
    """
    workday_social_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)

    time_itr = starttime
    while time_itr <= endtime

        n_groups_per_node = zeros(Int64, cmax)
        for node_itr = 1:cmax
            n_groups_per_node[node_itr] = Distributions.rand(rng, n_groups_per_day_distribution)
            workday_social_contacts_by_day[time_itr,node_itr] = Int64[]
        end

        while sum(n_groups_per_node) > 0

            nodes_remaining = findall(n_groups_per_node.>0)

            current_group = shuffle(nodes_remaining)[1:min(group_limit,length(nodes_remaining))]

            group_size = length(current_group)

            for ii = 1:(group_size-1)
                subject_id = current_group[ii]
                for jj = (ii+1):group_size
                    friend_id = current_group[jj]
                    if isassigned(workday_social_contacts_by_day,time_itr,subject_id)==false
                        workday_social_contacts_by_day[time_itr,subject_id] = Int64[]
                    end
                    if isassigned(workday_social_contacts_by_day,time_itr,friend_id)==false
                        workday_social_contacts_by_day[time_itr,friend_id] = Int64[]
                    end
                    push!(workday_social_contacts_by_day[time_itr,subject_id], friend_id)
                    push!(workday_social_contacts_by_day[time_itr,friend_id], subject_id)
                end
            end

            n_groups_per_node[current_group] .-= 1

        end

        if (dynamic_time_frame > 1) && (time_itr < endtime)
            days_left = endtime - time_itr
            if days_left >= dynamic_time_frame
                for time_rep = (time_itr+1):(time_itr+dynamic_time_frame-1)
                    workday_social_contacts_by_day[time_rep,:] = workday_social_contacts_by_day[time_itr,:]
                end
            else
                for time_rep = (time_itr+1):endtime
                    workday_social_contacts_by_day[time_rep,:] = workday_social_contacts_by_day[time_itr,:]
                end
            end
        end

        time_itr += dynamic_time_frame
    end

    nonworkday_social_contacts_by_day = workday_social_contacts_by_day

    return workday_social_contacts_by_day::Array{Array{Int64,1},2},
            nonworkday_social_contacts_by_day::Array{Array{Int64,1},2}
end


"""
    generate_social_contacts_each_day_fixed_daily_degree(args)

Generate daily, dynamic social contacts, from starttime to endtime, for both work and non-work days, according to specified daily degree.

Each day, all nodes (no social group structure) are split into groups according to a fixed group size and fixed number of daily contacts. Days can be repeated.
Used for rule of 6 analysis

Inputs: `rng` - random number generator,
        `cmax` - number of nodes in network,
        `starttime` - time from which to generate daily social contacts,
        `endtime` - length of simulation,
        `contacts_per_day` - number of social contacts made per day,
        `group_limit` - number of people in each group,
        `dynamic_time_frame` - number of days to repeat contacts for \n
Outputs: `workday_social_contacts_by_day` - 2D array{array} of social contacts made per day per node (workday),
         `nonworkday_social_contacts_by_day` - 2D array{array} of social contacts made per day per node (non-workday) \n
Location: network_generation_fns.jl
"""
function generate_social_contacts_each_day_fixed_daily_degree(rng::MersenneTwister,
                                            cmax::Int64,
                                            starttime::Int64,
                                            endtime::Int64,
                                            contacts_per_day::Int64,
                                            group_limit::Int64,
                                            dynamic_time_frame::Int64)

    """
    Initialise vector of vectors storing IDs of contacts for each node
    """
    workday_social_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)

    time_itr = starttime
    while time_itr <= endtime

        n_groups_per_node = ones(Int64, cmax)*round(Int64,contacts_per_day/group_limit)
        for node_itr = 1:cmax
            workday_social_contacts_by_day[time_itr,node_itr] = Int64[]
        end

        while sum(n_groups_per_node) > 0

            nodes_remaining = findall(n_groups_per_node.>0)

            current_group = shuffle(nodes_remaining)[1:min(group_limit,length(nodes_remaining))]

            group_size = length(current_group)

            for ii = 1:(group_size-1)
                subject_id = current_group[ii]
                for jj = (ii+1):group_size
                    friend_id = current_group[jj]
                    if isassigned(workday_social_contacts_by_day,time_itr,subject_id)==false
                        workday_social_contacts_by_day[time_itr,subject_id] = Int64[]
                    end
                    if isassigned(workday_social_contacts_by_day,time_itr,friend_id)==false
                        workday_social_contacts_by_day[time_itr,friend_id] = Int64[]
                    end
                    push!(workday_social_contacts_by_day[time_itr,subject_id], friend_id)
                    push!(workday_social_contacts_by_day[time_itr,friend_id], subject_id)
                end
            end

            n_groups_per_node[current_group] .-= 1

        end

        if (dynamic_time_frame > 1) && (time_itr < endtime)
            days_left = endtime - time_itr
            if days_left >= dynamic_time_frame
                for time_rep = (time_itr+1):(time_itr+dynamic_time_frame-1)
                    workday_social_contacts_by_day[time_rep,:] = workday_social_contacts_by_day[time_itr,:]
                end
            else
                for time_rep = (time_itr+1):endtime
                    workday_social_contacts_by_day[time_rep,:] = workday_social_contacts_by_day[time_itr,:]
                end
            end
        end

        time_itr += dynamic_time_frame
    end

    nonworkday_social_contacts_by_day = workday_social_contacts_by_day

    return workday_social_contacts_by_day::Array{Array{Int64,1},2},
            nonworkday_social_contacts_by_day::Array{Array{Int64,1},2}
end

"""
    generate_dynamic_worker_contacts_ER(args)

Generates daily dynamic workplace contacts, from day 1 to endtime, according to sector-specific Normal distribution.

Contacts are generated by iterating through each node and generating a random number of contacts from the appropriate distribution. Contacts can be made with any other node, except oneself.

Inputs: `rng` - random number generator,
        `cmax` - number of nodes in network,
        `endtime` - length of simulation,
        `worker_nodes` - array of WorkerParameters structures containing worker info,
        `dynamic_conts_mean` - mean number of dynamic contacts per day,
        `dynamic_conts_sd` - SD of number of dynamic contacts per day \n
Outputs: `dynamic_worker_contacts` - 2D array{array} of dynamic work contacts made per day per node \n
Location: network_generation_fns.jl
"""
function generate_dynamic_worker_contacts_ER(rng::MersenneTwister,
                                            cmax::Int64,
                                            endtime::Int64,
                                            worker_nodes::Array{WorkerParameters,1},
                                            dynamic_conts_mean::Array{Float64,1},
                                            dynamic_conts_sd::Array{Float64,1})

    """
    Initialise vector of vectors storing IDs of contacts for each node
    """
    dynamic_worker_contacts = Array{Array{Int64,1},2}(undef,endtime,cmax)

    """
    Iterate over all nodes
    For those returning to work in role with dynamic contacts,
    assign dynamic contacts for each timestep
    """
    for node_itr = 1:cmax
        # Get dynamic worker group type for node_itr
        node_dynamic_grp_ID = worker_nodes[node_itr].sector_ID

        for time_itr = 1:endtime
            # Generate number of dynamic contacts from appropriate distribution
            gg = round(Int64,abs(randn(rng)*dynamic_conts_sd[node_dynamic_grp_ID]+dynamic_conts_mean[node_dynamic_grp_ID]))

            # If dynamic worker contacts made, assign to output variable
            if gg > 0
                dynamic_worker_contacts[time_itr,node_itr] = zeros(gg)

                # Generate required number of contacts
                for contact_itr = 1:gg
                    gg1 = ceil(Int64,rand(rng)*cmax) # Get IDs of nodes connected by dynamic link on current timestep
                    while gg1 == node_itr # Redraw if returned the index node
                        gg1 = ceil(Int64,rand(rng)*cmax) # Get IDs of nodes connected by dynamic link on current timestep
                    end
                    dynamic_worker_contacts[time_itr,node_itr][contact_itr] = gg1
                end
            else # No dynamic worker contacts made on given timestep
                dynamic_worker_contacts[time_itr,node_itr] = Int64[]
            end
        end
    end

    return dynamic_worker_contacts::Array{Array{Int64,1},2}

end

"""
    generate_social_contacts_each_day_ER(args)

(BROKEN) Generate daily, dynamic social contacts, from starttime to endtime, for both work and non-work days, according to specified mean degree.

For each node each day, a Poisson distributed number of contacts are generated from their friend group, with repetition.

Inputs: `rng` - random number generator,
        `cmax` - number of nodes in network,
        `starttime` - time from which to generate daily social contacts,
        `endtime` - length of simulation,
        `social_contacts` - array{array} of social contacts per node,
        `social_contacts_per_node` - array of number of social contacts per node,
        `n_social_mean_workday` - mean number of social contacts per workday,
        `n_social_mean_workday` - mean number of social contacts per non-workday \n
Outputs: `workday_social_contacts_by_day` - 2D array{array} of social contacts made per day per node (workday),
         `nonworkday_social_contacts_by_day` - 2D array{array} of social contacts made per day per node (non-workday) \n
Location: network_generation_fns.jl
"""
function generate_social_contacts_each_day_ER(rng::MersenneTwister,
                                            cmax::Int64,
                                            starttime::Int64,
                                            endtime::Int64,
                                            social_contacts::Array{Array{Int64,1},1},
                                            social_contacts_per_node::Array{Int64,1},
                                            n_social_mean_workday::Int64,
                                            n_social_mean_nonworkday::Int64)

    """
    Initialise vector of vectors storing IDs of contacts for each node
    """
    workday_social_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)
    nonworkday_social_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)

    """
    Iterate over all nodes
    For those potentially have social contacts,
    assign a sample of those contacts for each timestep
    """
    for node_itr = 1:cmax

        # Add in social contacts if possible
        if social_contacts_per_node[node_itr] > 0

            # Social contacts may differ each day
            for time_itr = 1:endtime

                # Draw random number of social links to be made today
                n_social_workday = min(rand(rng,Poisson(n_social_mean_workday)), social_contacts_per_node[node_itr])
                n_social_nonworkday = min(rand(rng,Poisson(n_social_mean_nonworkday)), social_contacts_per_node[node_itr])

                if n_social_workday >0
                    workday_social_contacts_by_day[time_itr,node_itr] = zeros(n_social_workday)

                    for wd_social_it = 1:n_social_workday
                        wd_contact_idx = ceil(Int64, rand(rng)*social_contacts_per_node[node_itr])
                        workday_social_contacts_by_day[time_itr,node_itr][wd_social_it] = social_contacts[node_itr][wd_contact_idx]
                    end
                else
                    workday_social_contacts_by_day[time_itr,node_itr] = Int64[]
                end

                if n_social_nonworkday >0
                    nonworkday_social_contacts_by_day[time_itr,node_itr] = zeros(n_social_nonworkday)

                    for nwd_social_it = 1:n_social_nonworkday
                        nwd_contact_idx = ceil(Int64, rand(rng)*social_contacts_per_node[node_itr])
                        nonworkday_social_contacts_by_day[time_itr,node_itr][nwd_social_it] = social_contacts[node_itr][nwd_contact_idx]
                    end
                else
                    nonworkday_social_contacts_by_day[time_itr,node_itr] = Int64[]
                end
            end
        else
            for time_itr = 1:endtime
                workday_social_contacts_by_day[time_itr,node_itr] = Int64[]
                nonworkday_social_contacts_by_day[time_itr,node_itr] = Int64[]
            end
        end
    end

    return workday_social_contacts_by_day::Array{Array{Int64,1},2},
            nonworkday_social_contacts_by_day::Array{Array{Int64,1},2}
end

"""
    generate_random_contacts(args)

Generates daily random contacts, from day 1 to endtime, according to fixed probability of contact occurring between any two nodes (Erdos-Reyni).

Inputs: `rng` - random number generator,
        `cmax` - number of nodes in network,
        `endtime` - length of simulation,
        `prob_random_contact` - fixed probability of each contact occurring \n
Outputs: `random_contacts_by_day` - 2D array{array} of random contacts made per day per node \n
Location: network_generation_fns.jl
"""
function generate_random_contacts(rng::MersenneTwister,
                                    cmax::Int64,
                                    endtime::Int64,
                                    prob_random_contact::Float64)

    random_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)

    for time_itr = 1:endtime
        for node_id1 = 1:(cmax-1)
            if isassigned(random_contacts_by_day, time_itr, node_id1) == false
                random_contacts_by_day[time_itr,node_id1] = Int64[]
            end
            for node_id2 = node_id1:cmax
                if isassigned(random_contacts_by_day, time_itr, node_id2) == false
                    random_contacts_by_day[time_itr,node_id2] = Int64[]
                end
                if rand(rng) < prob_random_contact
                    push!(random_contacts_by_day[time_itr,node_id1], node_id2)
                    push!(random_contacts_by_day[time_itr,node_id2], node_id1)
                end
            end
        end
    end

    return random_contacts_by_day::Array{Array{Int64,1},2}
end

"""
    generate_social_contacts_each_day_household_groups(args)

Generate daily, dynamic social contacts, from starttime to endtime, for both work and non-work days, according to specified number of households per group.

Each day, groups are made up of a specified number of households. Size of group is fixed as specified. Days can be repeated.
Used for rule of 6 analysis.

Inputs: `rng` - random number generator,
        `cmax` - number of nodes in network,
        `starttime` - time from which to generate daily social contacts,
        `endtime` - length of simulation,
        `n_contacts_per_group` - number of contacts in each group,
        `n_groups_per_day` - number of meeting groups per day,
        `dynamic_time_frame` - number of days to repeat contacts for,
        `n_households_per_group` - number of households per meeting group,
        `contacts` - ContactStructure structure \n
Outputs: `workday_social_contacts_by_day` - 2D array{array} of social contacts made per day per node (workday),
         `nonworkday_social_contacts_by_day` - 2D array{array} of social contacts made per day per node (non-workday) \n
Location: network_generation_fns.jl
"""
function generate_social_contacts_each_day_household_groups(rng::MersenneTwister,
                                            cmax::Int64,
                                            starttime::Int64,
                                            endtime::Int64,
                                            n_contacts_per_group::Int64,
                                            n_groups_per_day::Int64,
                                            dynamic_time_frame::Int64,
                                            n_households_per_group::Int64,
                                            contacts::ContactStructure)

    @unpack nodes_by_household = contacts

    if isassigned(contacts.social_meeting_group_sizes) == false
        contacts.social_meeting_group_sizes = Int64[]
    end

    workday_social_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)

    time_itr = starttime
    while time_itr <= endtime

        n_groups_per_node = ones(Int64, cmax)*n_groups_per_day
        for node_itr = 1:cmax
            workday_social_contacts_by_day[time_itr, node_itr] = Int64[]
        end

        nodes_remaining_per_household = [length(nodes_by_household[hh_itr]) for hh_itr = 1:length(nodes_by_household)]

        while sum(nodes_remaining_per_household) > 0

            households_remaining = findall(nodes_remaining_per_household.>0)

            current_hh_group = shuffle(rng, households_remaining)[1:min(n_households_per_group,length(households_remaining))]

            current_hh_group_expanded = vcat(fill.(current_hh_group, nodes_remaining_per_household[current_hh_group])...)

            current_group_hh_ids = shuffle(current_hh_group_expanded)[1:min(n_contacts_per_group,length(current_hh_group_expanded))]

            current_group = Int64[]
            for node_itr = 1:length(current_group_hh_ids)
                hh_idx = current_group_hh_ids[node_itr]
                hh_nodes_remaining = nodes_remaining_per_household[hh_idx]
                push!(current_group, nodes_by_household[hh_idx][hh_nodes_remaining])
                nodes_remaining_per_household[hh_idx] -= 1
            end

            group_size = length(current_group)
            push!(contacts.social_meeting_group_sizes, group_size)

            for ii = 1:(group_size-1)
                subject_id = current_group[ii]
                for jj = (ii+1):group_size
                    friend_id = current_group[jj]
                    if isassigned(workday_social_contacts_by_day,time_itr,subject_id)==false
                        workday_social_contacts_by_day[time_itr,subject_id] = Int64[]
                    end
                    if isassigned(workday_social_contacts_by_day,time_itr,friend_id)==false
                        workday_social_contacts_by_day[time_itr,friend_id] = Int64[]
                    end
                    push!(workday_social_contacts_by_day[time_itr,subject_id], friend_id)
                    push!(workday_social_contacts_by_day[time_itr,friend_id], subject_id)
                end
            end
        end

        if (dynamic_time_frame > 1) && (time_itr < endtime)
            days_left = endtime - time_itr
            if days_left >= dynamic_time_frame
                for time_rep = (time_itr+1):(time_itr+dynamic_time_frame-1)
                    workday_social_contacts_by_day[time_rep,:] = workday_social_contacts_by_day[time_itr,:]
                end
            else
                for time_rep = (time_itr+1):endtime
                    workday_social_contacts_by_day[time_rep,:] = workday_social_contacts_by_day[time_itr,:]
                end
            end
        end

        time_itr += dynamic_time_frame
    end

    nonworkday_social_contacts_by_day = workday_social_contacts_by_day

    return workday_social_contacts_by_day::Array{Array{Int64,1},2},
            nonworkday_social_contacts_by_day::Array{Array{Int64,1},2}
end


"""
    generate_social_contacts_each_day_household_groups_stochastic(args)

Generate daily, dynamic social contacts, from starttime to endtime, for both work and non-work days, according to specified number of households per group.

Each day, groups are made up of a specified number of households. Size of group is uniformly sampled between 0 and specified max. Days can be repeated.
Used for rule of 6 analysis.

Inputs: `rng` - random number generator,
        `cmax` - number of nodes in network,
        `starttime` - time from which to generate daily social contacts,
        `endtime` - length of simulation,
        `max_contacts_per_group` - maximum number of contacts in each group (uniformly sampled),
        `n_groups_per_day` - number of meeting groups per day,
        `dynamic_time_frame` - number of days to repeat contacts for,
        `n_households_per_group` - number of households per meeting group,
        `contacts` - ContactStructure structure,
        `intervention_set_itr` - number of current intervention set \n
Outputs: `workday_social_contacts_by_day` - 2D array{array} of social contacts made per day per node (workday),
         `nonworkday_social_contacts_by_day` - 2D array{array} of social contacts made per day per node (non-workday) \n
Location: network_generation_fns.jl
"""
function generate_social_contacts_each_day_household_groups_stochastic(rng::MersenneTwister,
                                            cmax::Int64,
                                            starttime::Int64,
                                            endtime::Int64,
                                            max_contacts_per_group::Int64,
                                            n_groups_per_day::Int64,
                                            dynamic_time_frame::Int64,
                                            n_households_per_group::Int64,
                                            contacts::ContactStructure,
                                            intervention_set_itr::Int64)

    @unpack nodes_by_household = contacts

    if isassigned(contacts.social_meeting_group_sizes, intervention_set_itr+1) == false
        contacts.social_meeting_group_sizes[intervention_set_itr+1] = Int64[]
    end

    workday_social_contacts_by_day = Array{Array{Int64,1},2}(undef,endtime,cmax)

    time_itr = starttime
    while time_itr <= endtime

        n_groups_per_node = ones(Int64, cmax)*n_groups_per_day
        for node_itr = 1:cmax
            workday_social_contacts_by_day[time_itr, node_itr] = Int64[]
        end

        nodes_remaining_per_household = [length(nodes_by_household[hh_itr]) for hh_itr = 1:length(nodes_by_household)]

        while sum(nodes_remaining_per_household) > 0

            households_remaining = findall(nodes_remaining_per_household.>0)

            current_hh_group = shuffle(rng, households_remaining)[1:min(n_households_per_group,length(households_remaining))]

            current_hh_group_expanded = vcat(fill.(current_hh_group, nodes_remaining_per_household[current_hh_group])...)

            n_contacts_per_group = rand(rng, DiscreteUniform(0,max_contacts_per_group))
            current_group_hh_ids = shuffle(current_hh_group_expanded)[1:min(n_contacts_per_group,length(current_hh_group_expanded))]

            current_group = Int64[]
            for node_itr = 1:length(current_group_hh_ids)
                hh_idx = current_group_hh_ids[node_itr]
                hh_nodes_remaining = nodes_remaining_per_household[hh_idx]
                push!(current_group, nodes_by_household[hh_idx][hh_nodes_remaining])
                nodes_remaining_per_household[hh_idx] -= 1
            end

            group_size = length(current_group)
            push!(contacts.social_meeting_group_sizes[intervention_set_itr+1], group_size)

            for ii = 1:(group_size-1)
                subject_id = current_group[ii]
                for jj = (ii+1):group_size
                    friend_id = current_group[jj]
                    if isassigned(workday_social_contacts_by_day,time_itr,subject_id)==false
                        workday_social_contacts_by_day[time_itr,subject_id] = Int64[]
                    end
                    if isassigned(workday_social_contacts_by_day,time_itr,friend_id)==false
                        workday_social_contacts_by_day[time_itr,friend_id] = Int64[]
                    end
                    push!(workday_social_contacts_by_day[time_itr,subject_id], friend_id)
                    push!(workday_social_contacts_by_day[time_itr,friend_id], subject_id)
                end
            end
        end

        if (dynamic_time_frame > 1) && (time_itr < endtime)
            days_left = endtime - time_itr
            if days_left >= dynamic_time_frame
                for time_rep = (time_itr+1):(time_itr+dynamic_time_frame-1)
                    workday_social_contacts_by_day[time_rep,:] = workday_social_contacts_by_day[time_itr,:]
                end
            else
                for time_rep = (time_itr+1):endtime
                    workday_social_contacts_by_day[time_rep,:] = workday_social_contacts_by_day[time_itr,:]
                end
            end
        end

        time_itr += dynamic_time_frame
    end

    nonworkday_social_contacts_by_day = workday_social_contacts_by_day

    return workday_social_contacts_by_day::Array{Array{Int64,1},2},
            nonworkday_social_contacts_by_day::Array{Array{Int64,1},2}
end
