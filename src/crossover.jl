function cr_single_point(a::T, b::T)::T where T<:Chromosome
    if a.length == 2
        split = 1
    else
        split = rand(2:a.length-1)
    end
    return T(append!(a.sequence[1:split], b.sequence[split+1:end]))
end

function cr_single_point(p#= ::Population =#)#= ::Population =#
    @assert p.size %2 == 0

    p.population = [cr_single_point(p.population[i], p.population[i+1]) for i in 1:2:p.size]
    p.size = length(p.population)
    return p
end
function cr_single_point(p#= ::Population =#, out_n::Int)#= ::Population =#
    new_pop = Vector{RealChromosome}()
    for _ in 1:out_n
        a = rand(1:p.size)
        b = rand(1:p.size)
        push!(new_pop, cr_single_point(p.population[a], p.population[b]))
    end

    #p.population = [cr_single_point(p.population[i], p.population[i+1]) for i in 1:2:p.size]
    p.population = new_pop
    p.size = out_n
    return p
end

### TSP

function cr_ordered(a::T, b::T)::T where T<:Chromosome
    
    # select start and end
    idx = sort(shuffle(collect(1:a.length))[1:2])
	s = idx[1]
	e = idx[2]

    child = zeros(a.length)

    child[s:e] = copy(a.sequence[s:e])

    i = 1
    for ch in b.sequence
        if ch in child
            continue
        else
            while child[i] != 0
                i += 1
            end
            child[i] = ch
        end
    end

    return get_real_chromosome(Int.(child))
end

function cr_ordered(p::Population)::Population
    p.population = [cr_ordered(p.population[i], p.population[i+1]) for i in 1:p.size-1]
    return p
end

function cr_subtour(a::T, b::T)::Vector{T} where T<:Chromosome
    
    # select subtour in one parent
    idx = sort(shuffle(collect(1:a.length))[1:2])
	start_i = idx[1]
    end_i   = idx[2]

    c1 = copy(a.sequence)
    c2 = copy(b.sequence)
    subtour = copy(a.sequence[start_i:end_i])
    j = 0
    for i in 1:b.length
        if b.sequence[i] in subtour
            c1[start_i+j] = b.sequence[i]
            c2[i] = subtour[j+1]
            j += 1
        end
    end

#= 
    # select subtour in one parent
    idx = sort(shuffle(collect(1:a.length))[1:2])
    start_i = idx[1]
    end_i   = idx[2] =#


    return [get_real_chromosome(Int.(c1)),  get_real_chromosome(Int.(c2))]
end

"""1:1 size of input and output population"""
function cr_subtour(p::Population)::Population
    g = []
    for i in 1:p.size-1
        append!(g, cr_subtour(p.population[i], p.population[i+1]))
    end
    p.population = g
    return p
end

function cr_edge_recombination(a::T, b::T)::T where T<:Chromosome

    nbr_cnt(dict) = [length(dict[key]) for key in keys(dict)]
    function pick_nbr(v::Int, nbr_dict::Dict, used::Vector{Int})
        nbrs = shuffle(nbr_dict[v])
        for nbr in nbrs
            if used[nbr] == 0
                used[nbr] = 1
                return nbr
            end
        end

        for i in 1:length(used)
            if used[i] == 0
                used[i] = 1
                return i
            end
        end

        println("cr_edge_recombination broken")
        return idx[1]
    end
    function remove_nbr(nbr_dict::Dict, f::Int, s::Int)
        deleteat!(nbr_dict[f], findall(x->x==s, nbr_dict[f]))
        deleteat!(nbr_dict[s], findall(x->x==f, nbr_dict[s]))
    end

    neighbors = Int.(zeros(a.length, b.length)) #same lengts

    for i in 1:a.length

        if i > 1
            neighbors[a.sequence[i], a.sequence[i-1]] = 1
            neighbors[b.sequence[i], b.sequence[i-1]] = 1
            neighbors[a.sequence[i-1], a.sequence[i]] = 1
            neighbors[b.sequence[i-1], b.sequence[i]] = 1
        end

        if i < a.length
            neighbors[a.sequence[i], a.sequence[i+1]] = 1
            neighbors[b.sequence[i], b.sequence[i+1]] = 1
            neighbors[a.sequence[i+1], a.sequence[i]] = 1
            neighbors[b.sequence[i+1], b.sequence[i]] = 1
        end

        if i == 1
            neighbors[a.sequence[i], a.sequence[end]] = 1
            neighbors[b.sequence[i], b.sequence[end]] = 1
            neighbors[a.sequence[end], a.sequence[i]] = 1
            neighbors[b.sequence[end], b.sequence[i]] = 1
        end

        if i == a.length
            neighbors[a.sequence[i], a.sequence[1]] = 1
            neighbors[b.sequence[i], b.sequence[1]] = 1
            neighbors[a.sequence[1], a.sequence[i]] = 1
            neighbors[b.sequence[1], b.sequence[i]] = 1
        end

    end

    nbr_dict = Dict( i => Vector{Int}() for i in 1:a.length)

    for i in 1:a.length
        for j in 1:b.length
            if neighbors[i, j] == 1
                push!(nbr_dict[i], j)
            end
        end
    end

    used            = zeros(Int, a.length)
    child           = zeros(Int, a.length)      # the child sequence
    current         = rand(1:a.length)          # initial vertex
    used[current]   = 1
    child[1]        = current

    # pick a neigbor of vertex and move there
    for i in 2:a.length
        nbr = pick_nbr(current, nbr_dict, used)     # pick neighbor
        remove_nbr(nbr_dict, nbr, current)          # remove neighbor from list
        current = nbr                               # move to the neighbor
        child[i] = current
    end

    return get_real_chromosome(Int.(child))
end

function cr_edge_recombination(p::Population)::Population
    p.population = [cr_edge_recombination(p.population[i], p.population[i+1]) for i in 1:p.size-1]
    return p
end

function cr_parent_sum(a::T, b::T)::T where T<:Chromosome
    alphas = rand(a.length)*1.3 .- 0.3
    child = alphas.*a.sequence .+ (1 .-alphas).*b.sequence
    return get_real_chromosome(child)
end

function cr_parent_sum(p::Population, out_n::Int)::Population
    new_pop = Vector{RealChromosome}()
    for _ in 1:out_n
        a = rand(1:p.size)
        b = rand(1:p.size)
        push!(new_pop, cr_parent_sum(p.population[a], p.population[b]))
    end

    #p.population = [cr_single_point(p.population[i], p.population[i+1]) for i in 1:2:p.size]
    p.population = new_pop
    p.size = out_n
    return p
end

function cr_ASCHEA(a::T, b::T)::T where T<:Chromosome
    alphas = rand(a.length)*1.3 .- 0.3
    child = alphas.*a.sequence .+ (1 .-alphas).*b.sequence
    return get_real_chromosome(child)
end

function cr_ASCHEA(p::Population, out_n::Int, T::Float64, constraints::Vector{Constraint})::Population

    violations = G(p, constraints)
    good = count(v -> all(x -> x == 0, v), violations)

    # If we have too few feasible solutions, breed them with best infeasible ones
    if good/p.size < T && good > 0
        feasible_idx = findall(v -> all(x -> x == 0, v), violations)
        infeasible_idx = findall(v -> !all(x -> x == 0, v), violations)

        new_pop = Vector{RealChromosome}()
        
        # Sort infeasible solutions by fitness
        infeasible_fitness = (p.fitness+p.penalty)[infeasible_idx]
        sorted_infeasible  = infeasible_idx[sortperm(infeasible_fitness)]
        feasible_fitness   = p.fitness[feasible_idx]
        sorted_feasible    = feasible_idx[sortperm(feasible_fitness)]
        
        # Breed each feasible solution with top infeasible ones
        i = 1
        j = 1
        while length(new_pop) < out_n
            if i > length(sorted_infeasible)
                i = 1
            end
            if j > length(sorted_feasible)
                j = 1
            end
            a = sorted_infeasible[i]
            b = sorted_feasible[j]
            push!(new_pop, cr_ASCHEA(p.population[a], p.population[b]))
            i += 1
            j += 1
        end
        
        p.population = new_pop
        p.size = out_n
        return p
    end


    # else breed normally
    new_pop = Vector{RealChromosome}()
    for _ in 1:out_n
        a = rand(1:p.size)
        b = rand(1:p.size)
        push!(new_pop, cr_ASCHEA(p.population[a], p.population[b]))
    end

    #p.population = [cr_single_point(p.population[i], p.population[i+1]) for i in 1:2:p.size]
    p.population = new_pop
    p.size = out_n
    return p
end

function cr_subtree(a::ExprChromosome, b::ExprChromosome, basis_functions::Vector{Function}, basis_variables::Vector{Any})::Tuple{EO.ExprChromosome, EO.ExprChromosome}

    # copy parents
    new_a = copy(a)
    new_b = copy(b)

    if rand() < 0.9
        subexps_a = get_subexpressions(a.sequence)
        subexps_b = get_subexpressions(b.sequence)
        subexp_cnt_a = length(subexps_a)
        subexp_cnt_b = length(subexps_b)
        if subexp_cnt_a > 1 && subexp_cnt_b > 1
            root_a_i = rand(1:subexp_cnt_a)
            root_b_i = rand(1:subexp_cnt_b)

            root_a = copy(subexps_a[root_a_i])
            root_b = copy(subexps_b[root_b_i])

            new_a.sequence = replace_subexp!(new_b.sequence, root_b_i, root_a)
            new_b.sequence = replace_subexp!(new_a.sequence, root_a_i, root_b)
        else
            # swap only leaves
            leaves_a = flatten(get_leaves(new_a.sequence))
            leaves_b = flatten(get_leaves(new_b.sequence))

            rand_leaf_i_a = rand(1:length(leaves_a))
            rand_leaf_i_b = rand(1:length(leaves_b))

            rand_leaf_a = leaves_a[rand_leaf_i_a]
            rand_leaf_b = leaves_b[rand_leaf_i_b]

            new_a.sequence = replace_leaf!(new_a.sequence, rand_leaf_i_a, rand_leaf_b)
            new_b.sequence = replace_leaf!(new_b.sequence, rand_leaf_i_b, rand_leaf_a)
        end
    else
        # swap only leaves
		leaves_a = flatten(get_leaves(new_a.sequence))
		leaves_b = flatten(get_leaves(new_b.sequence))

		rand_leaf_i_a = rand(1:length(leaves_a))
		rand_leaf_i_b = rand(1:length(leaves_b))

		rand_leaf_a = leaves_a[rand_leaf_i_a]
		rand_leaf_b = leaves_b[rand_leaf_i_b]

		new_a.sequence = replace_leaf!(new_a.sequence, rand_leaf_i_a, rand_leaf_b)
		new_b.sequence = replace_leaf!(new_b.sequence, rand_leaf_i_b, rand_leaf_a)
    end

    return new_a, new_b
end

function cr_subtree(p::Population, out_n::Int, basis_functions::Vector{Function}, basis_variables::Vector{Any})::Population
    @assert out_n%2 == 0
    new_pop = Vector{ExprChromosome}()
    for _ in 1:out_n/2
        a = rand(1:p.size)
        b = rand(1:p.size)
        c1, c2 = cr_subtree(p.population[a], p.population[b], basis_functions, basis_variables)
        push!(new_pop, c1)
        push!(new_pop, c2)
    end

    p.population = new_pop
    p.size = out_n
    return p
end

function cr_GSGP(a::ExprChromosome, b::ExprChromosome, basis_functions::Vector{Function}, basis_variables::Vector{Any})::EO.ExprChromosome

    # copy parents
    new_a = copy(a)
    new_b = copy(b)

    #new_c = :(mean($(new_a.sequence), $(new_b.sequence)))
#=     new_c = :(($+)($(new_a.sequence), $(new_b.sequence)))
    new_c = :(($*)(0.5, $(new_c))) =#

    new_c = :(($mmean)($(new_a.sequence), $(new_b.sequence)))

    return get_expr_chromosome(new_c)
end

function cr_GSGP(p::Population, out_n::Int, basis_functions::Vector{Function}, basis_variables::Vector{Any})::Population
    new_pop = Vector{ExprChromosome}()
    for i in 1:out_n
        a = rand(1:p.size)
        b = rand(1:p.size)
        c = cr_GSGP(p.population[a], p.population[b], basis_functions, basis_variables)
        push!(new_pop, c)
    end

    p.population = new_pop
    p.size = out_n
    return p
end

mmean(x, y) = 0.5*(x+y)