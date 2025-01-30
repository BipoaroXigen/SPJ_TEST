r_replacement(old::Population, new::Population)::Population = new
r_merge(old::MultiObjPopulation, new::MultiObjPopulation)::MultiObjPopulation = MultiObjPopulation(append!(old.population, new.population), append!(old.fitness, new.fitness), append!(old.penalty, new.penalty), old.size+new.size)
r_merge(old::SingleObjPopulation, new::SingleObjPopulation)::SingleObjPopulation = SingleObjPopulation(append!(old.population, new.population), append!(old.fitness, new.fitness), append!(old.penalty, new.penalty), old.size+new.size)

function r_tournament(old::SingleObjPopulation, new::SingleObjPopulation, n::Int, poll_size::Int)::Population
    p = EO.r_merge(old, new)
    order = [findmin((p.fitness+p.penalty)[rand(1:p.size, poll_size)])[2] for i in 1:n]
    p.population = p.population[order][1:n]
    p.fitness = p.fitness[order][1:n]
    p.penalty = p.penalty[order][1:n]
    p.size = n
    return p
end

function r_best_n_diverse(old::Population, new::Population, n::Int, r::Float64)::Population
    @assert r > 0
    @assert r <= 1
    p = EO.r_merge(old, new)

    good_cnt = ceil(Int, r*n)
    bad_cnt = n-good_cnt

    i = 1
    for v in p.fitness
        if isnan(v)
            p.fitness[i] = Inf
        end
        i += 1
    end

    good_order = sortperm(p.fitness+p.penalty)[1:good_cnt]
    bad_order = rand(1:p.size, bad_cnt)

    order = vcat(good_order, bad_order)

    p.population = p.population[order]
    p.fitness = p.fitness[order]
    p.penalty = p.penalty[order]
    p.size = n
    return p
end

function r_keep_best_n(old::Population, new::Population, n::Int)::Population
    p = EO.r_merge(old, new)
    order = sortperm(p.fitness.+p.penalty)
    p.population = p.population[order][1:n]
    p.fitness = p.fitness[order][1:n]
    p.penalty = p.penalty[order][1:n]
    p.size = n
    return p
end

function r_keep_best_n_stoch(old::Population, new::Population, n::Int, constraints::Vector{Constraint}, P::Float64)::Population
    p = EO.r_merge(old, new)                  
    order = stochastic_ranking(p, constraints, P)
    p.population = p.population[order][1:n]
    p.fitness = p.fitness[order][1:n]
    p.penalty = p.penalty[order][1:n]
    p.size = n
    return p
end

function r_NSGA(old::T, new::T, n::Int)::T where {T<:Population}
    p = EO.r_merge(old, new)
    binary_tournament = [get_domination_count(p.fitness) -crowding_distance(p.fitness)] # minimize domination count and maximize crowding distance
    order = sortperm(eachslice(binary_tournament, dims=1))
    p.population = p.population[order][1:n]
    p.fitness = p.fitness[order][1:n]
    #p.penalty = p.penalty[order][1:n]      not needed for multiobjective
    p.size = n
    return p
end

"""constrained version of NSGA"""
function r_cNSGA(old::T, new::T, n::Int, constraints::Vector{Constraint})::T where {T<:Population}
    p = EO.r_merge(old, new)
    binary_tournament = [get_domination_count(p.fitness) G(p, constraints) -crowding_distance(p.fitness)] # minimize domination count and maximize crowding distance
    #binary_tournament = [get_domination_count(p.fitness) G(p, constraints)] # minimize domination count and maximize crowding distance
    order = sortperm(eachslice(binary_tournament, dims=1))
    p.population = p.population[order][1:n]
    p.fitness = p.fitness[order][1:n]
    #p.penalty = p.penalty[order][1:n]
    p.size = n
    return p
end

function r_ASCHEA(old::T, new::T, n::Int, tau::Float64, constraints::Vector{Constraint})::T where {T<:Population}
    p = EO.r_merge(old, new)
    count = p.size*tau/2    # the pop size doubled

    # Get feasible solutions sorted by fitness
    feasible_idx = findall(v -> all(x -> x == 0, v), G(p, constraints))
    feasible_fitness = p.fitness[feasible_idx]                              # the fitness is not penalized
    sorted_feasible = feasible_idx[sortperm(feasible_fitness)]
    # Select top count feasible solutions
    selected = sorted_feasible[1:min(length(sorted_feasible), ceil(Int, count))]

    # Get infeasible solutions sorted by fitness
    infeasible_idx = findall(v -> !all(x -> x == 0, v), G(p, constraints))
    infeasible_fitness = p.fitness[infeasible_idx] + p.penalty[infeasible_idx]  # penalized fitness here
    sorted_infeasible = infeasible_idx[sortperm(infeasible_fitness)]
    
    # Fill remaining slots with best infeasible solutions
    remaining = n - length(selected)
    remaining > 0 && append!(selected, sorted_infeasible[1:min(length(sorted_infeasible), remaining)])
    
    # Update population with selected solutions
    p.population = p.population[selected]
    p.fitness = p.fitness[selected]
    p.penalty = p.penalty[selected][1:n]
    p.size = length(selected)
    return p
end



function enclose_replacement(f::Function, a...)
    return (x...) -> f(x..., a...)
end