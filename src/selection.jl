s_tournament(p::SingleObjPopulation, out_n::Int, poll_size::Int)::Vector{Int} = [findmin((p.fitness+p.penalty)[rand(1:p.size, poll_size)])[2] for i in 1:out_n]
s_tournament(p::MultiObjPopulation, out_n::Int, poll_size::Int)::Vector{Int} = [findmin(p.fitness[rand(1:p.size, poll_size)])[2] for i in 1:out_n]
s_identity(p::Population) = collect(1:p.size)
s_identity(p::MultiObjPopulation) = collect(1:p.size)
s_identity(p::MultiObjPopulation, n::Int) = append!(collect(1:p.size), collect(1:p.size))
s_identity(p::Population, n::Int) = append!(collect(1:p.size), collect(1:p.size))

function s_stochastic_tournament(p::Population, out_n::Int, poll_size::Int, constraints::Vector{Constraint}, P::Float64)::Vector{Int}
    rank = sortperm(stochastic_ranking(p, constraints, P))
    [findmin(rank[rand(1:p.size, poll_size)])[2] for i in 1:out_n]
end

function s_greedy_overselection(p::SingleObjPopulation, out_n::Int, prc::Float64)::Vector{Int}
    total_fitness = sum(p.fitness)
    elite_fitness = 0.0

    indexes = sortperm(p.fitness)
    elite = Vector{Int}()
    pop = Vector{Int}()

    i = 1
    while elite_fitness < total_fitness*prc
        push!(elite, indexes[i])
        elite_fitness += p.fitness[indexes[i]]
        i += 1
    end

    if i > length(indexes)
        return indexes
    end

    j = 1
    while length(pop) < out_n*0.8
        if j > length(elite)
            j = 1
        end
        push!(pop, indexes[j])
        j += 1
    end

    k = i
    println(k)

    while length(pop) < out_n
        if i > length(indexes)
            i = k
        end
        push!(pop, indexes[i])
        i += 1
    end

    return pop
end