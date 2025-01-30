function solvink_hart(objective_function::Function, initialization::Function, selection::Function, crossover::Function, mutation::Function, replacement::Function, termination::Function; 
    constraints::Union{Vector{Constraint}, NConstraint}=NConstraint(), penalty::Union{Function, NFunction}=NFunction())

    population = initialization()                                       # initialize population (evaluated)

    top_idx = findmin(population.fitness)
    top_f::Float64 = top_idx[1]
    top_x = population.population[top_idx[2]]

    f_diff = 7  # not used for now
    time = 0    # not used for nowrg = 100.0
    rg = 100.0
    beta1 = 1.1
    beta2 = 1.2
    k = 10      # how many times the optimal solution has to be feasible/unfeasible to update rgfeasibles = round(Int, k/2)
    feasibles = round(Int, k/2)
    unfeasibles = k .- feasibles

    d = sqrt(population.population[1].length)

    #x_history = Vector{Vector{typeof(population.population[1].sequence.args[1])}}()
    x_history = Vector{Any}()
    f_history = Vector{Float64}()
    p_history = Vector{Float64}()

    pophist = Vector{typeof(population)}()

    iteration = 0
    while !termination(iteration, f_diff, time)
        #println("generation: $(iteration/population.size); fitness: $(top_f); diversity: $(std(population.fitness))")
        prune!.(population.population)

        # select pant
        parents = copy.(population.population[selection(population)])   # select parent chromosomes
        parents = Population(parents)

        # create new generation
        children = crossover(parents)

        # mutate children
        mutation.(children.population)
        prune!.(children.population)
        children.fitness = objective_function.(children.population)
        children.penalty = penalty(children, rg, constraints)

        population = replacement(population, children)                  # create new population from both generations

        # update output
        improving_pop_size = sum(population.fitness + population.penalty .< top_f + penalty(top_x, rg, constraints))
        improving_r = improving_pop_size/population.size

        top_idx = findmin(population.fitness + population.penalty)
        if top_idx[1] < top_f + penalty(top_x, rg, constraints) || penalty isa NFunction
            top_x = population.population[top_idx[2]]
            top_f = top_idx[1]
            #top_f = objective_function(top_x)
        end

        population.fitness = objective_function.(population.population)
        population.penalty = penalty(population, rg, constraints)

        push!(f_history, top_f)
        if iteration%100 == 0
            push!(x_history, copy(top_x.sequence))
            #push!(f_history, top_f)
            push!(pophist, population)
            push!(p_history, penalty(top_x, rg, constraints))
        end

        if sum(constraints(top_x)) == 0         # does satisfy all constraints
            if unfeasibles > 0
                unfeasibles -= 1
                feasibles += 1
            end
        else
            if feasibles > 0
                feasibles -= 1
                unfeasibles += 1
            end
        end

        if feasibles/k <= 0.6
            rg = 1.1*rg
        else
            rg = rg/1.1
        end

#=         if unfeasibles == k
            rg = beta2 * rg
        elseif feasibles == k
            rg = rg / beta1
        end =#

        #= p = mutation.a[1] * (exp(improving_r-0.2))^(1/d)
        mutation = enclose_arguments(mutation.f, p) =#

        iteration += 1*population.size
    end
    #println("was?")
    #return Result_benchmark(f_history, p_history)
    #return f_history
    return Result_big(top_x.sequence, top_f, x_history, f_history, pophist, p_history)
    #return top_x.sequence
end

function solvink_hart(objective_function::Function, initialization::Function, selection::Function, crossover::Function, mutation::Function, replacement::Function, termination::Function, ::String; 
    constraints::Union{Vector{Constraint}, NConstraint}=NConstraint(), penalty::Union{Function, NFunction}=NFunction())

    population = initialization()                                       # initialize population (evaluated)

    top_idx = findmin(population.fitness)
    top_f::Float64 = top_idx[1]
    top_x = population.population[top_idx[2]]

    f_diff = 7  # not used for now
    time = 0    # not used for now
    rg = ones(length(constraints))*100.0*10000
    #rg = 100.0
    beta1 = 1.1
    beta2 = 1.2
    k = 10      # how many times the optimal solution has to be feasible/unfeasible to update rg
    feasibles = ones(length(constraints))*round(Int, k/2)
    #feasibles = round(Int, k/2)
    unfeasibles = k .- feasibles

    d = sqrt(population.population[1].length)

    x_history = Vector{Vector{typeof(population.population[1].sequence[1])}}()
    f_history = Vector{Float64}()
    p_history = Vector{Float64}()

    pophist = Vector{typeof(population)}()

    iteration = 0
    while !termination(iteration, f_diff, time)
        # select pant
        parents = copy.(population.population[selection(population)])   # select parent chromosomes
        parents = Population(parents)

        # create new generation
        children = crossover(parents)

        # mutate children
        mutation.(children.population)
        children.fitness = objective_function.(children.population)
        children.penalty = penalty(children, rg, constraints)

        population = replacement(population, children)                  # create new population from both generations

        # update output
        improving_pop_size = sum(population.fitness + population.penalty .< top_f + penalty(top_x, rg, constraints))
        improving_r = improving_pop_size/population.size

        top_idx = findmin(population.fitness + population.penalty)
        #top_f, top_x = get_BSF(top_idx, penalty, top_x, rg, constraints)
        if top_idx[1] < top_f + penalty(top_x, rg, constraints) || penalty isa NFunction
            #top_f = top_idx[1]
            #top_f = objective_function(population.population[top_idx[2]])
            top_x = population.population[top_idx[2]]
            top_f = objective_function(top_x)
        end

        population.fitness = objective_function.(population.population)
        population.penalty = penalty(population, rg, constraints)

        if iteration%10*population.size== 0
            push!(x_history, copy(top_x.sequence))
            push!(f_history, top_f)
            push!(pophist, population)
            push!(p_history, penalty(top_x, rg, constraints))
        end

        violations = constraints(top_x)
        
        for i in 1:length(violations)
            if unfeasibles[i] > 0 && violations[i] > 0
                unfeasibles[i] -= 1
                feasibles[i] += 1
            elseif feasibles[i] > 0 && violations[i] <= 0
                unfeasibles[i] += 1
                feasibles[i] -= 1
            end
            if feasibles[i]/k <= 0.6
                rg[i] = 1.1*rg[i]
            else
                rg[i] = rg[i]/1.1
            end
        end

        p = mutation.a[1] * (exp(improving_r-0.2))^(1/d)
        mutation = enclose_arguments(mutation.f, p)

        iteration += 1*population.size
    end

    println(top_f)
    #return Result_benchmark(f_history, p_history)
    return Result_big(top_x.sequence, top_f, x_history, f_history, pophist, p_history)
end

"""multiobjective"""
function solvink_hart(objective_function::MultiObjFunction, initialization::Function, selection::Function, crossover::Function, mutation::Function, replacement::Function, termination::Function)
    population = initialization()                                       # initialize population (evaluated)

    pareto_front = findall(get_domination_count(population.fitness).==0)# indexes of non-dominated solutions

    f_diff = 7  # not used for now
    time = 0    # not used for now

    d = sqrt(population.population[1].length)

    #x_history = Vector{Vector{typeof(population.population[1].sequence[1])}}()
    x_history = Vector{Any}()
    pophist   = Vector{typeof(population)}()
    f_history = Vector{Float64}()
    p_history = Vector{Float64}()
    #= f_history = Vector{Vector{Float64}}() =#

    iteration = 0
    while !termination(iteration, f_diff, time)

        prune!.(population.population)
        i = 1
        for v in population.fitness
            if isnan(v[1])
                population.fitness[i][1] = Inf
            end
            i += 1
        end
        #println(population.fitness[1:5])

        # select pant
        parents = copy.(population.population[selection(population)])       # select parent chromosomes
        parents = MultiObjPopulation(parents)

        # create new generation
        children = crossover(parents)

        # mutate children
        mutation.(children.population)
        prune!.(children.population)
        children.fitness = objective_function.(children.population)
        i = 1
        for v in children.fitness
            if isnan(v[1])
                children.fitness[i][1] = Inf
            end
            i += 1
        end


        population = replacement(population, children)                      # create new population from both generations


        pareto_front = findall(get_domination_count(population.fitness).==0)# indexes of non-dominated solutions

#=
        current_f, top_idx = findmin(population.fitness)
        if current_f < top_f
            top_f = current_f
            top_x = population.population[top_idx]
        end =#

        if iteration%10*population.size== 0

            violations = [sum(f[2:end]) for f in population.fitness]
            
            # Find solution with minimum sum of other objectives
            min_idx = argmin(violations)
            min_sum = violations[min_idx]
            min_solution = population.population[min_idx]
            
            push!(x_history, copy(min_solution.sequence))
            push!(p_history, min_sum)  # Save sum of other objectives
            push!(f_history, population.fitness[min_idx][1])  # Save first objective
            #push!(x_history, copy(top_x.sequence))
            #push!(f_history, top_f)
            push!(pophist, population)
        end

        improving_pop_size = length(pareto_front)
        improving_r = improving_pop_size/population.size

        #p = mutation.a[1] * (exp(improving_r-0.2))^(1/d)
        #mutation = enclose_arguments(mutation.f, p)

        iteration += 1*population.size
    end

    #return Result_real(top_x.sequence, top_f, x_history, f_history)
    return Result_big(x_history[end], f_history[end], x_history, f_history, pophist, p_history)
    #return Result_multi(x_history[end], vcat(f_history[end], p_history[end]), x_history, f_history, pophist, p_history)
    #return Result_multi(x_history, pophist, f_history, p_history)
    #return Result_benchmark(f_history, p_history)
end

function solvink_hart(objective_function::Function, initialization::Function, selection::Function, crossover::Function, mutation::Function, replacement::Function, termination::Function,
    constants::Vector{Vector{Float64}}, variables::Vector{Vector{Float64}}, y::Vector{Float64}; 
    constraints::Union{Vector{Constraint}, NConstraint}=NConstraint(), penalty::Union{Function, NFunction}=NFunction())

    population = initialization()                                       # initialize population (evaluated)

    top_idx = findmin(population.fitness)
    top_f::Float64 = top_idx[1]
    top_x = population.population[top_idx[2]]

    f_diff = 7  # not used for now
    time = 0    # not used for nowrg = 100.0
    rg = 100.

    #x_history = Vector{Vector{typeof(population.population[1].sequence.args[1])}}()
    x_history = Vector{Any}()
    f_history = Vector{Float64}()
    p_history = Vector{Float64}()

    pophist = Vector{typeof(population)}()

    iteration = 0
    while !termination(iteration, f_diff, time)
#=         println("generation: $(iteration/population.size); fitness: $(top_f); diversity: $(std(population.fitness))")
        prune!.(population.population)

        println(population.fitness)
        println(population.penalty) =#

        # select pant
        parents = copy.(population.population[selection(population)])   # select parent chromosomes
        parents = Population(parents)

        # create new generation
        children = crossover(parents)

        # mutate children
        mutation.(children.population)
        prune!.(children.population)
        children.fitness = objective_function.(children.population)
        children.penalty = penalty(children, rg, constraints)

        population = replacement(population, children)                  # create new population from both generations

        # update output
        improving_pop_size = sum(population.fitness + population.penalty .< top_f + penalty(top_x, rg, constraints))
        improving_r = improving_pop_size/population.size

        top_idx = findmin(population.fitness + population.penalty)
        if top_idx[1] < top_f + penalty(top_x, rg, constraints) || penalty isa NFunction
            top_x = population.population[top_idx[2]]
            top_f = top_idx[1]
        end

        population.fitness = objective_function.(population.population)
        population.penalty = penalty(population, rg, constraints)

        if iteration%100 == 0
            push!(x_history, copy(top_x.sequence))
            push!(f_history, top_f)
            push!(pophist, population)
            push!(p_history, penalty(top_x, rg, constraints))
        end

        optimize_constants!(Population(copy.(population.population[selection(population)])), constants, variables, y)

        iteration += 1*population.size
    end
    return Result_big(top_x.sequence, top_f, x_history, f_history, pophist, p_history)
end