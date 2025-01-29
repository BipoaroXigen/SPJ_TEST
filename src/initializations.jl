function binary_initialization(dims::Int, size::Int, f::Function)<:Population
    population = Vector{BinaryChromosome}()
    for i in 1:size
        values = rand(0:1, dims)
        push!(population, get_binary_chromosome(values))
    end
    fitness = f.(population)
    return Population(population, fitness, size) 
end


function interval_real_initialization(dims::Int, size::Int, f::Function, lb::T, ub::T) where T<:Real
    function random_interval_vector(length, lb, ub)
        return rand(length)*(ub-lb) .+ lb
    end

    population = Vector{RealChromosome}()
    for i in 1:size
        push!(population, get_real_chromosome(random_interval_vector(dims, lb, ub)))
    end
    fitness = f.(population)
    return Population(population, fitness, size)
end

"""Dimensions, pop size"""
function TSP_initialization(dims::Int, size::Int, f::Function)<:Population
    population = Vector{RealChromosome}()
    for i in 1:size
        indexes = shuffle(collect(1:dims))
        push!(population, get_real_chromosome(indexes))
    end
    fitness = f.(population)
    return Population(population, fitness, size)
end

function TSP_NN_initialization(dims::Int, size::Int, f::Function, G::Matrix{Float64})<:Population
    population = Vector{RealChromosome}()
    for i in 1:size
        free = Bool.(ones(dims))
        free[i] = false
        indexes = [i]
        for _ in 1:dims-1
            closeness = sortperm(G[i, :])       # first is index of the closest vertex, last the farthest
            # pick the fisrt element of closeness for which free[element] == true
            nn_i = 0
            for j in 1:dims
                if free[closeness[j]] == true
                    nn_i = closeness[j]
                    break
                end
            end
            append!(indexes, nn_i)
            free[nn_i] = false
            i = nn_i
        end
        push!(population, get_real_chromosome(indexes))
    end
    fitness = f.(population)
    return Population(population, fitness, size)
end

function copy_initialization(dims::Int, size::Int, f::Function, seed::Vector{<:Real})<:Population
    population = Vector{RealChromosome}()
    for i in 1:size
        push!(population, get_real_chromosome(copy(seed)))
    end
    fitness = f.(population)
    return Population(population, fitness, size)
end

function one_expression_initialization(size::Int, f::Function, basis_functions::Vector{Function}, basis_variables::Vector{Any})#<:Population
    @assert size == 1
    population = Vector{ExprChromosome}()
    for i in 1:1
        d = 5
        root = prune!(random_expression_grow(basis_functions, basis_variables, d, 2))
        while !isa(root, Expr)
            root = prune!(random_expression_full(basis_functions, basis_variables, d))
        end
        push!(population, get_expr_chromosome(root))
    end
    fitness = f.(population)
    return Population(population, fitness, size) 
end

function expression_initialization(size::Int, f::Function, basis_functions::Vector{Function}, basis_variables::Vector{Any})#<:Population
    @assert size%2 == 0
    population = Vector{ExprChromosome}()
    for i in 1:size/2
        d = 5
        root = prune!(random_expression_full(basis_functions, basis_variables, d))
        while !isa(root, Expr)
            root = prune!(random_expression_full(basis_functions, basis_variables, d))
        end
        push!(population, get_expr_chromosome(root))
    end
    for i in 1:size/2
        d = 5
        root = prune!(random_expression_grow(basis_functions, basis_variables, d, 2))
        while !isa(root, Expr)
            root = prune!(random_expression_full(basis_functions, basis_variables, d))
        end
        push!(population, get_expr_chromosome(root))
    end
    fitness = f.(population)
    return Population(population, fitness, size) 
end

function random_expression_full(basis_functions::Vector{Function}, basis_variables::Vector{Any}, d::Int)
    if d == 0
        if rand() < 1/(length(basis_variables)+1)       # random initialization constants
            return round(rand()*10 -5, digits=3)
        end
        return rand(basis_variables)
    end

    operation = rand(basis_functions)
    if operation == sin || operation == cos || operation == square || operation == cube || operation == logaritmus
        a = random_expression_full(basis_functions, basis_variables, d-1)
        return :($(operation)($a))
    else
        a = random_expression_full(basis_functions, basis_variables, d-1)
        b = random_expression_full(basis_functions, basis_variables, d-1)
        return :($(operation)($a, $b))
    end
end

function random_expression_grow(basis_functions::Vector{Function}, basis_variables::Vector{Any}, d::Int, m::Int)
    if (rand() < 0.5 && m <= 0) || d == 0        
        if rand() < 1/(length(basis_variables)+1)       # random initialization constants
            return round(rand()*10 -5, digits=3)
        end
        return rand(basis_variables)
    end

    operation = rand(basis_functions)
    if operation == sin || operation == cos || operation == square || operation == cube || operation == logaritmus
        a = random_expression_grow(basis_functions, basis_variables, d-1, m-1)
        return :($(operation)($a))
    else
        a = random_expression_grow(basis_functions, basis_variables, d-1, m-1)
        b = random_expression_grow(basis_functions, basis_variables, d-1, m-1)
        return :($(operation)($a, $b))
    end
end

function square(x)
    return x^2
end
function cube(x)
    return x^3
end
logaritmus(x) = log(abs(x))