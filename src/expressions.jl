
### EVALUATION

# TODO replace with mine
function getvars(ex::Expr)
    syms = Symbol[]
    for e in ex.args
        isa(e, Symbol) && push!(syms, e)
        isa(e, Expr) && append!(syms, getvars(e))
    end
    unique!(syms)
end

#= getvars(e::Expr)::Vector{Symbol} = Vector{Symbol}(unique(filter(x -> !isnothing(x), [getvars(s) for s in e.args])))
getvars(e::Symbol) = e
getvars(e::Any) = nothing =#

function evaluate(e::Expr, symbols::Dict{Symbol,Int}, values::T...)::T where T
	operation = e.args[1]
	arguments = [evaluate(subexpr, symbols, values...) for subexpr in e.args[2:end]]
    inf_idxs = findall(x->(x==Inf)||(x==-Inf)||(x==NaN), arguments)
    if length(inf_idxs) > 0

        # TODO replace the subtree returning nonsense by sense
        #consts = round.(rand(length(inf_idxs))*100, digits=3)
        consts = round.(rand(length(inf_idxs))*10 .-5, digits=3)
        #consts = zeros(length(inf_idxs))

        e.args[inf_idxs.+1] = consts
        arguments = [evaluate(subexpr, symbols, values...) for subexpr in e.args[2:end]]
    end
    return operation(arguments...)
end

function evaluate(e::Symbol, symbols::Dict{Symbol,Int}, values::T...)::T where T
    symbol_id = get(symbols, e, -1)
    return T(values[symbol_id])
end

function evaluate(e::Any, symbols::Dict{Symbol,Int}, values::T...)::T where T
    return e
end

### EDITING

get_subexp_cnt(e::Expr) = length(get_subexpressions(e))
get_subexpressions(e::Any) = nothing
function get_subexpressions(e::Expr)
    subexps = [e]
    for arg in e.args
        subexp = get_subexpressions(arg)
        if !isnothing(subexp)
            append!(subexps, subexp)
        end
    end
    return subexps
end

# by AI
#= function flatten(expr)
    result = []
    if isa(expr, Vector)
        for item in expr
            append!(result, flatten(item))
        end
    else
        push!(result, expr)
    end
    return result
end =#

flatten(expr) = isa(expr, Vector) ? vcat(flatten.(expr)...) : [expr]
# by AI

get_leaves(e::Expr)::Vector{Any} = Vector{Any}(filter(x -> !isnothing(x), [get_leaves(s) for s in e.args[2:end]]))	# args[1] is a function
get_leaves(e::Symbol) = e
get_leaves(e::T) where T<:Real = e
get_leaves(e::Any) = nothing

#= function replace_subexp_helper!(expr::Expr, index::Vector{Int}, replacement::Any, target_index::Int)
    index[1] += 1
    
	println("index: ", index[1])
	println("target_index: ", target_index)

	if index[1] > target_index
        return expr
	end
	if index[1] == target_index
        return replacement
    else
        for i in 1:length(expr.args)
            expr.args[i] = replace_subexp_helper!(expr.args[i], index, replacement, target_index)
        end
	end
	return expr
end

function replace_subexp_helper!(expr::Any, index::Vector{Int}, replacement::Any, target_index::Int)
#=     for i in 1:length(expr.args)
        expr.args[i] = replace_subexp_helper!(expr.args[i], index, replacement)
    end =#
	return expr
end =#
function replace_subexp!(e::Expr, target_index::Int, replacement::Any)
    index = [0]

    # Define a recursive helper function
    function replace!(expr::Expr, index::Vector{Int})
        index[1] += 1
        
        if index[1] == target_index
            return replacement  # Replace the specific occurrence
        else
            for i in 1:length(expr.args)
                if isa(expr.args[i], Expr)
                    expr.args[i] = replace!(expr.args[i], index)  # Recursive replacement
                end
				if index[1] >= target_index
					break
				end
            end
        end
        return expr
    end
    
    return replace!(e, index)
end

# Function to replace a specific subexpression instance
#= function replace_subexp!(e::Expr, target_index::Int, replacement::Any)
	#println("started subex")
    index = [0]
	r = replace_subexp_helper!(e, index, replacement, target_index)
	#println("finished subex")
	return r
end =#


function replace_leaf_helper!(expr::Expr, index::Vector{Int}, replacement::Any, target_index::Int)
	for i in 2:length(expr.args)
		expr.args[i] = replace_leaf_helper!(expr.args[i], index, replacement, target_index)
	end
	return expr
end

function replace_leaf_helper!(expr::Any, index::Vector{Int}, replacement::Any, target_index::Int)
	index[1] += 1
	if index[1] == target_index
		return replacement
	else
		return expr
	end
end

# Function to replace a specific subexpression instance
function replace_leaf!(e::Expr, target_index::Int, replacement::Any)
    index = [0]
	return replace_leaf_helper!(e, index, replacement, target_index)
	return r
end

### PRUNING

prune!(e::Any) = e
function prune!(e::Expr)
    e.args = prune!.(e.args)
    return apply_rule!(e)
end

function apply_rule!(e::Expr)

#=     # nullyfy
    if e.args[1] == :- && e.args[2] == e.args[3]
        return 0
    end
    if e.args[1] == :* && (e.args[3] == 0 || e.args[2] == 0)
        return 0
    end
    if e.args[1] == :/ && e.args[2] == 0
        return 0
    end

    # identity
    if e.args[1] == :* && e.args[2] == 1
        return e.args[3]
    end
    if e.args[1] == :* && e.args[3] == 1
        return e.args[2]
    end
    if (e.args[1] == :+ || e.args[1] == :- ) && e.args[3] == 0
        return e.args[2]
    end
    if (e.args[1] == :+ ) && e.args[2] == 0
        return e.args[3]
    end
    if e.args[1] == :/ && e.args[3] == 1
        return e.args[2]
    end

    # unify
    if e.args[1] == :/ && e.args[2] == e.args[3]
        return 1
    end    
    if e.args[1] == :cos && e.args[2] == 0
        return pi
    end
    if e.args[1] == :sin && e.args[2] == 0
        return 0
    end

    # undefined
    if e.args[1] == :/ && e.args[3] == 0
        return 0
    end =#

        # nullyfy
        if e.args[1] == (-) && e.args[2] == e.args[3]
            return 0
        end
        if e.args[1] == (*) && (e.args[3] == 0 || e.args[2] == 0)
            return 0
        end
        if e.args[1] == (/) && e.args[2] == 0
            return 0
        end
    
        # identity
        if e.args[1] == (*) && e.args[2] == 1
            return e.args[3]
        end
        if e.args[1] == (*) && e.args[3] == 1
            return e.args[2]
        end
        if (e.args[1] == (+) || e.args[1] == (-) ) && e.args[3] == 0
            return e.args[2]
        end
        if (e.args[1] == (+) ) && e.args[2] == 0
            return e.args[3]
        end
        if e.args[1] == (/) && e.args[3] == 1
            return e.args[2]
        end
    
        # unify
        if e.args[1] == (/) && e.args[2] == e.args[3]
            return 1
        end    
        if e.args[1] == (cos) && e.args[2] == 0
            return 1
        end
        if e.args[1] == (sin) && e.args[2] == 0
            return 0
        end
    
        # undefined
        if e.args[1] == (/) && e.args[3] == 0
            return 0
        end

    return e
end

function optimize_constants!(p::Population, constants::Vector{Vector{Float64}}, x::Vector{Vector{Float64}}, y::Vector{Float64})

    println("in")

    function f(x::RealChromosome, p::Population, y::Vector{Float64}, params::Vector{Float64}...)
        res = 0.0
        vars = append!(params..., x.sequence)
        for some in p.population
            res += sum((y .- Expr_parser(some.sequence).(vars...)).^2)
        end
        return res
    end

    initial_fit = f(get_real_chromosome([c[1] for c in constants]), p, y, x...)
    println(initial_fit)

    println("in")
    pop_size = 10
    dimension = length(constants)
    min = minimum(constants)[1]
    max = maximum(constants)[1]

    objective_function  = enclose_arguments(f, p, y, x...)
    initialization      = enclose_noargs(interval_real_initialization, dimension, pop_size, objective_function, min, max)
    selection           = enclose_arguments(s_tournament, pop_size, 5)
    crossover           = enclose_arguments(cr_parent_sum, pop_size)
    mutation            = enclose_arguments(gaussian_preturbation!, 10.)
    replacement         = enclose_replacement(r_best_n_diverse, pop_size, 0.7)
    termination         = enclose_argument(iteration_termination, pop_size*10)

    solution = solvink_hart(objective_function, initialization, selection, crossover, mutation, replacement, termination)

    println(solution.top_value)
    println([c[1] for c in constants])

    if solution.top_value >= initial_fit
        return
    end

    for i in 1:length(solution)
        constants[i] = solution.top_coords[i]
    end
end