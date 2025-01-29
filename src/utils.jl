function bool_to_units(x::Vector{Bool})::Vector{Int}
	ret = ones(Int, length(x))
	mask = x .== false
	ret[mask] .= -1
	return ret
end

function sequence_correlation(x::Vector{Int}, k::Int)::Float64
	c = 0
	D = length(x)
	for i in 1:D-k
		c += x[i]*x[i+k]
	end
	return c
end

function bin2real(lb::Real, ub::Real, x::BinaryChromosome)::Vector{Float64}
	L = parse(Int, join(ones(Int, x.length)), base=2)+1
	range = LinRange(lb, ub, L)
	index = parse(Int, join(Int.(x.sequence)), base=2)

	[range[index+1]]
end

function bin2real(lb::Vector{<:Real}, ub::Vector{<:Real}, x::BinaryChromosome)::Vector{Float64}
    if x.length % length(lb) != 0
        display("The binary vector length is not divisible by the dimensionality of the target vector space.")
        return []
    end

    step = Int(x.length / length(lb))
    result = Vector{Float64}()

	L = parse(Int, join(ones(Int, Int(x.length/length(lb)))), base=2)+1
    for i in 1:length(lb)
	    range = LinRange(lb[i], ub[i], L)
	    index = parse(Int, join(Int.(x.sequence[(i-1)*step+1:(i)*step])), base=2)
        append!(result, range[index+1])
    end
    result
end

"""
Returns indices of non-dominated solutions from a set of solutions and their objective scores.
A solution is non-dominated if no other solution is better in all objectives.

Parameters:
- solutions: Vector of solutions (any type)
- scores: Vector of objective score vectors for each solution
Returns vector of indices of non-dominated solutions

WRITTEN BY AI. It wrote functional code on first try
"""
function get_domination_count(scores::Vector{Vector{Float64}})
    n = length(scores)
    domination_count = zeros(Int, n)
    
    for i in 1:n
        for j in 1:n
            if i != j
                # Check if solution j dominates solution i
                if all(scores[j] .<= scores[i]) && any(scores[j] .< scores[i])
                    domination_count[i] += 1
                end
            end
        end
    end
    
	return domination_count
end

"""
Computes crowding distances for a set of solutions based on their multi-objective fitness values.
Crowding distance measures how far a solution is from its neighbors in objective space.
Used in NSGA-II for maintaining diversity.

Parameters:
- fitness: Vector of fitness vectors, where each inner vector contains multiple objective values
Returns vector of crowding distances

Reference: Deb, K., et al. "A fast and elitist multiobjective genetic algorithm: NSGA-II." 
IEEE transactions on evolutionary computation 6.2 (2002): 182-197.

Again written by AI, I was curious how well it works. It messed up indexing, but that aside it wrote a working function...
"""
function crowding_distance(fitness::Vector{Vector{Float64}})::Vector{Float64}
    n = length(fitness)  # Number of solutions
    if n < 3
        return fill(Inf, n) # Edge case - assign infinite distance if too few solutions
    end
    
    #m = length(fitness[1]) # Number of objectives
    m = length(fitness[1])
    distances = zeros(n)
    
    for obj_i in 1:m
        # Extract values for current objective
        obj_values = [f[obj_i] for f in fitness]  #it came up with very unintuitive interpretation of fitness vector
        #obj_values = fitness[obj_i]
        
        # Sort indices by current objective
        sorted_indices = sortperm(obj_values)
        
        # Set boundary points to infinity
        distances[sorted_indices[1]] = Inf
        distances[sorted_indices[end]] = Inf
        
        # Get normalization factor for this objective
        f_min = obj_values[sorted_indices[1]]
        f_max = obj_values[sorted_indices[end]]
        norm = max(f_max - f_min, 1e-10) # Avoid division by zero
        
        # Calculate distances for intermediate points
        for i in 2:n-1
            curr_idx = sorted_indices[i]
            if distances[curr_idx] != Inf # Only update if not already a boundary point # unnecessary condition
                prev_idx = sorted_indices[i-1]
                next_idx = sorted_indices[i+1]
                # Add normalized distance between neighbors for this objective
                distances[curr_idx] += abs(obj_values[next_idx] - obj_values[prev_idx]) / norm  # it forgot abs here
            end
        end
    end
    
    return distances
end

"""supply a path to the german tsp thing, without extensions ;)"""
function parse_TSP_problem(file_path::String)#::Matrix{Float64}, Vector{Vector{Float64}}
	problem_path = file_path*".tsp"

	vertices = Vector{Vector{Float64}}()

	edging = false
	open(problem_path) do f
		while !eof(f)
			line = readline(f)
			if occursin("DIMENSION", line)
				d = parse(Int64, split(line)[end])
				global G = zeros(d, d)
			end
			if occursin("NODE_COORD_SECTION", line)
				edging = true
				continue
			end
			if occursin("EOF", line)
				break
			end
			if edging == true
				vertex = split(line)[2:end]
				vertex = parse.(Float64, vertex)
				push!(vertices, vertex)
			end
		end
	end

	# get all distances
	for i in 1:length(vertices)
		for j in 1:length(vertices)
			G[i, j] = sqrt(sum((vertices[i]-vertices[j]).^2))
			G[j, i] = G[i, j]
		end
	end
	
	return G, vertices
end

function benchmark(objective_function::Function, initialization::Function, selection::Function, crossover::Function, mutation::Function, replacement::Function, termination::Function, K::Int, prefix::String)
	for k in 1:K
		run_output = solvink_hart(objective_function, initialization, selection, crossover, mutation, replacement, termination)
		open("benchmarks/"*prefix*"_"*string(k), "w") do f
			println(f, "COORDS")
			for coord in run_output.coords_history
				println(f, string(coord)[2:end-1])
			end
			println(f, "VALUES")
			for value in run_output.value_history
				println(f, string(value))
			end
		end
	end
end

function setup_simple_problems()
	simple_problems = Vector{Problem}()
	f_1(x) = (x.sequence[1]-10)^3 + (x.sequence[2]-20)^3
	constraints = [EO.constraint(x -> -(x[1]-5)^2 - (x[2]-5)^2 + 100, .<=, 0.), EO.constraint(x -> (x[1]-6)^2 + (x[2]-5)^2 − 82.81, .<=, 0.),
					EO.constraint(x -> -x[1]+13, <=, 0.), EO.constraint(x -> x[1]-100, <=, 0.), 
					EO.constraint(x -> -x[2], <=, 0.), EO.constraint(x -> x[2]-100, <=, 0.)]
	push!(simple_problems, Problem("g06", f_1, constraints, 2, [14.09500000000000064, 0.8429607892154795668], −6961.8138755801))
	f_2(x) = -(sin(2*pi*x.sequence[1])^3 * sin(2*pi*x.sequence[2]))/(x.sequence[1]^3*(x.sequence[1] + x.sequence[2]))
	constraints = [EO.constraint(x -> x[1]^2 - x[2] + 1, <=, 0.), EO.constraint(x -> 1 - x[1] + (x[2]-4)^2, <=, 0.),
					EO.constraint(x -> -x[1], <=, 0.), EO.constraint(x -> x[1] - 10, <=, 0.), 
					EO.constraint(x -> -x[2], <=, 0.), EO.constraint(x -> x[2] - 10, <=, 0.)]
	push!(simple_problems, Problem("g08", f_2, constraints, 2, [1.22797135260752599, 4.24537336612274885], −0.0958250414180359))
	f_3(x) = x.sequence[1]^2 + (x.sequence[2]-1)^2
	constraints = [EO.constraint(x -> x[2] - x[1]^2 - 1e-7, <=, 0.), EO.constraint(x -> -x[2] + x[1]^2 + 1e-7, <=, 0.),
					EO.constraint(x -> -x[1]-1, <=, 0.), EO.constraint(x -> x[1]-1, <=, 0.), 
					EO.constraint(x -> -x[2]-1, <=, 0.), EO.constraint(x -> x[2]-1, <=, 0.)]
	push!(simple_problems, Problem("g11", f_3, constraints, 2, [-0.707036070037170616, 0.500000004333606807], 0.7499))
	f_4(x) = - x.sequence[1] - x.sequence[2]
	constraints = [EO.constraint(x -> -2*x[1]^4 + 8*x[1]^3 - 8*x[1]^2 + x[2] - 2, <=, 0.), EO.constraint(x -> -4*x[1]^4 + 32*x[1]^3 - 88*x[1]^2 + 96*x[1] + x[2] - 36, <=, 0.),
					EO.constraint(x -> -x[1], <=, 0.), EO.constraint(x -> x[1]-3, <=, 0.), 
					EO.constraint(x -> -x[2], <=, 0.), EO.constraint(x -> x[2]-4, <=, 0.)]
	push!(simple_problems, Problem("g24", f_4, constraints, 2, [2.329520197477623, 3.17849307411774], −5.5080132715953))
	return simple_problems
end

function setup_unsimple_problems()
	simple_problems = Vector{Problem}()
	f_1(x) = 5.3578547*x.sequence[3]^2 + x.sequence[1]*x.sequence[5]*0.8356891 + x.sequence[1]*37.293239 − 40792.141
	constraints = [EO.constraint(x -> 85.334407 + 0.0056858*x[2]*x[5] + 0.0006262*x[1]*x[4] - 0.0022053*x[3]*x[5] - 92, <=, 0.)]
	push!(constraints, EO.constraint(x -> -85.334407 - 0.0056858*x[2]*x[5] - 0.0006262*x[1]*x[4] + 0.0022053*x[3]*x[5], <=, 0.))
	push!(constraints, EO.constraint(x -> 80.51249 + 0.0071317*x[2]*x[5] + 0.0029955*x[1]*x[2] + 0.0021813*x[3]^2 - 110, <=, 0.))
	push!(constraints, EO.constraint(x -> -80.51249 - 0.0071317*x[2]*x[5] - 0.0029955*x[1]*x[2] - 0.0021813*x[3]^2 + 90, <=, 0.))
	push!(constraints, EO.constraint(x -> 9.300961 + 0.0047026*x[3]*x[5] + 0.0012547*x[1]*x[3] + 0.0019085*x[3]*x[4] - 25, <=, 0.))
	push!(constraints, EO.constraint(x -> -9.300961 - 0.0047026*x[3]*x[5] - 0.0012547*x[1]*x[3] - 0.0019085*x[3]*x[4] + 20, <=, 0.))
	push!(constraints, EO.constraint(x -> -x[1] + 78, <=, 0.))
	push!(constraints, EO.constraint(x -> x[1] - 102, <=, 0.))
	push!(constraints, EO.constraint(x -> -x[2] + 33, <=, 0.))
	push!(constraints, EO.constraint(x -> x[2] - 45, <=, 0.))
	for i in 3:5
		push!(constraints, EO.constraint(x -> -x[i] + 27, <=, 0.))
		push!(constraints, EO.constraint(x -> x[i] - 45, <=, 0.))
	end
	push!(simple_problems, Problem("g04", f_1, constraints, 5, [78.0, 33.0, 29.9952560256815985, 45.0, 36.7758129057882073], −3.066553867178332e4))

	f_2(x) = 3*x.sequence[1] + 0.000001*x.sequence[1]^3 + 2*x.sequence[2] + (0.000002/3)*x.sequence[2]^3
	constraints = [
		EO.constraint(x -> -x[4] + x[3] - 0.55, <=, 0.),
		EO.constraint(x -> -x[3] + x[4] - 0.55, <=, 0.),
		EO.constraint(x -> 1000*sin(-x[3] - 0.25) + 1000*sin(-x[4] - 0.25) + 894.8 - x[1], <=, 0.),
		EO.constraint(x -> -(1000*sin(-x[3] - 0.25) + 1000*sin(-x[4] - 0.25) + 894.8 - x[1]), <=, 0.),
		EO.constraint(x -> 1000*sin(x[3] - 0.25) + 1000*sin(x[3] - x[4] - 0.25) + 894.8 - x[2], <=, 0.),
		EO.constraint(x -> -(1000*sin(x[3] - 0.25) + 1000*sin(x[3] - x[4] - 0.25) + 894.8 - x[2]), <=, 0.),
		EO.constraint(x -> 1000*sin(x[4] - 0.25) + 1000*sin(x[4] - x[3] - 0.25) + 1294.8, <=, 0.),
		EO.constraint(x -> -(1000*sin(x[4] - 0.25) + 1000*sin(x[4] - x[3] - 0.25) + 1294.8), <=, 0.),
		EO.constraint(x -> -x[1], <=, 0.),
		EO.constraint(x -> x[1] - 1200, <=, 0.),
		EO.constraint(x -> -x[2], <=, 0.),
		EO.constraint(x -> x[2] - 1200, <=, 0.),
		EO.constraint(x -> -x[3] - 0.55, <=, 0.),
		EO.constraint(x -> x[3] - 0.55, <=, 0.),
		EO.constraint(x -> -x[4] - 0.55, <=, 0.),
		EO.constraint(x -> x[4] - 0.55, <=, 0.)
	]
	push!(simple_problems, Problem("g05", f_2, constraints, 4, [679.945148297028709, 1026.06697600004691, 0.118876369094410433, -0.39623348521517826], 5126.4967140071))

	f_3(x) = (x.sequence[1]-10)^2 + 5*(x.sequence[2]-12)^2 + x.sequence[3]^4 + 3*(x.sequence[4]-11)^2 + 10*x.sequence[5]^6 + 7*x.sequence[6]^2 + x.sequence[7]^4 - 4*x.sequence[6]*x.sequence[7] - 10*x.sequence[6] - 8*x.sequence[7]
	constraints = [
		EO.constraint(x -> -127 + 2*x[1]^2 + 3*x[2]^4 + x[3] + 4*x[4]^2 + 5*x[5], <=, 0.),
		EO.constraint(x -> -282 + 7*x[1] + 3*x[2] + 10*x[3]^2 + x[4] - x[5], <=, 0.),
		EO.constraint(x -> -196 + 23*x[1] + x[2]^2 + 6*x[6]^2 - 8*x[7], <=, 0.),
		EO.constraint(x -> 4*x[1]^2 + x[2]^2 - 3*x[1]*x[2] + 2*x[3]^2 + 5*x[6] - 11*x[7], <=, 0.)
	]
	for i in 1:7
		push!(constraints, EO.constraint(x -> -x[i] - 10, <=, 0.))
		push!(constraints, EO.constraint(x -> x[i] - 10, <=, 0.))
	end
	push!(simple_problems, Problem("g09", f_3, constraints, 7, [2.33049935147405174, 1.95137236847114592, -0.477541399510615805, 4.36572624923625874, -0.624486959100388983, 1.03813099410962173, 1.5942266780671519], 680.630057374402))

	f_4(x) = x.sequence[1]
	constraints = [
		EO.constraint(x -> -x[1] + 35*real(complex(x[2])^0.6) + 35*real(complex(x[3])^0.6), <=, 0.),
		EO.constraint(x -> -300*x[3] + 7500*x[5] - 7500*x[6] - 25*x[4]*x[5] + 25*x[4]*x[6] + x[3]*x[4], <=, 0.),
		EO.constraint(x -> 300*x[3] - 7500*x[5] + 7500*x[6] + 25*x[4]*x[5] - 25*x[4]*x[6] - x[3]*x[4], <=, 0.),
		EO.constraint(x -> 100*x[2] + 155.365*x[4] + 2500*x[7] - x[2]*x[4] - 25*x[4]*x[7] - 15536.5, <=, 0.),
		EO.constraint(x -> -100*x[2] - 155.365*x[4] - 2500*x[7] + x[2]*x[4] + 25*x[4]*x[7] + 15536.5, <=, 0.),
		EO.constraint(x -> -x[5] + real(log(complex(-x[4] + 900))), <=, 0.),
		EO.constraint(x -> x[5] - real(log(complex(-x[4] + 900))), <=, 0.),
		EO.constraint(x -> -x[6] + real(log(complex(x[4] + 300))), <=, 0.),
		EO.constraint(x -> x[6] - real(log(complex(x[4] + 300))), <=, 0.),
		EO.constraint(x -> -x[7] + real(log(complex(-2*x[4] + 700))), <=, 0.),
		EO.constraint(x -> x[7] - real(log(complex(-2*x[4] + 700))), <=, 0.),
		EO.constraint(x -> -x[1], <=, 0.),
		EO.constraint(x -> x[1] - 1000, <=, 0.),
		EO.constraint(x -> -x[2], <=, 0.),
		EO.constraint(x -> x[2] - 40, <=, 0.),
		EO.constraint(x -> -x[3], <=, 0.),
		EO.constraint(x -> x[3] - 40, <=, 0.),
		EO.constraint(x -> -x[4] + 100, <=, 0.),
		EO.constraint(x -> x[4] - 300, <=, 0.),
		EO.constraint(x -> -x[5] + 6.3, <=, 0.),
		EO.constraint(x -> x[5] - 6.7, <=, 0.),
		EO.constraint(x -> -x[6] + 5.9, <=, 0.),
		EO.constraint(x -> x[6] - 6.4, <=, 0.),
		EO.constraint(x -> -x[7] + 4.5, <=, 0.),
		EO.constraint(x -> x[7] - 6.25, <=, 0.)
	]
	push!(simple_problems, Problem("g21", f_4, constraints, 7, [193.724510070034967, 5.56944131553368433e-27, 17.3191887294084914, 100.047897801386839, 6.68445185362377892, 5.99168428444264833, 6.21451648886070451], 93.724510070035))

	return simple_problems
end

function stochastic_ranking(p::T, constraints::Vector{Constraint}, P::Float64) where T<:Population
	gs = Gs(p, constraints)
	fs = p.fitness
	rank = collect(1:p.size)

	for _ in 1:p.size
		swapped = false
		for i in 1:p.size-1
			#swapped = false
			if (gs[rank[i]] == gs[rank[i+1]] && gs[rank[i]] == 0) || rand() < P         # compare by fitness
				if fs[rank[i]] > fs[rank[i+1]]
					rank[i], rank[i+1] = rank[i+1], rank[i]
					swapped = true
				end
			else                                                                        # compare by violation
				if gs[rank[i]] > gs[rank[i+1]]
					rank[i], rank[i+1] = rank[i+1], rank[i]
					swapped = true
				end
			end
		end
		if swapped == false
			break
		end
	end

	return rank
end

function save_result(path::String, name::String, i::Int, r::Result_big)
	mkpath(path)
	open(path*name*"_"*string(i), "w") do txt
		println(txt, string(r.top_coords))
		println(txt, r.top_value)
		println(txt, r.coords_history)
		println(txt, r.value_history)
		println(txt, r.pop_history)
		println(txt, r.penalty_history)
	end
end

function save_result(path::String, name::String, i::Int, r::Result_multi)
	mkpath(path)
	open(path*name*"_"*string(i), "w") do txt
		println(txt, r.coords_history)
		println(txt, r.pop_history)
	end
end

function denoise(xx::Vector{Float64}; λ=0.1, β=10.0)
    x = copy(xx)
    mean = signal_mean(x; λ=λ)
    diff = abs.(x-mean)
    x[diff .>= β] = mean[diff .>= β]
    return x
end

function signal_mean(x::Vector{Float64}; λ=0.1)
	mean = similar(x)
	mean[1] = x[1]
	for i in 2:length(x)
		if isnan(x[i])
			if i == 1
				x[i] = 0.0
			else
				x[i] = x[i-1]
			end
		end
		mean[i] = (1-λ)*mean[i-1] + λ*x[i]
	end
	return mean
end

function protected_div(x, y)
    if y == 0
        return x/1e-14
    end
    return x/y
end
function square(x)
    return x^2
end
function cube(x)
    return x^3
end
logaritmus(x) = log(abs(x))
exponenciala(x) = exp(abs(x))
power(x, y) = min(abs(x)^y, 1e10)

function enclose_arguments(f::Function, a...)::Function
	return x -> f(x, a...)
end

function enclose_argument(f::Function, a)::Function
	return (x...) -> f(a, x...)
end

function enclose_noargs(f::Function, a...)::Function
    return ()-> f(a...)
end

identity(x) = x