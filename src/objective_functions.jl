### BINARY
f_onemax(x::BinaryChromosome)::Int = sum(x.sequence)

function f_LABS(x::BinaryChromosome)::Float64
	fitness = 0
	s = bool_to_units(x.sequence)
	for k in 1:x.length-1
		fitness += sequence_correlation(s, k)^2
	end
	return fitness
end

f_weighted_sum(x::BinaryChromosome, w::Vector{T}) where {T<:Real} = Float64(sum(w[x.sequence]))

### REAL

f_sphere(x::RealChromosome, o::Vector{<:Real})::Float64 = sum((x.sequence.-o).^2)

function f_rosenbrock(x::RealChromosome)::Float64
	fitness = 0
	for i in 1:x.length-1
		fitness += 100*(x.sequence[i+1]−x.sequence[i]^2)^2+(1−x.sequence[i])^2
	end
	return fitness
end

f_linear(x::RealChromosome, a::Vector{<:Real})::Float64 = a[1] + sum(a[2:end].*x.sequence)

f_step(x::RealChromosome, a::Vector{<:Real})::Float64 = a[1] + sum(floor.(a[2:end].*x.sequence))

f_rastrigin(x::RealChromosome)::Float64 = 10*x.length + sum(x.sequence.^2 .- 10*cos.(2*pi.*x.sequence))

f_griewank(x::RealChromosome)::Float64 = 1 + 1/4000*sum(x.sequence.^2) - prod(cos.(x.sequence./sqrt.(collect(1:x.length))))

f_schwefel(x::RealChromosome)::Float64 = -sum(x.sequence.*sin.(sqrt.(abs.(x.sequence))))

### TSP

f_dist_sum(x::RealChromosome, G::Matrix{Float64})::Float64 = sum( [G[x.sequence[i], x.sequence[i+1]] for i in 1:x.length-1] ) + G[x.sequence[1], x.sequence[end]]

### MULTIOBjECTIVE

f_multi(x::RealChromosome, fs::Vector{Function})::vector{Float64} = [f(x.sequence) for f in fs]


mutable struct MultiObjFunction{T} <: Function
	#f::Vector{Any}
	f::Vector{>:T}
end
#get_multiobjective_function(f::Vector{Function}) = MultiObjFunction(f)
#get_multiobjective_function(f::Vector{Any}) = MultiObjFunction(f)
(fs::MultiObjFunction)(x::RealChromosome) = [f(x) for f in fs.f]
(fs::MultiObjFunction)(x::BinaryChromosome) = [f(x) for f in fs.f]
(fs::MultiObjFunction)(x::ExprChromosome)::Vector{Float64} = [f(x) for f in fs.f]

### symbolic regression

#f_function_diff(f::ExprChromosome, x::Vector{<:Real}, y::Vector{<:Real}) = abs(eval(f.sequence).(x) - y)
f_function_diff_squared(f::ExprChromosome, y::Vector{<:Real}, x::Vector{<:Real}...) = sum((Expr_parser(f.sequence).(x...) .- y).^2)

f_function_diff_abs(f::ExprChromosome, y::Vector{<:Real}, x::Vector{<:Real}) = sum(abs.(Expr_parser(f.sequence).(x...) .- y))

function f_function_diff_subset(f::ExprChromosome, y::Vector{<:Real}, r::Float64, x::Vector{<:Real}...)
	selection = rand(1:length(x[1]), max(ceil(Int, rand()*length(x[1])), ceil(Int, r*length(x[1]))))
	return sum((Expr_parser(f.sequence).(x[1][selection]) .- y[selection]).^2)
end

function f_function_diff_derivative(f::ExprChromosome, y::Vector{<:Real}, x::Vector{<:Real}...)
	res = 0.0
	for i in 1:length(x[1])-1
		dif = y[i+1]-y[i]
		decrease = dif >= 0
		dif_hat = Expr_parser(f.sequence).(x[i+1]...) - Expr_parser(f.sequence).(x[i]...)
		decrease_hat = dif_hat >= 0
		if decrease_hat != decrease
			res += 1.0
		end
	end
	return res
end

function f_function_diff_derivative_w(f::ExprChromosome, y::Vector{<:Real}, x::Vector{<:Real}...)

	error =  (Expr_parser(f.sequence).(x...) .- y).^2

	res = 0.0
	for i in 1:length(x[1])-1
		dif = y[i+1]-y[i]
		decrease = dif >= 0
		dif_hat = Expr_parser(f.sequence).(x[i+1]...) - Expr_parser(f.sequence).(x[i]...)
		decrease_hat = dif_hat >= 0
		if decrease_hat != decrease
			res += 1.0 * error[i]
		end
	end
	return res
end