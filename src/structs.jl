
abstract type Chromosome end
abstract type Result end
abstract type Population end

mutable struct BinaryChromosome <: Chromosome
    sequence::Vector{Bool}
    length::Int
end
mutable struct RealChromosome <: Chromosome
    sequence::Vector{Real}
    length::Int
end
mutable struct ExprChromosome <: Chromosome
    sequence::Expr
    length::Int
end

mutable struct SingleObjPopulation{T<:Chromosome} <: Population
    population::Vector{T}
    fitness::Vector{<:Real}
    penalty::Vector{Float64}
    size::Int
end

mutable struct MultiObjPopulation{T<:Chromosome} <: Population
    population::Vector{T}
    fitness::Vector{Vector{<:Real}}
    penalty::Vector{Float64}
    size::Int
end

mutable struct Result_real <: Result
    top_coords::Vector{<:Real}
    top_value::Float64
    coords_history::Vector{Vector{<:Real}}
    value_history::Vector{Float64}
end

mutable struct Result_binary <: Result
    top_coords::Vector{Bool}
    top_value::Float64
    coords_history::Vector{Vector{Bool}}
    value_history::Vector{Float64}
end

struct Constraint
    leftside::Function
    operator::Function
    constr_value::Float64
end

mutable struct Result_big <: Result
   # top_coords::Vector{<:Real}
    top_coords::Any
    top_value::Float64
    #coords_history::Vector{Vector{<:Real}}
    coords_history::Vector{Any}
    value_history::Vector{Float64}
    pop_history::Vector{<:Population}
    penalty_history::Vector{Float64}
end

mutable struct Result_multi <: Result
    #top_coords::Vector{<:Real}
    top_coords::Any
    top_value::Vector{Float64}
    #coords_history::Vector{Vector{<:Real}}
    coords_history::Vector{Any}
    value_history::Vector{Float64}
    pop_history::Vector{<:Population}
    penalty_history::Vector{Float64}
end

struct Result_benchmark <: Result
    f_history::Vector{Float64}
    p_history::Vector{Float64}
end

struct NConstraint end

struct NFunction end

struct Problem
	name::String
    objective::Function
    constraints::Vector{Constraint}
	dimension::Int
	solution::Vector{Float64}
	solution_v::Float64
end

# allows user to instantiate chromosome
function get_binary_chromosome(seq::Vector{})::BinaryChromosome
    return BinaryChromosome(seq, length(seq))
end

function get_binary_chromosome(seq::BitVector)::BinaryChromosome
    return BinaryChromosome(seq, length(seq))
end

function get_binary_chromosome(seq::Matrix{})::BinaryChromosome
    return BinaryChromosome(vec(seq), length(seq))
end


function get_real_chromosome(seq::Vector{<:Real})::RealChromosome
    return RealChromosome(seq, length(seq))
end

function get_real_chromosome(seq::Matrix{<:Real})::RealChromosome
    return RealChromosome(vec(seq), length(seq))
end

function get_expr_chromosome(seq::Expr)::ExprChromosome
    return ExprChromosome(seq, length(seq.args))
end

function get_expr_chromosome(seq::Int)::ExprChromosome
    return ExprChromosome(:($seq), 1)
end

# for automatical instantiation by functions

function BinaryChromosome(seq::Vector{<:Real})::BinaryChromosome
    return BinaryChromosome(Bool.(seq), length(seq))
end

function RealChromosome(seq::Vector{<:Real})::RealChromosome
    return RealChromosome(seq, length(seq))
end

function ExprChromosome(seq::Expr)::ExprChromosome
    return ExprChromosome(seq, length(seq.args))
end

function Result(x::RealChromosome, f, x_history, f_history)
    return Result_real(x_best.sequence, f, x_history, f_history)
end
function Result(x::BinaryChromosome, f, x_history, f_history)
    return Result_binary(x_best.sequence, f, x_history, f_history)
end

#= function Population(pop::Vector{T}, f::Function)::Population where T<:Chromosome
    return Population(pop, f.(pop), length(pop))
end =#

"""do I even use this?"""
#""" yes i do yes i do yes i do"""
function MultiObjPopulation(pop::Vector{T})::MultiObjPopulation where T<:Chromosome
    return MultiObjPopulation{T}(pop, [zeros(length(pop))], zeros(length(pop)), length(pop))
end

function Population(pop::Vector{T}, fitness::Vector{Vector{<:Real}}, size::Int)::MultiObjPopulation{T} where T<:Chromosome
    return MultiObjPopulation{T}(pop, fitness, zeros(size), size)
end

function Population(pop::Vector{T}, fitness::Vector{Vector{Float64}}, size::Int)::MultiObjPopulation{T} where T<:Chromosome
    return MultiObjPopulation{T}(pop, fitness, zeros(size), size)
end

function Population(pop::Vector{T}, fitness::Vector{Float64}, size::Int)::SingleObjPopulation{T} where T<:Chromosome
    return SingleObjPopulation{T}(pop, fitness, zeros(size), size)
end

function Population(pop::Vector{T}, fitness::Vector{<:Real}, size::Int)::SingleObjPopulation{T} where T<:Chromosome
    return SingleObjPopulation{T}(pop, fitness, zeros(size), size)
end

function Population(pop::Vector{T})::SingleObjPopulation{T} where T<:Chromosome
    return SingleObjPopulation{T}(pop, zeros(length(pop)), zeros(length(pop)), length(pop))
end

function Base.copy(x::BinaryChromosome)::BinaryChromosome
    return BinaryChromosome(Base.copy(x.sequence), length(x.sequence))
end

function Base.copy(x::RealChromosome)::RealChromosome
    return RealChromosome(Base.copy(x.sequence), length(x.sequence))
end

function Base.copy(x::ExprChromosome)::ExprChromosome
    return ExprChromosome(Base.copy(x.sequence), x.length)
end

get_ast_len(e::Expr) = sum(get_ast_len.(e.args))
get_ast_len(e::Any) = 1

#constraint(leftside<:Function, operator::Function, constr_value::Float64) = Constraint(leftside, operator, constr_value)
constraint(leftside#= <:Function =#, operator#= ::Function =#, constr_value::Float64) = Constraint(leftside, operator, constr_value)

"""return 1 if any constraint is violated"""
#(c::Constraint)(x::T) where T<:Chromosome = !all(c.operator(c.leftside(x.sequence), c.constr_value))
(c::Constraint)(x::T) where T<:Chromosome = max(c.leftside(x.sequence), 0)
(c::Constraint)(x::Vector{<:Chromosome}) = c.(x)
(cs::Vector{Constraint})(x::Vector{<:Chromosome}) = [c(x) for c in cs]
(cs::Vector{Constraint})(x::T) where T<:Chromosome = [c(x) for c in cs]

"""placeholder for unconstrained problems"""
(c::NConstraint)(x::T) where T<:Chromosome = 0
(c::NConstraint)(x::Vector{<:Chromosome}) = zeros(x.length)

"""placeholder for unpenalized problems"""
(c::NFunction)(::T, ::Float64, ::Union{NConstraint, Vector{Constraint}}) where T<:Chromosome = 0
(c::NFunction)(x::P, ::Float64, ::Union{NConstraint, Vector{Constraint}}) where P<:Population = zeros(x.size)
(c::NFunction)(::T, ::Vector{Float64}, ::Union{NConstraint, Vector{Constraint}}) where T<:Chromosome = 0
(c::NFunction)(x::P, ::Vector{Float64}, ::Union{NConstraint, Vector{Constraint}}) where P<:Population = zeros(x.size)

function Base.show(io::IO, x::Result_big)
    println(io, "top_x: ", x.top_coords)
    println(io, "top_f: ", x.top_value)
    println(io, "penalty: ", x.penalty_history[end])
end

function Base.show(io::IO, x::Result_multi)
    println(io, "top_x: ", x.top_coords)
    println(io, "top_f: ", x.top_value)
    println(io, "penalty: ", x.penalty_history[end])
end

function print_solution(x::Result_multi, cs::Vector{Constraint})
    println("top_x: ", x.top_coords)
    println("top_f: ", x.top_value)
    println("violations: ", G(get_real_chromosome(x.top_coords), cs))
end

function print_solution(x::Result_big, cs::Vector{Constraint})
    println("top_x: ", x.top_coords)
    println("top_f: ", x.top_value)
    println("violations: ", G(get_real_chromosome(x.top_coords), cs))
end

function print_solution(x::Result_multi, cs::Vector{Constraint}, p::Problem)
    println("top_x: ", x.top_coords)
    println("top_f: ", x.top_value)
    println("violations: ", G(get_real_chromosome(x.top_coords), cs))
end

function compare(x::T, p::Problem) where T<:Result
    println("top_x: ", x.top_coords)
    println("top_f: ", x.top_value)
    println("violations: ", G(get_real_chromosome(x.top_coords), p.constraints))
    println("--------------------------------")
    println("opt_x: ", p.solution_v)
    println("opt_f: ", p.solution)
end

### manual expression evaluation

struct Expr_parser
    e::Expr
    symbols::Dict{Symbol,Int}
end

function Expr_parser(e::Expr)
    symbols = Dict( s=>i for (i,s) in pairs(sort!(getvars(e))) )
    Expr_parser(e, symbols)
end

(e::Expr_parser)(values::T...) where T = evaluate(e.e, e.symbols, values...)