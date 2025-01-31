# Data types and structures

The EO library relies on data structures defined in `src/structs.jl`

### Chromosome
Is a structure for the solutions in the population.
The solver expects subtypes of this type, allowing the use for defining his custom format of the solutions.
The library currently has its own implementationf of
```
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
```

### Population
Is a structure for handling the populations.
The solver expects subtypes of this type, allowing the use for defining his custom format of the population.
The library currently has its own implementationf of
```
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
```

### Result
A data structure returned by the solver
The solver expects subtypes of this type, allowing the use for defining his custom format of the result.

```
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

mutable struct Result_big <: Result
    top_coords::Any
    top_value::Float64
    coords_history::Vector{Any}
    value_history::Vector{Float64}
    pop_history::Vector{<:Population}
    penalty_history::Vector{Float64}
end
```

### Constraint
A data structure for handling constraints of objectives 

```
struct Constraint
    leftside::Function
    operator::Function
    constr_value::Float64
end
```

Examples of setting up the constraints can be seen in `src/utils.jl`.