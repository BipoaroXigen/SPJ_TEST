# Introduction and Operator functions

The solver evolves a population of solutions.
The solutions are initialized, modified and evolved accordingly to functions the solver obtains as agruments.
The specific types of functions are explained later.

The solvers behaviour can be flexibly changed by the operator functions provided by the user.
Operator functions in this context mean functions defining certain strategies for evolutionary algorithms.
The library contains several implemented operator functions, however the library is made with custom operator functions in mind, cince it is expected, the user will have ideas for functions defining his specific problem better.

---

## Initialization

Function which (usually) randomly creates a population of solutions.
Only requirement on the function is for it to return a population of a type which is subtype of `EO.Population`

Some initializations already present in the library are:

```
EO.interval_real_initialization(dims::Int, size::Int, f::Function, lb::T, ub::T) where T<:Real
EO.TSP_initialization(dims::Int, size::Int, f::Function)
EO.expression_initialization(size::Int, f::Function, basis_functions::Vector{Function}, basis_variables::Vector{Any})

```

## Selection

Function which selects from the current population parents of the next genertion.
Expected output is a vector indexing the selected parents from the population.

Some selections already present in the library are:

```
EO.s_tournament(p::SingleObjPopulation, out_n::Int, poll_size::Int)::Vector{Int}
EO.s_stochastic_tournament(p::Population, out_n::Int, poll_size::Int, constraints::Vector{Constraint}, P::Float64)::Vector{Int}
EO.s_greedy_overselection(p::SingleObjPopulation, out_n::Int, prc::Float64)::Vector{Int}

```

## Crossover

Function combining parents into new solutions.
The solver expects two functions, one implementing the creation of new solution from parent solutions, and other calling that function on the population `<:EO.Population` of parents
Expected output is population `<:EO.Population` of children

Some crossovers already present in the library are:

```
EO.cr_single_point(a::T, b::T)::T where T<:Chromosome
EO.cr_single_point(p<:Population)<:Population

EO.cr_ordered(a::T, b::T)::T where T<:Chromosome
EO.cr_ordered(p::Population)::Population

EO.cr_subtree(a::ExprChromosome, b::ExprChromosome, basis_functions::Vector{Function}, basis_variables::Vector{Any})::Tuple{EO.ExprChromosome, EO.ExprChromosome}
EO.cr_subtree(p::Population, out_n::Int, basis_functions::Vector{Function}, basis_variables::Vector{Any})::Population

```

## Replacement

Function combining current population with the children, into new population.
The function is expected to accept two populations and return one population

Some replacements already present in the library are:

```
EO.r_replacement(old::Population, new::Population)::Population = new
EO.r_merge(old<:Population, new<:Population)<:Population
r_keep_best_n(old::Population, new::Population, n::Int)::Population
```

## Mutation

Function randomly mutating single solution
Accepts solution as subtype of `Chromosome` and returns the same

Some mutations already present in the library are:

```
EO.basic_preturbation!(x::BinaryChromosome, p::Real)
EO.gaussian_preturbation!(x::RealChromosome, p::Real)
EO.pair_switch!(x::RealChromosome, G::Matrix{Float64})
```

## Objective function

Function returning a scalar or vector of real numbers, which the solver minimizes

Some objectives already present in the library are:

```
EO.f_onemax(x::BinaryChromosome)::Int = sum(x.sequence)
EO.f_sphere(x::RealChromosome, o::Vector{<:Real})::Float64 = sum((x.sequence.-o).^2)
EO.f_dist_sum(x::RealChromosome, G::Matrix{Float64})::Float64 = sum( [G[x.sequence[i], x.sequence[i+1]] for i in 1:x.length-1] ) + G[x.sequence[1], x.sequence[end]]
```