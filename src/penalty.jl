"""function which combines the used objective function and penalty function"""
penalized_objective(x::RealChromosome, f::Function, p::Function) = f.(x) + p.(x)

"""funciton which detects the violation of constraints"""
function Gs(x::T, constraints::Vector{Constraint}) where T<:Chromosome
    return sum([constraint(x) for constraint in constraints])       # sum of constraint violations for a solution
end
function G(x::T, constraints::Vector{Constraint}) where T<:Chromosome
    return [constraint(x) for constraint in constraints]            # vector of constraint violations for a solution
end

function Gs(p::Population, constraints::Vector{Constraint})
    return map(x->Gs(x, constraints), p.population)                 # sum of constraint violations for each solution
end
function G(p::Population, constraints::Vector{Constraint})
    return map(x->G(x, constraints), p.population)                  # vector of constraint violations for each solution
end

"""nonlinear penalty function"""
function p_non_linear(p::T, rg::Float64, constraints::Vector{Constraint})::Vector{Float64} where{T<:Population}
    return rg*Gs(p, constraints).^2
end

function p_non_linear(x::T, rg::Float64, constraints::Vector{Constraint}) where T<:Chromosome
    return rg*Gs(x, constraints)^2
end

function p_non_linear(p::T, rg::Vector{Float64}, constraints::Vector{Constraint})::Vector{Float64} where{T<:Population}
    return vec(sum(rg' * (reduce(hcat, G(p, constraints))).^2, dims=1))
end

function p_non_linear(x::T, rg::Vector{Float64}, constraints::Vector{Constraint}) where T<:Chromosome
    return sum(rg.*(G(x, constraints)).^2)
end

"""death penalty function"""
function death_penalty(p::T, rg::Float64, constraints::Vector{Constraint}, pencle)::Vector{Float64} where{T<:Population}
    return pencle.*(Gs(p, constraints) .> 0)
end

function death_penalty(x::T, rg::Float64, constraints::Vector{Constraint}, pencle) where T<:Chromosome
    return pencle*(Gs(x, constraints) > 0)
end

function death_penalty(p::T, rg::Vector{Float64}, constraints::Vector{Constraint}, pencle)::Vector{Float64} where{T<:Population}
    return pencle.*(Gs(p, constraints) .> 0)
end

function death_penalty(x::T, rg::Vector{Float64}, constraints::Vector{Constraint}, pencle) where T<:Chromosome
    return pencle*(Gs(x, constraints) > 0)
end


