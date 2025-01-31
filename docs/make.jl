using Documenter
using EO

makedocs(
    sitename = "EO.jl",
    #modules = [EO],
    #format = Documenter.HTML(),  # Output format
    pages = [
        "Home" => "index.md",   # Main page
        "Introduction and Operator functions" => "operators.md",
        "Data types and structures" => "types.md",
        "Simple solver examples" => "solver.md",
        "Traveling salesman problem" => "tsp.md",
        "Solvers for constrained optimization" => "constrains.md"
    ]
)