using Documenter
using EO

makedocs(
    sitename = "SPJ_TEST",
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

deploydocs(
    #repo = "https://github.com/BipoaroXigen/SPJ_TEST",  # Replace with your repo URL
    #branch = "gh-pages",  # This is the branch where the documentation will be hosted
    target = "build"      # The build directory where documentation is generated
)