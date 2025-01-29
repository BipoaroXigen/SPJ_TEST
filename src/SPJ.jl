module EO

    using Random
    using Plots
    using Statistics
    #using CairoMakie

    include("structs.jl")
    include("utils.jl")
    include("plotting.jl")

    include("objective_functions.jl")
    include("initializations.jl")
    include("optimizers.jl")
    include("preturbations.jl")
    include("terminal_conditions.jl")
    include("selection.jl")
    include("crossover.jl")
    include("replacement.jl")
    include("penalty.jl")
    include("expressions.jl")
    include("FS.jl")

    export get_binary_chromosome, get_real_chromosome
    export enclose_arguments, enclose_argument, enclose_noargs, enclose_replacement
    export f_onemax, f_LABS, f_sphere, f_rosenbrock, f_linear, f_step, f_rastrigin, f_griewank, f_schwefel
    export binary_initialization, interval_real_initialization, TSP_initialization
    export first_improving_local_search, solvink_hart
    export basic_preturbation!, gaussian_preturbation!
    export time_termination, iteration_termination, progress_termination
    export s_tournament
    export cr_single_point
    export r_replacement, r_merge, r_keep_best_n
    export plot_TSP_benchmark

    export Chromosome, BinaryChromosome, RealChromosome, Constraint

end # module EO

println("compiled")