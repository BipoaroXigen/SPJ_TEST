function plot_results(res::T) where T<:Result
    p = plot(layout=(2,1), legend=false)
    plot!(p, res.value_history, subplot=1, title="Objective")
    x = map(x->x[2][1], enumerate(res.coords_history))
    y = map(x->x[2][2], enumerate(res.coords_history))
    plot!(p, x, y, subplot=2, title="Coordinates")
end

function plot_results(res::T, v::Vector{Vector{Float64}}) where T<:Result

    p = plot(layout=(1,3), legend=false)

    order = res.coords_history[1]
    x = map(x->x[2][1], enumerate(v[order]))
    y = map(x->x[2][2], enumerate(v[order]))

    plot!(p, x, y, lc = :buda, subplot=1, title="initialization")

    order = res.coords_history[round(Int, length(res.coords_history)*0.5)]
    x = map(x->x[2][1], enumerate(v[order]))
    y = map(x->x[2][2], enumerate(v[order]))

    plot!(p, x, y, lc = :buda, subplot=2, title="middle iteration")

    order = res.top_coords
    x = map(x->x[2][1], enumerate(v[order]))
    y = map(x->x[2][2], enumerate(v[order]))

    plot!(p, x, y , lc = :buda, subplot=3, size=(1200, 400), title="final")
    plot!(p, [x[1], x[end]], [y[1], y[end]], lc = :green, ls=:dot, lw=3, subplot=3, size=(1200, 400), title="final")
end

function plot_graph_sequence(order#=::Vector{Int}=#, v::Vector{Vector{Float64}})
    x = map(x->x[2][1], enumerate(v[order]))
    y = map(x->x[2][2], enumerate(v[order]))
    p = plot(x, y, lc = :buda)
    plot!(p, [x[1], x[end]], [y[1], y[end]], lc = :green, ls=:dot, lw=3)
end
#=
function parse_benchmark_file(name::String)::Tuple{Vector{Vector{Float64}}, Vector{Float64}}
    coord_history = Vector{Vector{Float64}}()
    value_history = Vector{Float64}()
    open(name, "r") do f
        coords = false
        values = false
        while !eof(f)
            line = readline(f)
            if occursin("COORDS", line)
                coords = true
                values = false
                continue
            end
            if occursin("VALUES", line)
                coords = false
                values = true
                continue
            end
            if coords == true
                coord = parse.(Float64, (split(line, ", ")))
                push!(coord_history, coord)
            end
            if values == true
                value = parse.(Float64, line)
                push!(value_history, value)
            end
        end
    end
    return coord_history, value_history
end


"""Average performance of LS on each graph"""
function plot_TSP_benchmark(prefix::String, N_graphs::Int, N_runs::Int)
    """for given LS/ES display average performance on each TSP graph with std"""

    avgs = Vector{Vector{Float64}}()
    stds = Vector{Vector{Float64}}()

    _, values = parse_benchmark_file(prefix*"_G"*string(1)*"_"*string(1))
    value_history = copy(values)
    for run_i in 2:N_runs
        _, v_h = parse_benchmark_file(prefix*"_G"*string(1)*"_"*string(run_i))
        values = hcat(values, v_h)
    end
    value_history = hcat(value_history, copy(values))
    push!(avgs, vec(mean(values, dims=2)))
    push!(stds, vec(std(values,  dims=2)))

    for graph_i in 2:N_graphs
        _, values = parse_benchmark_file(prefix*"_G"*string(graph_i)*"_"*string(1))
        for run_i in 2:N_runs
            _, v_h = parse_benchmark_file(prefix*"_G"*string(graph_i)*"_"*string(run_i))
            values = hcat(values, v_h)
        end
        value_history = hcat(value_history, copy(values))
        push!(avgs, vec(mean(values, dims=2)))
        push!(stds, vec(std(values,  dims=2)))
    end

    average_run = mean(value_history, dims=2)
    std_run     = std(value_history, dims=2)
    return avgs, stds, average_run, std_run
end

function strahov_plot(x, y, i::Int)
    beetle_string = "M100 40 C90 20, 110 20, 100 40 M80 50 C70 70, 130 70, 120 50 M70 70 C50 100, 150 100, 130 70M65 100 C65 140, 135 140, 135 100 M80 140 L70 170 M120 140 L130 170M80 50 C70 30, 65 10, 85 10 M120 50 C130 30, 135 10, 115 10"

    batsymbol = BezierPath(beetle_string, fit = true, flipy = true)
    f = Figure()
    ax = Axis(f[1, 1], xlabel = "x_1", ylabel = "x_2", title = "Strahov plot of population in iteration "*string(i))
    CairoMakie.scatter!(ax, x, y, marker = batsymbol, color = :red,
                    markersize = range(1, 20, length = length(x)),
                    rotation = range(0, 2pi, length = length(x)+1)[1:end-1], strokewidth = 1, strokecolor = :black)
    return f
end

function strahov_plot(solution, i::Int; x_i=1, y_i=2)
    enum = enumerate(solution[i].population)
    x = map(x->x[2].sequence[x_i], enum)
    y = map(x->x[2].sequence[y_i], enum)
    beetle_string = "M100 40 C90 20, 110 20, 100 40 M80 50 C70 70, 130 70, 120 50 M70 70 C50 100, 150 100, 130 70M65 100 C65 140, 135 140, 135 100 M80 140 L70 170 M120 140 L130 170M80 50 C70 30, 65 10, 85 10 M120 50 C130 30, 135 10, 115 10"

    batsymbol = BezierPath(beetle_string, fit = true, flipy = true)
    f = Figure()
    ax = Axis(f[1, 1], xlabel = "x_1", ylabel = "x_2", title = "Strahov plot of population in iteration "*string(i))
    CairoMakie.scatter!(ax, x, y, marker = batsymbol, color = :red,
                    markersize = range(1, 20, length = length(x)),
                    rotation = range(0, 2pi, length = length(x)+1)[1:end-1], strokewidth = 1, strokecolor = :black)
    return f
end

function plot_pareto_fronts(solution::T) where T<:Population
    p = Plots.scatter(xlabel="objective", ylabel="violation", title="pareto front, and other fronts")

    for i in 1:maximum(EO.get_domination_count(solution.fitness))
        idx = findall(EO.get_domination_count(solution.fitness).==i)
        x = map(x->x[2][1], enumerate(solution.fitness[idx]))
        y = map(x->x[2][2], enumerate(solution.fitness[idx]))
        p = Plots.scatter!(p, x, y, label="losers", color=:red)
    end
    
    idx = findall(EO.get_domination_count(solution.fitness).==0)
    x = map(x->x[2][1], enumerate(solution.fitness[idx]))
    y = map(x->x[2][2], enumerate(solution.fitness[idx]))
    p = Plots.scatter!(p, x, y, label="Pareto front", color=:yellow)
    
    p
end

function plot_results(solution::Result_big; x_i=1, y_i=2)
    display("here")
    beetle_string = "M100 40 C90 20, 110 20, 100 40 M80 50 C70 70, 130 70, 120 50 M70 70 C50 100, 150 100, 130 70M65 100 C65 140, 135 140, 135 100 M80 140 L70 170 M120 140 L130 170M80 50 C70 30, 65 10, 85 10 M120 50 C130 30, 135 10, 115 10"

    batsymbol = BezierPath(beetle_string, fit = true, flipy = true)
    f = Figure(size = (1200, 900))

    enum = enumerate(solution.pop_history[1].population)
    x = map(x->x[2].sequence[x_i], enum)
    y = map(x->x[2].sequence[y_i], enum)
    ax = Axis(f[1, 1], xlabel = "x_1", ylabel = "x_2", title = "Initial population")
    CairoMakie.scatter!(ax, x, y, marker = batsymbol, color = :red,
                    markersize = range(1, 20, length = length(x)),
                    rotation = range(0, 2pi, length = length(x)+1)[1:end-1], strokewidth = 1, strokecolor = :black)

    enum = enumerate(solution.pop_history[round(Int, length(solution.pop_history)/2)].population)
    x = map(x->x[2].sequence[x_i], enum)
    y = map(x->x[2].sequence[y_i], enum)
    ax = Axis(f[1, 2], xlabel = "x_1", ylabel = "x_2", title = "Population in iteration "*string(round(Int, length(solution.pop_history)/2)))
    CairoMakie.scatter!(ax, x, y, marker = batsymbol, color = :red,
                    markersize = range(1, 20, length = length(x)),
                    rotation = range(0, 2pi, length = length(x)+1)[1:end-1], strokewidth = 1, strokecolor = :black)

    enum = enumerate(solution.pop_history[end].population)
    x = map(x->x[2].sequence[x_i], enum)
    y = map(x->x[2].sequence[y_i], enum)
    ax = Axis(f[1, 3], xlabel = "x_1", ylabel = "x_2", title = "Final population")
    CairoMakie.scatter!(ax, x, y, marker = batsymbol, color = :red,
                    markersize = range(1, 20, length = length(x)),
                    rotation = range(0, 2pi, length = length(x)+1)[1:end-1], strokewidth = 1, strokecolor = :black)

    if maximum(solution.value_history) >= 1e20
        ax = Axis(f[2, 1:3], xlabel = "obj. f. calls", ylabel = "f(x)", title = "Best solution so far", yscale = CairoMakie.log)
        CairoMakie.lines!(ax, solution.value_history, label="f(x)", color=:green)
    else
        ax = Axis(f[2, 1:3], xlabel = "obj. f. calls", ylabel = "f(x)", title = "Best solution so far"#= , yscale = CairoMakie.log =#)
        CairoMakie.lines!(ax, solution.value_history, label="f(x)", color=:green)
    end
    if maximum(solution.penalty_history) >= 1e20
        ax = Axis(f[3, 1:3], xlabel = "obj. f. calls", ylabel = "f(x)", title = "Best solution so far with penalty", yscale = CairoMakie.log)
        CairoMakie.lines!(ax, solution.value_history, label="f(x)", color=:green)
        CairoMakie.lines!(ax, solution.penalty_history #= + solution.value_history =#, label="f(x) + penalty", color=:red)
    else
        ax = Axis(f[3, 1:3], xlabel = "obj. f. calls", ylabel = "f(x)", title = "Best solution so far with penalty"#= , yscale = CairoMakie.log =#)
        CairoMakie.lines!(ax, solution.value_history, label="f(x)", color=:green)
        CairoMakie.lines!(ax, solution.penalty_history #= + solution.value_history =#, label="f(x) + penalty", color=:red)
    end

    return f
end

function plot_results(r::Result_multi)
    plot_pareto_fronts(r.pop_history[end])
end

using FileIO, JLD2

"""plot results on g06,g08,g11,g24"""
function plot_easy_benchmark(algs::Vector{String}, solutions::Vector{Float64})
    f = Figure(size = (1200, 900))

    j = 1
    ax = Axis(f[1, 1], title = "g06 fitness")
    for alg in algs
        average_value_history = zeros(101)
        for i in 1:100
            print("\r", alg, " ", i)
            solution = load(joinpath(@__DIR__, "..", "benchmarks", "constr", alg, "ez", "problem_"*string(j)*"_"*string(i)*".jld2"))["solution"]
            average_value_history += solution.f_history 
        end
        average_value_history /= 100
        CairoMakie.lines!(ax, average_value_history)
    end
    CairoMakie.hlines!(ax, [solutions[j]], label="solution", linewidth=2, color=:black)

    j = 2
    ax = Axis(f[2, 1], title = "g08 fitness")
    CairoMakie.ylims!(ax, -1, 1)
    for alg in algs
        average_value_history = zeros(101)
        for i in 1:100
            solution = load(joinpath(@__DIR__, "..", "benchmarks", "constr", alg, "ez", "problem_"*string(j)*"_"*string(i)*".jld2"))["solution"]
            average_value_history += solution.f_history 
        end
        average_value_history /= 100
        CairoMakie.lines!(ax, average_value_history)
    end
    CairoMakie.hlines!(ax, [solutions[j]], label="solution", linewidth=2, color=:black)

    j = 3
    ax = Axis(f[3, 1], title = "g11 fitness")
    for alg in algs
        average_value_history = zeros(101)
        for i in 1:100
            solution = load(joinpath(@__DIR__, "..", "benchmarks", "constr", alg, "ez", "problem_"*string(j)*"_"*string(i)*".jld2"))["solution"]
            average_value_history += solution.f_history 
        end
        average_value_history /= 100
        CairoMakie.lines!(ax, average_value_history)
    end
    CairoMakie.hlines!(ax, [solutions[j]], label="solution", linewidth=2, color=:black)

    j = 4
    ax = Axis(f[4, 1], xlabel = "generation", title = "g24 fitness")
    for alg in algs
        average_value_history = zeros(101)
        for i in 1:100
            solution = load(joinpath(@__DIR__, "..", "benchmarks", "constr", alg, "ez", "problem_"*string(j)*"_"*string(i)*".jld2"))["solution"]
            average_value_history += solution.f_history 
        end
        average_value_history /= 100
        CairoMakie.lines!(ax, average_value_history)
    end
    CairoMakie.hlines!(ax, [solutions[j]], label="solution", linewidth=2, color=:black)

    j = 1
    ax = Axis(f[1, 2], title = "g06 penalty", yscale = CairoMakie.log)
    for alg in algs
        average_value_history = zeros(101)
        for i in 1:100
            print("\r", alg, " ", i)
            solution = load(joinpath(@__DIR__, "..", "benchmarks", "constr", alg, "ez", "problem_"*string(j)*"_"*string(i)*".jld2"))["solution"]
            average_value_history += solution.p_history 
        end
        average_value_history /= 100
        CairoMakie.lines!(ax, average_value_history)
    end

    j = 2
    ax = Axis(f[2, 2], title = "g08 penalty", yscale = CairoMakie.log)
    for alg in algs
        average_value_history = zeros(101)
        for i in 1:100
            solution = load(joinpath(@__DIR__, "..", "benchmarks", "constr", alg, "ez", "problem_"*string(j)*"_"*string(i)*".jld2"))["solution"]
            average_value_history += solution.p_history 
        end
        average_value_history /= 100
        CairoMakie.lines!(ax, average_value_history)
    end

    j = 3
    ax = Axis(f[3, 2], title = "g11 penalty", yscale = CairoMakie.log)
    for alg in algs
        average_value_history = zeros(101)
        for i in 1:100
            solution = load(joinpath(@__DIR__, "..", "benchmarks", "constr", alg, "ez", "problem_"*string(j)*"_"*string(i)*".jld2"))["solution"]
            average_value_history += solution.p_history 
        end
        average_value_history /= 100
        CairoMakie.lines!(ax, average_value_history)
    end

    j = 4
    ax = Axis(f[4, 2], xlabel = "generation", title = "g24 penalty", yscale = CairoMakie.log)
    for alg in algs
        average_value_history = zeros(101)
        for i in 1:100
            solution = load(joinpath(@__DIR__, "..", "benchmarks", "constr", alg, "ez", "problem_"*string(j)*"_"*string(i)*".jld2"))["solution"]
            average_value_history += solution.p_history 
        end
        average_value_history /= 100
        CairoMakie.lines!(ax, average_value_history, label=alg)
    end

    f[1:4, 3] = Legend(f, ax, "Algorithms", framevisible = false)


    return f
end

"""plot results on g04,g05,g09,g21"""
function plot_hard_benchmark(algs::Vector{String}, solutions::Vector{Float64})
    f = Figure(size = (1200, 900))


    j = 1
    ax = Axis(f[1, 1], title = "g04 fitness")
    for alg in algs
        average_value_history = zeros(101)
        for i in 1:100
            solution = load(joinpath(@__DIR__, "..", "benchmarks", "constr", alg, "hard", "problem_"*string(j)*"_"*string(i)*".jld2"))["solution"]
            average_value_history += solution.f_history 
        end
        average_value_history /= 100
        CairoMakie.lines!(ax, average_value_history)
    end
    CairoMakie.hlines!(ax, [solutions[j]], label="solution", linewidth=2, color=:black)

    j = 2
    ax = Axis(f[2, 1], title = "g05 fitness")
    CairoMakie.ylims!(ax, 0, 1000)
    for alg in algs
        average_value_history = zeros(101)
        for i in 1:100
            solution = load(joinpath(@__DIR__, "..", "benchmarks", "constr", alg, "hard", "problem_"*string(j)*"_"*string(i)*".jld2"))["solution"]
            average_value_history += solution.f_history 
        end
        average_value_history /= 100
        CairoMakie.lines!(ax, average_value_history)
    end
    CairoMakie.hlines!(ax, [solutions[j]], label="solution", linewidth=2, color=:black)

    j = 3
    ax = Axis(f[3, 1], title = "g09 fitness")
    for alg in algs
        average_value_history = zeros(101)
        for i in 1:100
            solution = load(joinpath(@__DIR__, "..", "benchmarks", "constr", alg, "hard", "problem_"*string(j)*"_"*string(i)*".jld2"))["solution"]
            average_value_history += solution.f_history 
        end
        average_value_history /= 100
        CairoMakie.lines!(ax, average_value_history)
    end
    CairoMakie.hlines!(ax, [solutions[j]], label="solution", linewidth=2, color=:black)

    j = 4
    ax = Axis(f[4, 1], xlabel = "generation", title = "g21 fitness")
    for alg in algs
        average_value_history = zeros(101)
        for i in 1:100
            solution = load(joinpath(@__DIR__, "..", "benchmarks", "constr", alg, "hard", "problem_"*string(j)*"_"*string(i)*".jld2"))["solution"]
            average_value_history += solution.f_history 
        end
        average_value_history /= 100
        CairoMakie.lines!(ax, average_value_history)
    end
    CairoMakie.hlines!(ax, [solutions[j]], label="solution", linewidth=2, color=:black)

    j = 1
    ax = Axis(f[1, 2], title = "g04 penalty", yscale = CairoMakie.log)
    for alg in algs
        average_value_history = zeros(101)
        for i in 1:100
            solution = load(joinpath(@__DIR__, "..", "benchmarks", "constr", alg, "hard", "problem_"*string(j)*"_"*string(i)*".jld2"))["solution"]
            average_value_history += solution.p_history 
        end
        average_value_history /= 100
        CairoMakie.lines!(ax, average_value_history)
    end

    j = 2
    ax = Axis(f[2, 2], title = "g05 penalty", yscale = CairoMakie.log)
    for alg in algs
        average_value_history = zeros(101)
        for i in 1:100
            solution = load(joinpath(@__DIR__, "..", "benchmarks", "constr", alg, "hard", "problem_"*string(j)*"_"*string(i)*".jld2"))["solution"]
            average_value_history += solution.p_history 
        end
        average_value_history /= 100
        CairoMakie.lines!(ax, average_value_history)
    end

    j = 3
    ax = Axis(f[3, 2], title = "g09 penalty", yscale = CairoMakie.log)
    for alg in algs
        average_value_history = zeros(101)
        for i in 1:100
            solution = load(joinpath(@__DIR__, "..", "benchmarks", "constr", alg, "hard", "problem_"*string(j)*"_"*string(i)*".jld2"))["solution"]
            average_value_history += solution.p_history 
        end
        average_value_history /= 100
        CairoMakie.lines!(ax, average_value_history)
    end

    j = 4
    ax = Axis(f[4, 2], xlabel = "generation", title = "g21 penalty", yscale = CairoMakie.log)
    for alg in algs
        average_value_history = zeros(101)
        for i in 1:100
            solution = load(joinpath(@__DIR__, "..", "benchmarks", "constr", alg, "hard", "problem_"*string(j)*"_"*string(i)*".jld2"))["solution"]
            average_value_history += solution.p_history 
        end
        average_value_history /= 100
        CairoMakie.lines!(ax, average_value_history, label=alg)
    end
    
    f[1:4, 3] = Legend(f, ax, "Algorithms", framevisible = false)

    return f
end

"""plot results on g06,g08,g11,g24"""
function plot_easy_benchmark(j::Int, alg::String, s::Float64)
    f = Figure(size=(900, 450))
    
    ax = Axis(f[1, 1], xlabel = "generation", title = "g06 fitness")
    t = range(0, 101, length=101)
    avg = zeros(101, 100)
    for i in 1:100
        print("\r", alg, " ", i)
        solution = load(joinpath(@__DIR__, "..", "benchmarks", "constr", alg, "ez", "problem_"*string(j)*"_"*string(i)*".jld2"))["solution"]
        avg[:, i] = solution.f_history 
    end
    st = vec(std(avg, dims=2))
    avg = vec(mean(avg, dims=2))
    CairoMakie.lines!(ax, t, avg, label="found fitness")
    CairoMakie.hlines!(ax, [s], color=:black, label="solution")
    CairoMakie.band!(ax, t, avg - st, avg + st)
    axislegend(ax, framevisible = false)

    ax = Axis(f[1, 2], xlabel = "generation", title = "g06 penalty")
    avg = zeros(101, 100)
    for i in 1:100
        print("\r", alg, " ", i)
        solution = load(joinpath(@__DIR__, "..", "benchmarks", "constr", alg, "ez", "problem_"*string(j)*"_"*string(i)*".jld2"))["solution"]
        avg[:, i] = solution.p_history 
    end
    st = vec(std(avg, dims=2))
    avg = vec(mean(avg, dims=2))
    CairoMakie.lines!(ax, t, avg, label="penalty")
    CairoMakie.band!(ax, t, avg - st, avg + st)

    return f
end

"""plot results on g06,g08,g11,g24"""
function plot_hard_benchmark(j::Int, alg::String, s::Float64)
    f = Figure(size=(900, 450))
    
    ax = Axis(f[1, 1], xlabel = "generation", title = "g06 fitness")
    t = range(0, 101, length=101)
    avg = zeros(101, 100)
    for i in 1:100
        print("\r", alg, " ", i)
        solution = load(joinpath(@__DIR__, "..", "benchmarks", "constr", alg, "hard", "problem_"*string(j)*"_"*string(i)*".jld2"))["solution"]
        avg[:, i] = solution.f_history 
    end
    st = vec(std(avg, dims=2))
    avg = vec(mean(avg, dims=2))
    CairoMakie.lines!(ax, t, avg, label="found fitness")
    CairoMakie.hlines!(ax, [s], color=:black, label="solution")
    CairoMakie.band!(ax, t, avg - st, avg + st)
    axislegend(ax, framevisible = false)

    ax = Axis(f[1, 2], xlabel = "generation", title = "g06 penalty")
    avg = zeros(101, 100)
    for i in 1:100
        print("\r", alg, " ", i)
        solution = load(joinpath(@__DIR__, "..", "benchmarks", "constr", alg, "hard", "problem_"*string(j)*"_"*string(i)*".jld2"))["solution"]
        avg[:, i] = solution.p_history 
    end
    st = vec(std(avg, dims=2))
    avg = vec(mean(avg, dims=2))
    CairoMakie.lines!(ax, t, avg, label="penalty")
    CairoMakie.band!(ax, t, avg - st, avg + st)

    return f
end

=#