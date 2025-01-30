using EO
using Test

f = x->0

#@testset "EO.jl" begin
@testset "basics" begin
    # Write your tests here.
    @test 1 == 1

    # create chromosomes
    @test all(EO.get_binary_chromosome(ones(10)).sequence .= 1)
    @test EO.get_real_chromosome(ones(10)).length == 10

    # objective functions
    @test EO.f_sphere(EO.get_real_chromosome(zeros(10)), zeros(10)) == 0

    # population initialization
    @test EO.binary_initialization(2, 100, f).size == 100
    @test EO.interval_real_initialization(2, 100, f, 99, 100).size == 100
end

dims = 10
pop_size = 100
objective_function  = enclose_arguments(f_sphere, zeros(dims))  # min at origin
initialization      = enclose_noargs(interval_real_initialization, 10, pop_size, objective_function, 0, 100)
selection           = enclose_arguments(EO.s_tournament, pop_size, 3)
crossover           = enclose_arguments(EO.cr_parent_sum, pop_size)
mutation            = enclose_arguments(gaussian_preturbation!, 0.25)
replacement         = EO.enclose_replacement(EO.r_keep_best_n, pop_size)
termination         = enclose_argument(iteration_termination, pop_size*100)
solution = solvink_hart(objective_function, initialization, selection, crossover, mutation, replacement, termination)

@testset "LS" begin
    
    @test solution.top_value < 100

end
