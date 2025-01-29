using EO
using Test

f = x->0

@testset "EO.jl" begin
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
