using EO
using Test

@testset "EO.jl" begin
    # Write your tests here.
    @test 1 == 1

    # create chromosomes
    @test isa(EO.BinaryChromosome, EO.get_binary_chromosome(ones(10)))

end
