using OrbitalElements
using Test
using BenchmarkTools

deg2rad = π/180.0
μ = 398600.4418

# Vallado test case (Example 2-6)
p = 11067.790
e = 0.83285
coe = [p/(1-e^2), e, 87.87*deg2rad, 53.38*deg2rad, 227.89*deg2rad, 92.335*deg2rad]
rv = [6525.36812098609, 6861.531834896053, 6449.11861416016,
      4.902278646418963, 5.533139568361491, -1.975710099535108]

@testset "OrbitalElements.jl" begin

    @test coe2rv([1.0, 0.0, 0.0, 0.0, 0.0, 0.0]) == [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
    @test isapprox(coe2rv(coe, μ), rv, rtol=1e-15*μ)

    @test rv2coe([1.0, 0.0, 0.0, 0.0, 1.0, 0.0]) == [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    @test isapprox(rv2coe(rv, μ), coe, rtol=1e-15*μ)

    println("benchmarking coe2rv...")
    @btime coe2rv($coe, $μ)
    println("benchmarking rv2coe...")
    @btime rv2coe($rv, $μ)
end
