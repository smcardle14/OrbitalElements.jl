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

    # coe2rv
    # Direct equatorial case (i=0)
    @test coe2rv([1.0, 0.0, 0.0, 0.0, 0.0, 0.0]) == [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]

    # Retrograde equatorial case (i=π)
    @test isapprox(coe2rv([1.0, 0.0, π, 0.0, 0.0, 0.0]), [1.0, 0.0, 0.0, 0.0, -1.0, 0.0], rtol=1e-15)

    # Direct equatorial case (i=2π)
    @test isapprox(coe2rv([1.0, 0.0, 2π, 0.0, 0.0, 0.0]), [1.0, 0.0, 0.0, 0.0, 1.0, 0.0], rtol=1e-15)

    # Vallado test case
    @test isapprox(coe2rv(coe, μ), rv, rtol=1e-15*μ)


    # rv2coe
    # Circular, polar case (i=90)
    @test isapprox(rv2coe([1.0, 0.0, 0.0, 0.0, 0.0, 1.0]), [1.0, 0.0, π/2, 0.0, 0.0, 0.0], rtol=1e-15)

    # Direct circular, equatorial case (i=0)
    @test rv2coe([1.0, 0.0, 0.0, 0.0, 1.0, 0.0]) == [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    # Direct eccentric, equatorial case (i=0)
    @test isapprox(rv2coe([1.0, 0.0, 0.0, 0.0, 0.9, 0.0]), [0.8403361344537815, 0.19, 0.0, 0.0, 0.0, π], rtol=1e-15)

    # Retrograde circular, equatorial case (i=π)
    @test isapprox(rv2coe([1.0, 0.0, 0.0, 0.0, -1.0, 0.0]), [1.0, 0.0, π, 0.0, 0.0, 0.0], rtol=1e-15)

    # Retrograde eccentric, equatorial case (i=π)
    @test isapprox(rv2coe([1.0, 0.0, 0.0, 0.0, -0.9, 0.0]), [0.8403361344537815, 0.19, π, 0.0, 0.0, π], rtol=1e-15)

    # Vallado test case
    @test isapprox(rv2coe(rv, μ), coe, rtol=1e-15*μ)


    # Benchmarking
    println("benchmarking coe2rv...")
    @btime coe2rv($coe, $μ)

    println("benchmarking rv2coe...")
    @btime rv2coe($rv, $μ)

end
