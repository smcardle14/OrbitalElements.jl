using OrbitalElements
using Test

@testset "OrbitalElements.jl" begin
    @test coe2rv([1.0, 0.0, 0.0, 0.0, 0.0, 0.0]) == [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]

    # # Vallado test case (Example 2-6)
    # deg2rad = π/180.0
    # p = 11067.790
    # e = 0.83285
    # a = p/(1-e^2)
    # i = 87.87*deg2rad
    # ω = 53.38*deg2rad
    # Ω = 227.89*deg2rad
    # ν = 92.335*deg2rad
    # μ = 398600.4418
    # # More precision than Vallado's result
    # @test coe2rv([a, e, i, ω, Ω, ν], μ) == 
    #     [6525.36812098609, 6861.531834896053, 6449.11861416016,
    #      4.902278646418963, 5.533139568361491, -1.975710099535108]
end
