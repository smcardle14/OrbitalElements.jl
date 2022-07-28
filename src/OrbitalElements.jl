module OrbitalElements

using StaticArrays

"""
    coe2rv(coe, μ=1)

Convert classical orbital elements 'coe' to Cartesian position and velocity 'rv'.

Order classical orbital elements like so: [a, e, i, ω, Ω, ν]. 
See Algorithm 9 from Vallado, D. "Fundamentals of Astrodynamics and Applications." 
4th Edition. 2013.

See also [`rv2coe`](@ref)

"""
function coe2rv(coe, μ=1)
    a, e, i, ω, Ω, ν = coe

    # Special cases
    circular = abs(e) < 1.0e-16
    equatorial = abs(i) < 1.0e-16
    if circular && equatorial
        ω = 0.0
        Ω = 0.0
    elseif circular && !equatorial
        ω = 0.0
    elseif !circular && equatorial
        Ω = 0.0
    end

    p = a*(1-e^2)

    rpqw = @SVector([ cos(ν),   sin(ν), 0])*p/(1+e*(cos(ν)))
    vpqw = @SVector([-sin(ν), e+cos(ν), 0])*sqrt(μ/p)

    R = rot3(-Ω)*rot1(-i)*rot3(-ω)

    r = R*rpqw
    v = R*vpqw

    rv = vcat(r, v)

    return rv
end

function rv2coe(rv, μ=1)
    coe=rv
    return coe
end

function rot1(θ)
    return @SMatrix(
    [1       0      0;
     0  cos(θ) sin(θ);
     0 -sin(θ) cos(θ)])
end

function rot3(θ)
    return @SMatrix(
    [ cos(θ)  sin(θ)  0;
     -sin(θ)  cos(θ)  0;
          0       0   1])
end

export coe2rv

end
