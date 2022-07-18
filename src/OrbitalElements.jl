module OrbitalElements

using StaticArrays

function coe2rv(coe, μ=1)
    a, e, i, ω, Ω, ν = coe

    p = a*(1-e^2)

    rpqw = @SVector([ cos(ν),   sin(ν), 0])*p/(1+e*(cos(ν)))
    vpqw = @SVector([-sin(ν), e+cos(ν), 0])*sqrt(μ/p)

    R = rot3(-Ω)*rot1(-i)*rot3(-ω)

    r = R*rpqw
    v = R*vpqw

    rv = vcat(r, v)

    return rv
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
