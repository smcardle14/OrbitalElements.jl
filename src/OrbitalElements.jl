module OrbitalElements

using StaticArrays


"""
    coe2rv(coe, μ=1)

Convert classical orbital elements 'coe' to Cartesian position and velocity 'rv'.

Input classical orbital elements in order: [a, e, i, ω, Ω, ν]. 
See Algorithm 10 from Vallado, D. "Fundamentals of Astrodynamics and Applications." 
4th Edition. 2013.

See also [`rv2coe`](@ref)

"""
function coe2rv(coe, μ=1)
    a, e, i, ω, Ω, ν = coe

    # Special cases
    circular = abs(e) < 1.0e-16
    equatorial = (abs(i) < 1.0e-16) || (abs(i-π) < 1.0e-16)
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


"""
    rv2coe(rv, μ=1)

Convert Cartesian position and velocity to classical orbital elements.


Output classical orbital elements in order: [a, e, i, ω, Ω, ν]. 
See Algorithm 9 from Vallado, D. "Fundamentals of Astrodynamics and Applications." 
4th Edition. 2013.

See also [`coe2rv`](@ref)

"""
function rv2coe(rv, μ=1)

    xvec = @SVector([1,0,0])
    zvec = @SVector([0,0,1])

    x, y, z, u, v, w = rv

    rvec = @SVector([x, y, z])
    vvec = @SVector([u, v, w])
    hvec = cross3(rvec, vvec)
    nvec = cross3(zvec, hvec)

    r = sqrt(sum(rvec.^2))
    v = sqrt(sum(vvec.^2))
    h = sqrt(sum(hvec.^2))
    n = sqrt(sum(nvec.^2))

    rdotv = sum(rvec.*vvec)
    evec = ((v^2-μ/r)*rvec - rdotv*vvec)/μ

    e = sqrt(sum(evec.^2))

    ξ = 0.5*v^2 - μ/r

    if (abs(e - 1.0) > 1.0e-16)
        a = -0.5*μ/ξ
        # p = a*(1-e^2)
    else
        # p = h^2/μ
        a = Inf
    end

    i = angle_vectors(hvec, zvec)

    Ω = angle_vectors(nvec, xvec)
    if nvec[2] < 0.0; Ω = 2.0*π - Ω; end

    ω = angle_vectors(evec, nvec)
    if evec[3] < 0.0; ω = 2.0*π - ω; end

    ν = angle_vectors(rvec, evec)
    if rdotv < 0.0; ν = 2.0*π - ν; end

    # Special cases
    circular = abs(e) < 1.0e-16
    equatorial = (abs(i) < 1.0e-16) || (abs(i-π) < 1.0e-16)
    if !circular && equatorial # ν is now true longitude of periapsis
        ω = 0.0
        ν = angle_vectors(rvec, evec)
        if evec[2]*cos(i) < 0.0; ν = 2.0*π - ν; end
    elseif circular && !equatorial # ν is now argument of latitude
        Ω = 0.0
        ν = angle_vectors(rvec, nvec)
        if rvec[3] < 0.0; ν = 2.0*π - ν; end
    elseif circular && equatorial # ν is now true longitude
        ω = 0.0
        Ω = 0.0
        ν = angle_vectors(rvec, xvec)
        if rvec[2]*cos(i) < 0.0; ν = 2.0*π - ν; end
    end

    coe = @SVector([a, e, i, ω, Ω, ν])

    return coe
end


"""
    cross3(a, b)

Compute cross product of 3-element vectors a and b.
"""
function cross3(a, b)
    return @SVector([a[2]*b[3]-a[3]*b[2], 
                     a[3]*b[1]-a[1]*b[3], 
                     a[1]*b[2]-a[2]*b[1]])
end


"""
    rot1(θ)

Compute matrix for rotation of θ about 1-axis.
"""
function rot1(θ)
    return @SMatrix(
    [1       0      0;
     0  cos(θ) sin(θ);
     0 -sin(θ) cos(θ)])
end


"""
    rot3(θ)

Compute matrix for rotation of θ about 3-axis.
"""
function rot3(θ)
    return @SMatrix(
    [ cos(θ)  sin(θ)  0;
     -sin(θ)  cos(θ)  0;
          0       0   1])
end


"""
    angle_vectors(a, b)

Compute angle between two vectors with high precision.
Credit to Alseidon for suggesting this as an alternative to using acos.
"""
function angle_vectors(a, b)
    angle = atan(sqrt(sum(cross3(a, b).^2)), sum(a.*b))
    return mod(angle, 2π)
end


export coe2rv
export rv2coe


end
