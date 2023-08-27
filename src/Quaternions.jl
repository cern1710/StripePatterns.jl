using LinearAlgebra

module QuaternionModule

export Quaternion, slerp

struct Quaternion
    s::Float64
    v::Vector{Float64}
end

function Quaternion(s::Float64=0.0)
    return Quaternion(s, [0.0, 0.0, 0.0])
end

function Quaternion(v::Vector{Float64})
    return Quaternion(0.0, v)
end

function Quaternion(z::Complex{Float64})
    return Quaternion(real(z), [imag(z), 0.0, 0.0])
end

# addition and subtraction
Base.:+(q1::Quaternion, q2::Quaternion) = Quaternion(q1.s + q2.s, q1.v .+ q2.v)
Base.:-(q1::Quaternion, q2::Quaternion) = Quaternion(q1.s - q2.s, q1.v .- q2.v)

# unary
Base.:-(q::Quaternion) = Quaternion(-q.s, -q.v)

# multiplication with scalar
Base.:*(q::Quaternion, c::Float64) = Quaternion(q.s * c, q.v .* c)
Base.:*(c::Float64, q::Quaternion) = q * c

# division by sclar
Base.:/(q::Quaternion, c::Float64) = Quaternion(q.s / c, q.v ./ c)

function *(q1::Quaternion, q2::Quaternion)
    s1 = q1.s
    s2 = q2.s
    v1 = q1.v
    v2 = q2.v
    return Quaternion(s1*s2 - dot(v1, v2), s1.*v2 .+ s2.*v1 .+ cross(v1, v2))
end

conj(q::Quaternion) = Quaternion(q.s, -q.v)
inv(q::Quaternion) = conj(q) / norm2(q)
norm2(q::Quaternion) = q.s^2 + dot(q.v, q.v)
normalize(q::Quaternion) = q / norm(q)

function slerp(q0::Quaternion, q1::Quaternion, t::Float64)
    m0 = norm(q0)
    m1 = norm(q1)
    m = (1.0-t)*m0 + t*m1
    
    p0 = q0 / m0
    p1 = q1 / m1
    theta = acos(real(conj(p0)*p1))
    p = (sin((1.0-t)*theta)*p0 + sin(t*theta)*p1) / sin(theta)
    
    return m * p
end

Base.show(io::IO, q::Quaternion) = print(io, "( ", q.s, ", ", q.v, " )")

end