"""
    RBiQuaternion{T<:Real} <: Number

Reduced BiQuaternion number type with real and imaginary parts of type `T`.

`RBiQuaternionF16`, `RBiQuaternionF32`, and `RBiQuaternionF64` are aliases for
`RBiQuaternion{Float16}`, `RBiQuaternion{Float32}`, and `RBiQuaternion{Float64}`, respectively.

See also: [`rbiquat`](@ref), [`real`](@ref), [`imag`](@ref).
"""
struct RBiQuaternion{T<:Real} <: Number
	s::T
	v1::T
	v2::T
	v3::T
end

using Random 

const RBiQuaternionF16 = RBiQuaternion{Float16}
const RBiQuaternionF32 = RBiQuaternion{Float32}
const RBiQuaternionF64 = RBiQuaternion{Float64}

RBiQuaternion{T}(x::Real) where {T<:Real} = RBiQuaternion(convert(T, x))
RBiQuaternion{T}(q::RBiQuaternion) where {T<:Real} = RBiQuaternion{T}(q.s, q.v1, q.v2, q.v3)
RBiQuaternion(s::Real, v1::Real, v2::Real, v3::Real) = RBiQuaternion(promote(s, v1, v2, v3)...)
RBiQuaternion(x::Real) = RBiQuaternion(x, zero(x), zero(x), zero(x))

Base.promote_rule(::Type{RBiQuaternion{T}}, ::Type{S}) where {T <: Real, S <: Real} = RBiQuaternion{promote_type(T, S)}
Base.promote_rule(::Type{RBiQuaternion{T}}, ::Type{RBiQuaternion{S}}) where {T <: Real, S <: Real} = RBiQuaternion{promote_type(T, S)}

# Added
Base.isless(x::RBiQuaternion{T}, y::RBiQuaternion{T}) where T=Base.isless(abs(x),abs(y))
LinearAlgebra.norm(x::Dual{T}) where T = Dual(norm(realpart(x)),norm(dualpart(x)))
"""
    rbiquat(w, [x, y, z])

Convert real numbers or arrays to reduced biquaternion quaternion. `x, y, z` defaults to zero.

# Examples
```jldoctest
julia> rbiquat(7)
RBiQuaternion{Int64}(7, 0, 0, 0)

julia> rbiquat(1.0, 2, 3, 4)
RBiQuaternionF64(1.0, 2.0, 3.0, 4.0)

julia> rbiquat([1, 2, 3])
3-element Vector{Quaternion{Int64}}:
 RBiQuaternion{Int64}(1, 0, 0, 0)
 RBiQuaternion{Int64}(2, 0, 0, 0)
 RBiQuaternion{Int64}(3, 0, 0, 0)
```
"""
rbiquat

rbiquat(q::RBiQuaternion) = q
rbiquat(s::Real) = RBiQuaternion(s)
rbiquat(s::Real, v1::Real, v2::Real, v3::Real) = RBiQuaternion(s, v1, v2, v3)

## Array operations on quaternions ##
rbiquat(A::AbstractArray{<:RBiQuaternion}) = A
function rbiquat(A::AbstractArray{T}) where T
    if !isconcretetype(T)
        error("`quat` not defined on abstractly-typed arrays; please convert to a more specific type")
    end
    convert(AbstractArray{typeof(rbiquat(zero(T)))}, A)
end

rbiquat(p, v1, v2, v3) = RBiQuaternion(p, v1, v2, v3)



"""
    rbiquat(T::Type)

Return an appropriate type that can represent a value of type `T` as a reduced biquaternion.
Equivalent to `typeof(rbiquat(zero(T)))`.

# Examples
```jldoctest
julia> rbiquat(RBiQuaternion{Int})
RBiQuaternion{Int64}

julia> rbiquat(Int)
RBiQuaternion{Int64}
```
"""
rbiquat(::Type{T}) where {T<:Real} = RBiQuaternion{T}
rbiquat(::Type{RBiQuaternion{T}}) where {T<:Real} = RBiQuaternion{T}

rbiquat(::Missing) = missing
# Same definitioin as in Base: https://github.com/JuliaLang/julia/blob/v1.9.3/base/missing.jl#L111-L115
rbiquat(::Type{Missing}) = Missing
function rbiquat(::Type{Union{T, Missing}}) where T
    T === Any && throw(MethodError(rbiquat, (Any,)))  # To prevent StackOverflowError
    Union{rbiquat(T), Missing}
end

"""
    real(T::Type{<:Quaternion})

Return the type that represents the real part of a value of type `T`.
e.g: for `T == RBiQuaternion{R}`, returns `RBiQuaternion(R)`.
Equivalent to `typeof(real(zero(T)))`.

# Examples
```jldoctest
julia> real(RBiQuaternion{Int})
RBiQuaternion{Int64}
```
"""
Base.real(::Type{RBiQuaternion{T}}) where {T} = RBiQuaternion{T}

"""
    real(q::RBiQuaternion)

Return the real part of the quaternion `q` defined as `real(q)=RBiQuaternion(q.s,zero(q.v1),q.v2,zero(q.v3))=(q+conj(q))/2)`

Note that this function is different from `Base.imag`, which returns `Real` for complex numbers and quaternions.

See also: [`imag`](@ref), [`rbiquat`](@ref)

# Examples
```jldoctest
julia> real(rbiquat(1,2,3,4))
RBiQuaternion(1,0,3,4)
```
"""
Base.real(q::RBiQuaternion) = RBiQuaternion(q.s,0,q.v2,0)
# Base.real(::AbstractArray{<:RBiQuaternion})

# For testing
Base.eps(::RBiQuaternion{T}) where T=eps(T)
Base.rtoldefault(::Type{RBiQuaternion{T}}) where T=Base.rtoldefault(T)

"""
    imag(q::RBiQuaternion{T})

Return the imaginary part of the reduced biquaternion `q` defined as `imag(q)=RBiQuaternion(zero(q.s),q.v1,zero(q.v2),q.v3)=(q-conj(q))/2)`

Note that this function is different from `Base.imag`, which returns `Real` for complex numbers and a triple for quaternions.

See also: [`real`](@ref), [`conj`](@ref).

# Examples
```jldoctest
julia> imag(RBiQuaternion(1,2,3,4))
RBiQuaternion(0,2, 0, 4)
```
"""
Base.imag(q::RBiQuaternion)=RBiQuaternion(0,q.v1,0,q.v3)

# imag_part(q::RBiQuaternion) = (q.v1, q.v2, q.v3)

Base.:/(q::RBiQuaternion, x::Real) = RBiQuaternion(q.s / x, q.v1 / x, q.v2 / x, q.v3 / x)
Base.:*(q::RBiQuaternion, x::Real) = RBiQuaternion(q.s * x, q.v1 * x, q.v2 * x, q.v3 * x)
Base.:*(x::Real, q::RBiQuaternion) = q * x
Base.:*(q::RBiQuaternion, x::Complex) = q*RBiQuaternion(real(x),imag(x),0.0,0.0)
Base.:*(x::Complex, q::RBiQuaternion) = q*x

e₁=RBiQuaternion(0.5,0.0,0.5,0.0)
e₂=RBiQuaternion(0.5,0.0,-0.5,0.0)

"""
    conj(q::RBiQuaternion)

Compute the reduced biquaternion conjugate of a reduced biquaternion `q`.

# Examples
```jldoctest
julia> conj(RBiQuaternion(1,2,3,4))
RBiQuaternion{Int64}(1, -2, 3, -4)
```
"""
Base.conj(q::RBiQuaternion) = RBiQuaternion(q.s, -q.v1, q.v2, -q.v3)

function Base.:*(q::RBiQuaternion, w::RBiQuaternion)
    s  = (q.s * w.s + q.v2 * w.v2) - (q.v1 * w.v1 + q.v3 * w.v3)
    v1 = (q.s * w.v1 + q.v1 * w.s) + (q.v2 * w.v3 + q.v3 * w.v2)
    v2 = (q.s * w.v2 + q.v2 * w.s) - (q.v3 * w.v1 + q.v1 * w.v3)
    v3 = (q.s * w.v3 + q.v3 * w.s) + (q.v1 * w.v2 + q.v2 * w.v1)
    return RBiQuaternion(s, v1, v2, v3)
end

function Base.abs(q::RBiQuaternion)
    a = max(abs(q.s), abs(q.v1), abs(q.v2), abs(q.v3))
    if isnan(a) && isinf(q)
        return typeof(a)(Inf)
    elseif iszero(a) || isinf(a)
        return a
    else
        # return sqrt(abs2(q / a)) * a
        return sqrt((q.s/a)^2+(q.v1/a)^2+(q.v2/a)^2+(q.v3/a)^2)*a
    end
end

Base.float(q::RBiQuaternion{T}) where T = convert(RBiQuaternion{float(T)}, q)

function abs_imag(q::RBiQuaternion)
    a = max(abs(q.v1), abs(q.v2), abs(q.v3))
    if isnan(a) && (isinf(q.v1) | isinf(q.v2) | isinf(q.v3))
        return oftype(a, Inf)
    elseif iszero(a) || isinf(a)
        return a
    else
        return sqrt((q.v1 / a)^2 + (q.v2 / a)^2 + (q.v3 / a)^2) * a
    end
end

Base.abs2(q::RBiQuaternion) = abs(q)^2

"""
    inv(q::RBiQuaternion)

Return the multiplicative inverse of `q::RBiQuaternion`, such that `q*inv(q)` or `inv(q)*q`
yields `one(q)` (the multiplicative identity) up to roundoff errors.
It is computed using splitting into pair of complex numbers.

# Examples
```jldoctest
julia> inv(rbiquat(1))
RBiQuaternionF64(1.0, -0.0, -0.0, -0.0)

julia> inv(rbiquat(1, 2, 0, 0))
It is computed using splitting into pair of complex numbers.
RBiQuaternionF64(0.2, -0.4, 0.0, 0.0)

julia> inv(rbiquat(2//3))
RBiQuaternionF64(1.5, 0.0, 0.0, 0.0)
```
"""
Base.inv(q::RBiQuaternion)=inv(splitc(q).c1)*e₁+inv(splitc(q).c2)*e₂

"""
    sqrt(q::RBiQuaternion)

Return the square root `r` of `q::RBiQuaternion`, such that `r*r=q` up to roundoff errors.
It is computed using splitting into pair of complex numbers.

# Examples
```jldoctest
julia> sqrt(rbiquat(1))
RBiQuaternionF64(1.0, 0.0, 0.0, 0.0)

julia> sqrt(rbiquat(1, 2, 0, 0))
RBiQuaternionF64(1.272019649514069, 0.7861513777574233, 0.0, 0.0)
```
"""
Base.sqrt(q::RBiQuaternion)=sqrt(splitc(q).c1)*e₁+sqrt(splitc(q).c2)*e₂

# Base.isreal(q::RBiQuaternion) = iszero(q.v1) & iszero(q.v2) & iszero(q.v3)
Base.isreal(q::RBiQuaternion) = iszero(q.v1) & iszero(q.v3)
Base.isfinite(q::RBiQuaternion) = isfinite(q.s) & isfinite(q.v1) & isfinite(q.v2) & isfinite(q.v3)
Base.iszero(q::RBiQuaternion) = iszero(q.s) & iszero(q.v1) & iszero(q.v2) & iszero(q.v3)
Base.isnan(q::RBiQuaternion) = isnan(q.s) | isnan(q.v1) | isnan(q.v2) | isnan(q.v3)
Base.isinf(q::RBiQuaternion) = isinf(q.s) | isinf(q.v1) | isinf(q.v2) | isinf(q.v3)
Base.isinteger(q::RBiQuaternion) = isinteger(q.s) & integer(q.v1) & isinteger(q.v2) & isinteger(q.v3)

# included strictly for documentation; the base implementation is sufficient
"""
    sign(q::RBiQuaternion) -> RBiQuaternion

Return zero if `q==0` and ``q/|q|`` otherwise.

# Examples
```jldoctest
julia> sign(RBiQuaternion(4, 0, 0, 0))
QuaternionF64(1.0, 0.0, 0.0, 0.0)

julia> sign(RBiQuaternion(1, 0, 1, 0))
RBiQuaternionF64(0.7071067811865475, 0.0, 0.7071067811865475, 0.0)
```
"""
sign(::RBiQuaternion)

# Base.sign(q::RBiQuaternion)=sign(q.s)

Base.:-(q::RBiQuaternion) = RBiQuaternion(-q.s, -q.v1, -q.v2, -q.v3)

Base.:+(q::RBiQuaternion, w::RBiQuaternion) =
    RBiQuaternion(q.s + w.s, q.v1 + w.v1, q.v2 + w.v2, q.v3 + w.v3)

Base.:-(q::RBiQuaternion, w::RBiQuaternion) =
    RBiQuaternion(q.s - w.s, q.v1 - w.v1, q.v2 - w.v2, q.v3 - w.v3)

Base.:/(q::RBiQuaternion{T}, w::RBiQuaternion{T}) where T = q*inv(w)
Base.:\(q::RBiQuaternion{T}, w::RBiQuaternion{T}) where T = inv(q)*w
Base.://(x::RBiQuaternion, y::Real) = biquat(x.s//y, x.v1//y, x.v2//y, x.v3//y)
Base.://(x::Number, y::RBiQuaternion) = x*conj(y)//abs2(y)

Base.:(==)(q::RBiQuaternion, w::RBiQuaternion) = (q.s == w.s) & (q.v1 == w.v1) & (q.v2 == w.v2) & (q.v3 == w.v3)
function Base.isequal(q::RBiQuaternion, w::RBiQuaternion)
    isequal(q.s, w.s) & isequal(q.v1, w.v1) & isequal(q.v2, w.v2) & isequal(q.v3, w.v3)
end

Base.widen(::Type{RBiQuaternion{T}}) where {T} = RBiQuaternion{widen(T)}

Base.flipsign(x::RBiQuaternion, y::Real) = ifelse(signbit(y), -x, x)

function Base.read(io::IO, ::Type{RBiQuaternion{T}}) where T<:Real
    return RBiQuaternion{T}(ntuple(_ -> read(io, T), Val(4))...)
end
Base.write(io::IO, q::RBiQuaternion) = write(io, q.s, q.v1, q.v2, q.v3...)

Base.big(::Type{RBiQuaternion{T}}) where {T<:Real} = RBiQuaternion{big(T)}
Base.big(z::RBiQuaternion{T}) where {T<:Real} = RBiQuaternion{big(T)}(z)

function Base.rand(rng::AbstractRNG, ::Random.SamplerType{RBiQuaternion{T}}) where {T<:Real}
    RBiQuaternion{T}(rand(rng, T), rand(rng, T), rand(rng, T), rand(rng, T))
end

function Base.randn(rng::AbstractRNG, ::Type{RBiQuaternion{T}}) where {T<:AbstractFloat}
    RBiQuaternion{T}(
        randn(rng, T) / 2,
        randn(rng, T) / 2,
        randn(rng, T) / 2,
        randn(rng, T) / 2,
    )
end

"""
    round(q::RBiQuaternion, RoundingModeReal,
          RoundingModeImaginary1, RoundingModeImaginary2, RoundingModeImaginary3; kwargs...)

Return the nearest integral value of the same type as the reduced biquaternion-valued `q` to `q`,
breaking ties using the specified `RoundingMode`s.

The first `RoundingMode` is used for rounding the real part while the second is used
for rounding the imaginary parts. Alternatively, a `RoundingMode` may be provided for each
part.

The `kwargs` are the same as those for `round(::Real[, RoundingMode]; kwargs...)`.

# Example
```jldoctest
julia> round(rbiquat(3.14, 4.5, 8.3, -2.8))
RBiQuaternionF64(3.0, 4.0, 8.0, -3.0)
```
"""
function Base.round(
    q::RBiQuaternion,
    rs::RoundingMode=RoundNearest,
    rv::RoundingMode=rs;
    kwargs...,
)
    return round(q, rs, rv, rv, rv; kwargs...)
end

function Base.round(
    q::RBiQuaternion,
    rs::RoundingMode,
    rv1::RoundingMode,
    rv2::RoundingMode,
    rv3::RoundingMode;
    kwargs...,
)
    return RBiQuaternion(
        round(q.s, rs; kwargs...),
        round(q.v1, rv1; kwargs...),
        round(q.v2, rv2; kwargs...),
        round(q.v3, rv3; kwargs...),
    )
end

"""
Splittings into real and complex parts.

	splitr(q::RBiQuaternion) -> SplitR(q.s, q.v1, q.v2, q.v3)
	splitc(q::RBiQuaternion{T}) -> SplitC{Complex{T}}(c1,c2)

splitr() returns the spltting of reduced biquaternion q into 4-tuple SplitR. Also works for AbstractArrays.

splitc() returns spltting of reduced biquaternion q into complex tuple SplitC(c1,c2) such that 
q=c1*e₁+c2*e₂.
Also works for AbstractArrays.

# Examples
```jldoctest
julia> splitr(rbiquat(1,2,3,4))
SplitR{Int64}(1, 2, 3, 4)

julia> splitc(rbiquat(1,2,3,4))
SplitC{Complex{Int64}}(4 + 6im, -2 - 2im)

splitc(rbiquat([1 2;3 4]))
SplitC{Matrix{Complex{Int64}}}(Complex{Int64}[1 + 0im 2 + 0im; 3 + 0im 4 + 0im], Complex{Int64}[1 + 0im 2 + 0im; 3 + 0im 4 + 0im])
```
"""
struct SplitR{T} 
	s::T
	v1::T
	v2::T
	v3::T
end

struct SplitC{T}
	c1::T
	c2::T
end

function splitr(A::AbstractArray{RBiQuaternion{T}}) where T
	if ndims(A)==1
		n=length(A)
		As=[A[i].s for i=1:n]
		Av1=[A[i].v1 for i=1:n]
		Av2=[A[i].v2 for i=1:n]
		Av3=[A[i].v3 for i=1:n]
	else
		m,n=size(A)
    	As=[A[i,j].s for i=1:m, j=1:n]
		Av1=[A[i,j].v1 for i=1:m, j=1:n]
		Av2=[A[i,j].v2 for i=1:m, j=1:n]
		Av3=[A[i,j].v3 for i=1:m, j=1:n]
	end
	return SplitR(As,Av1,Av2,Av3)
end

function splitc(A::AbstractArray{RBiQuaternion{T}}) where T
	Q=splitr(A)
	Ac1=Q.s+Q.v2+im*(Q.v1+Q.v3)
	Ac2=Q.s-Q.v2+im*(Q.v1-Q.v3)
	return SplitC(Ac1,Ac2)
end

function splitc(a::RBiQuaternion{T}) where T
	Ac1=a.s+a.v2+im*(a.v1+a.v3)
	Ac2=a.s-a.v2+im*(a.v1-a.v3)
	return SplitC(Ac1,Ac2)
end

function splitr(a::RBiQuaternion{T}) where T
	return SplitR(a.s, a.v1, a.v2, a.v3)
end

