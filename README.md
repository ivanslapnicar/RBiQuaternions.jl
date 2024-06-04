# RBiQuaternions
Julia package of reduced biquaternions. 

[![Build Status](https://github.com/ivanslapnicar/RBiQuaternions.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ivanslapnicar/RBiQuaternions.jl/actions/workflows/CI.yml?query=branch%3Amain)

A reduced biquaternion $a \in\mathbb{Q}_\mathbb{R}$ can be uniquely expressed as

$$
a = a_0 + a_1 \mathbf{i} + a_2 \mathbf{j} + a_3 \mathbf{k},$$

where $a_i \in \mathbb{R}$ for
$i = 0, 1, 2, 3$, and

$$
\mathbf{i}^2 = \mathbf{k}^2 = −1,\quad \mathbf{j}^2 = 1,\quad 
\mathbf{i}\mathbf{j} = \mathbf{j}\mathbf{i} = \mathbf{k},\quad  
\mathbf{j}\mathbf{k} = \mathbf{k}\mathbf{j} = \mathbf{i},\quad 
\mathbf{k}\mathbf{i} = \mathbf{i}\mathbf{k} = −\mathbf{j}.$$

The absolute value of $a$ is $|a| = \sqrt{a^2_0 + a^2_1+a^2_2+a^2_3}$

The conjugate of $a$ is $\bar{a}=a_0 - a_1 \mathbf{i} + a_2 \mathbf{j} - a_3 \mathbf{k}$.

Two basic singular (non-invertible) reduced biquaternions are

$$
e_1=\frac{1}{2}+\frac{1}{2}\mathbf{j}\quad \textrm{and} \quad e_2=\frac{1}{2}+\frac{1}{2}\mathbf{j}.$$

## Instalation

```jldoctest
pkg> add https://github.com/ivanslapnicar/RBiQuaternions.jl
```
or, after registering,
```jldoctest
pkg> add RBiQuaternions
```

## Examples

Basic properties
```jldoctest
julia> i=rbiquat(0,1,0,0)
RBiQuaternion{Int64}(0, 1, 0, 0)

julia> j=rbiquat(0,0,1,0)
RBiQuaternion{Int64}(0, 0, 1, 0)

julia> k=rbiquat(0,0,0,1)
RBiQuaternion{Int64}(0, 0, 0, 1)

julia> i^2 == -j^2 == k^2 == i*j*k == -1
true
```
Multiplication is commutative

```jldoctest
julia> a=rbiquat(1,2,3,4)
RBiQuaternion{Int64}(1, 2, 3, 4)

julia> b=rbiquat(5,6,7,8)
RBiQuaternion{Int64}(5, 6, 7, 8)

julia> a*b
RBiQuaternion{Int64}(-18, 68, -18, 60)

julia> b*a
RBiQuaternion{Int64}(-18, 68, -18, 60)
```

Sqare root and absolute value are well defined
```jldoctest
julia> c=√a
RBiQuaternionF64(1.5055993983149454, -0.14333523785320212, 0.8620051454093629, 1.4104387361768351)

julia> abs(c*c-a)
1.0877919644084146e-15
```
Basic singular quaternions are idempotent and mutually orthogonal

```jldoctest
julia> e₁
RBiQuaternionF64(0.5, 0.0, 0.5, 0.0)

julia> e₂
RBiQuaternionF64(0.5, 0.0, -0.5, 0.0)

julia> √e₁==e₁==e₁^2, e₂^2==e₂, e₁*e₂==0
(true, true, true)
```
Splitting into two complex parts
```
julia> d=splitc(a)
RBiQuaternions.SplitC{Complex{Int64}}(4 + 6im, -2 - 2im)

julia> d.c1*e₁+d.c2*e₂==a
true
```