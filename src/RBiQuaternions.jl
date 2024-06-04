module RBiQuaternions

using Random
using LinearAlgebra

include("RBiQuaternion.jl")

export RBiQuaternion
export RBiQuaternionF16, RBiQuaternionF32, RBiQuaternionF64
export rbiquat
export splitr
export splitc
export e₁
export e₂

end # module
