using Gmsh:gmsh
using LinearAlgebra
using Printf
using SparseArrays
using Arpack
using FunctionZeros   # for the validation example

using FEMfunctions

include("makeDrumMesh.jl")
include("fitBSpline2Fourier.jl")
include("computeFEMmatrices.jl")
include("inputOutput.jl")


#  Here you generate a sample to call the inputOutput function
validation = true
if ( validation )
  a0 = 1.0
  a = zeros(Float64,2)
  b = zeros(Float64,2)
end

λ = inputOutput(a0,a,b)
