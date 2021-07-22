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

N = 40  # number of BSplines used to represent the drum shape

validation = true
if ( validation )
  a0 = 1.0
  a = zeros(Float64,2)
  b = zeros(Float64,2)
else
  a0,a,b = sampleInnerGeometry()
end

r,err = fitBSpline2Fourier(a0,a,b,N)
@printf("B-spline approximation error (%d B-splines) is: %12.8f\n",N,err)

x,eConn,boundaryNodes = makeDrumMesh(r)

nBoundary = length(boundaryNodes)
dBoundary = zeros(Float64,nBoundary,1)
A,M = computeFEMmatrices(x,eConn,boundaryNodes,dBoundary)

nNodes = size(x,1)
ef1 = zeros(Float64,nNodes,1)

if validation
  nev = 13
  λ, ϕ = eigs(A,M; which=:SM, nev)

  @printf("  Numerical Eigenvalue    Analytical Solution\n")
  @printf("---------------------------------------------\n")
  tmp = real(λ[1])
  @printf(" 1:  %g     %g\n", tmp,besselj_zero(0,1)^2)

  tmp = real(λ[2])
  @printf(" 2:  %g     %g\n", tmp,besselj_zero(1,1)^2)

  tmp = real(λ[3])
  @printf(" 3:  %g     %g\n", tmp,besselj_zero(1,1)^2)

  tmp = real(λ[4])
  @printf(" 4:  %g     %g\n", tmp,besselj_zero(2,1)^2)

  tmp = real(λ[5])
  @printf(" 5:  %g     %g\n", tmp,besselj_zero(2,1)^2)

  tmp = real(λ[6])
  @printf(" 6:  %g     %g\n", tmp,besselj_zero(0,2)^2)

  tmp = real(λ[7])
  @printf(" 7:  %g     %g\n", tmp,besselj_zero(3,1)^2)

  tmp = real(λ[8])
  @printf(" 8:  %g     %g\n", tmp,besselj_zero(3,1)^2)
  
  tmp = real(λ[9])
  @printf(" 9:  %g     %g\n", tmp,besselj_zero(1,2)^2)
  
  tmp = real(λ[10])
  @printf("10:  %g     %g\n", tmp,besselj_zero(1,2)^2)
  
  tmp = real(λ[11])
  @printf("11:  %g     %g\n", tmp,besselj_zero(4,1)^2)

  tmp = real(λ[12])
  @printf("11:  %g     %g\n", tmp,besselj_zero(4,1)^2)

  tmp = real(λ[13])
  @printf("11:  %g     %g\n", tmp,besselj_zero(2,2)^2)

end


plotEigenfunctions = false
if ( plotEigenfunctions )
  λ, ϕ = eigs(A,M; which=:LM, nev=10)  # if we want to plot some eigenfunctions

  interiorNodes = convert(Array{Int64,1},1:nNodes)
  setdiff!(interiorNodes,boundaryNodes)

  ef1 = zeros(Float64,nNodes,5)
  ef1[interiorNodes,1:5] = ϕ[:,1:5]
  saveFEMasVTK("drum",x,eConn,["ef1","ef2","ef3","ef4","ef5"],ef1,[],[])
end

