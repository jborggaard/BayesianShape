using Gmsh:gmsh
using LinearAlgebra
using Printf
using SparseArrays
using Arpack

using FEMfunctions

include("makeTriangleMesh.jl")
include("computeFEMmatrices.jl")

a = zeros(Float64,2)
a[1] = 0.635#0.5
a[2] = 0.275#1.0
x,eConn,boundaryNodes = makeTriangleMesh(a)

nBoundary = length(boundaryNodes)
dBoundary = zeros(Float64,nBoundary,1)
A,M = computeFEMmatrices(x,eConn,boundaryNodes,dBoundary)

nNodes = size(x,1)
ef1 = zeros(Float64,nNodes,1)

nev = 30
λ, ϕ = eigs(A,M; which=:SM, nev)

for i=1:nev
  tmp = real(λ[i])
  @printf(" %d:  %15.10g \n", i, tmp)
end


plotEigenfunctions = true
if ( plotEigenfunctions )
  interiorNodes = convert(Array{Int64,1},1:nNodes)
  setdiff!(interiorNodes,boundaryNodes)

  ef1 = zeros(Float64,nNodes,5)
  ef1[interiorNodes,1:5] = ϕ[:,1:5]
  saveFEMasVTK("triangle",x,eConn,["ef1","ef2","ef3","ef4","ef5"],ef1,[],[])
end

