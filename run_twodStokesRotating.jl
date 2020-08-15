## Creates the input/output map for an optimization problem

#* From samples of a sinusoidal description, find an approximating B-Spline
#* Build a finite element mesh using Gmsh functions
#* Approximate Stokes equations using Taylor-Hood elements and a penalty method
#* Approximate the advection-diffusion equations
#* Compute output variables

using Gmsh:gmsh
using LinearAlgebra
using Makie
using AbstractPlotting
using SparseArrays
#using Plots
using Printf
using Random
using WriteVTK

include("makeMesh.jl")
include("saveFEMasVTK.jl")
include("twodQuadratureRule.jl")
include("twodShape.jl")
include("twodBilinear.jl")
include("twodLinForm.jl")
include("twodStokesRotating.jl")
include("twodAdvectionDiffusion.jl")


ω = 0.0;    # rotational velocity

### define the 40 parameters that describe the inner boundary

#### First draw 40 parameters from a normal distribution
# for sane reproducibility in debugging
Random.seed!(1);
param = randn(Float64,40);
#param = ones(Float64,40);
#param = zeros(40);
#param[1] = 0; param[2] = 1;

#### Then map them to a distribution between 0.5 and 1.5 using the arctan function, small alpha values (e.g. 0.1) cluster the results of the B-spline parameters around 1

α = 1;
r = 1.0 .+ atan.(α*param)/π;

### Generate the finite element mesh using Gmsh (implemented in makeMesh)

x,eConn, innerNodes,innerX, outerNodes,outerX = makeMesh(r)

sort!(innerNodes);
sort!(outerNodes);


#   Let's look at the mesh...
nNodes = size(x,2)
nElements = size(eConn,2)

xT = zeros(Float64,nNodes,2)
eC = zeros(Int64,nElements,6)
for i=1:nNodes
  xT[i,1] = x[1,i]
  xT[i,2] = x[2,i]
end
for i=1:nElements
  for j=1:6
    eC[i,j] = convert(Int64,eConn[j,i])
  end
end

#Plots.plot([x[1,innerNodes],x[1,outerNodes]],[x[2,innerNodes],x[2,outerNodes]],seriestype = :scatter)


velocity, pressure = twodStokesRotating(xT,eC,innerNodes,outerNodes,ω)

temperature = twodAdvectionDiffusion(xT,eC,innerNodes,outerNodes,velocity)

scalarLabels = ["temperature"]
vectorLabels = ["velocity"]
saveFEMasVTK("mixing",xT,eC,scalarLabels,temperature,vectorLabels,velocity)

velMag = sqrt.( velocity[:,1].*velocity[:,1] + velocity[:,2].*velocity[:,2] )
#poly(xT, eC[:,1:3], color = velocity[:,1], strokecolor = (:black, 0.6), strokewidth = .2)
#poly(xT, eC[:,1:3], color = velMag, strokecolor = (:black, 0.6), strokewidth = .3)
poly(xT, eC[:,1:3], color = temperature[:,1], strokecolor = (:black, 0.6), strokewidth = .2)
