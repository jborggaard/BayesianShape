## Creates the input/output map for an optimization problem

#* 1.) Sample coefficients of a sinusoidal description of the inner boundary
#* 2.) Find coefficients for an approximating B-Spline
#* 3.) Build a finite element mesh using Gmsh routines
#* 4.) Approximate Stokes equations using Taylor-Hood elements and a penalty method
#* 5.) Approximate the advection-diffusion equations using quadratic finite elements
#* 6.) Compute outputs

using Gmsh:gmsh
using LinearAlgebra
using Makie
using AbstractPlotting
using SparseArrays
using SpecialMatrices
#using Plots
using Polynomials
using Printf
using Random
using WriteVTK

include("makeMesh.jl")
include("fitBSpline2Fourier.jl")
include("saveFEMasVTK.jl")
include("sampleInnerGeometry.jl")
include("twodQuadratureRule.jl")
include("twodShape.jl")
include("twodBilinear.jl")
include("twodLinForm.jl")
include("twodStokesRotating.jl")
include("twodAdvectionDiffusion.jl")
include("twodProjectDerivatives.jl")
include("computeC.jl")
include("computeVorticity.jl")

###  Define parameters for the simulation
ω = 10.0;    # rotational velocity

a0,a,b = sampleInnerGeometry()

### define the 40 parameters that describe the inner boundary
N = 40;      # number of BSplines used to represent the inner boundary

#### First draw 40 parameters from a normal distribution
#Random.seed!(1);                        # for sane reproducibility in debugging
#param = randn(Float64,N);
#param = ones(Float64,N);
#param = zeros(N);
#param[1] = 0; param[2] = 1;

#### Then map them to a distribution between 0.5 and 1.5 using the arctan function, small alpha values (e.g. 0.1) cluster the results of the B-spline parameters around 1

#α = 1.0;
#r = 1.0 .+ atan.(α*param)/π;
### Generate the finite element mesh using Gmsh (implemented in makeMesh)

r = fitBSpline2Fourier(a0,a,b,N)
#r = ones(N,1); # uncomment for some degubbing
x,eConn,eConn2, innerNodes,innerX, outerNodes,outerX = makeMesh(r)

sort!(innerNodes);
sort!(outerNodes);


#   Let's look at the mesh...
nNodes = size(x,2)
nElements = size(eConn,2)
nElementsSensor = size(eConn2,2)

xT = zeros(Float64,nNodes,2)   # xT = transpose(x[1:2,:])
eC = zeros(Int64,nElements,6)  # eC = transpose(eConn), converted to Int64
for i=1:nNodes
  xT[i,1] = x[1,i]
  xT[i,2] = x[2,i]
end
for i=1:nElements#-nElementsSensor # the sensor elements are orientated correctly (why?)
#  for j=1:6
#    eC[i,j] = convert(Int64,eConn[j,i])
#  end
   eC[i,1] = convert(Int64,eConn[1,i])
   eC[i,2] = convert(Int64,eConn[3,i])
   eC[i,3] = convert(Int64,eConn[2,i])
   eC[i,4] = convert(Int64,eConn[6,i])
   eC[i,5] = convert(Int64,eConn[5,i])
   eC[i,6] = convert(Int64,eConn[4,i])
end

eC2  = zeros(Int64,nElementsSensor,6)
for i=1:nElementsSensor
  for j=1:6
    eC2[i,j] = convert(Int64,eConn2[j,i])
  end
end  # elements are positively orientated
C    = computeC(xT,eC2)

Call = computeC(xT,eC)

#for i=nElements-nElementsSensor+1 : nElements
#  for j=1:6
#    eC[i,j] = convert(Int64,eConn[j,i])
#  end
#end
#Plots.plot([x[1,innerNodes],x[1,outerNodes]],[x[2,innerNodes],x[2,outerNodes]],seriestype = :scatter)


velocity, pressure = twodStokesRotating(xT,eC,innerNodes,outerNodes,ω)

# rotate the entire field the opposite direction (-ω)
for j=1:nNodes
  velocity[j,1] = velocity[j,1] + ω*xT[j,2]
  velocity[j,2] = velocity[j,2] - ω*xT[j,1]
end

temperature = twodAdvectionDiffusion(xT,eC,innerNodes,outerNodes,velocity)

#  Output the solution or visualize
scalarLabels = ["temperature"]
vectorLabels = ["velocity"]
saveFEMasVTK("mixing",xT,eC,scalarLabels,temperature,vectorLabels,velocity)

#poly(xT, eC[:,1:3], color = velocity[:,2], strokecolor = (:black, 0.6), strokewidth = .2)

vorticity = computeVorticity(xT,eC,velocity)
#saveFEMasVTK("mixing",xT,eC,scalarLabels,vorticity,vectorLabels,velocity)
poly(xT, eC[:,1:3], color = vorticity[:,1], strokecolor = (:black, 0.6), strokewidth = 0.2)

velMag = sqrt.( velocity[:,1].*velocity[:,1] + velocity[:,2].*velocity[:,2] )
#poly(xT, eC[:,1:3], color = velMag, strokecolor = (:black, 0.6), strokewidth = .3)

#poly(xT, eC[:,1:3], color = temperature[:,1], strokecolor = (:black, 0.6), strokewidth = .2)

# compute the average temperature, x- and y-components of velocity, and vorticity.

V2 = sum(C)
Volume = sum(Call)

res = (C*temperature)/V2
@printf("the average temperature over the subregion is %g\n",res[1])
res = (Call*temperature)/Volume
@printf("                    and over the entire region is %g\n\n",res[1])

res = (C*vorticity)/V2
@printf("the average vorticity over the subregion is %g\n",res[1])
res = (Call*vorticity)/Volume
@printf("                  and over the entire region is %g\n\n\n",res[1])
