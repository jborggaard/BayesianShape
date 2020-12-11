## Given a Fourier parametrization of the inner boundary, compute the solutions to the Stokes and Advection-Diffusion problems
#  Inputs:
#    a        coefficients of cosines
#    b        coefficients of sines
#    a0       mean inner boundary
#    N        number of B splines used to approximate Fourier representation
#    ω        angular velocity
#    verbose  verbosity
#
#  Outputs:
#    sa       structure describing the mesh and the velocity, pressure, and temperature fields

#using Gmsh:gmsh
#using LinearAlgebra
#using Makie
#using AbstractPlotting
#using SparseArrays
#using SpecialMatrices
##using Plots
#using Polynomials
#using Printf
#using Random
#using WriteVTK
#using Distributions
#
#include("makeMesh.jl")
#include("fitBSpline2Fourier.jl")
#include("saveFEMasVTK.jl")
#include("sampleInnerGeometry.jl")
#include("twodQuadratureRule.jl")
#include("twodShape.jl")
#include("twodBilinear.jl")
#include("twodLinForm.jl")
#include("twodStokesRotating.jl")
#include("twodAdvectionDiffusion.jl")
#include("twodProjectDerivatives.jl")
#include("computeC.jl")
#include("computeVorticity.jl")
#include("solutionArray.jl")

function twodStokesAD(a,b,a0,N; ω = 10.0, κ = 1.0, verbose=true)
  ### Generate the finite element mesh using Gmsh (implemented in makeMesh)
  
  r,err = fitBSpline2Fourier(a0,a,b,N)
  #r = ones(N,1); # uncomment for some degubbing
  verbose && @printf("B-spline approximation error (%d B-splines) is: %12.8f\n",N,err);
  
  ### Then map them to a distribution between 0.5 and 1.5 using the arctan function, small alpha values (e.g. 0.1) cluster the results of the B-spline parameters around 1
  α = 1.0;
  r = 1.0 .+ atan.(α*r)/π;
  
  x,eConn,eConn2, innerNodes,innerX, outerNodes,outerX = makeMesh(r);
  
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
  
  #compute Stokes flow
  velocity, pressure = twodStokesRotatingOuter(xT,eC,innerNodes,outerNodes,ω)
  
  #solve steady Advection-Diffusion equation
  temperature, massMat = twodAdvectionDiffusion(xT,eC,innerNodes,outerNodes,velocity,κ)
  
  #squish everything we might need into a structure
  sa = solutionArray();
  sa.xT  = xT;
  sa.eC  = eC;
  sa.eC2 = eC2;
  sa.velocity = velocity;
  sa.pressure = pressure;
  sa.temperature = temperature;
  sa.massMat     = massMat;

  return sa;
end
