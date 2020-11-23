using Gmsh:gmsh
using LinearAlgebra
#using Makie
#using AbstractPlotting
using SparseArrays
using SpecialMatrices
#using Plots
using Polynomials
using Printf
using Random
#using WriteVTK
using Distributions
#using HDF5

include("makeMesh.jl")
include("fitBSpline2Fourier.jl")
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

### Number of samples to draw
nSamples = 10;

### Dimension of unknowns
unkDims = [5; 10; 20; 40; 80; 160; 320];

# number of BSplines used to represent the inner boundary
Ns = [40; 80; 100; 120; 160];      

#@printf("Running with Fourier dim=%d and B-spline dim=%d\n",unkDim,N);
errs = zeros(length(unkDims),length(Ns));

@printf("%4s: "," ");
for j=1:length(Ns)
  @printf("%12d ",Ns[j]);
end
for i=1:length(unkDims)
  unkDim=unkDims[i];
  #prior
  prior = MvNormal(zeros(unkDim),sqrt(sqrt(2)).^-(0:unkDim-1));

  @printf("\n%4d: ",unkDim);
  for j=1:length(Ns)
    N=Ns[j];
    for k=1:nSamples
      #draw samples
      a0 = 1.0;
      a  = rand(prior);
      b  = rand(prior);
    
      #approximate with B-splines
      r,err = fitBSpline2Fourier(a0,a,b,N)
      #@printf("B-spline approximation error is: %12.8f\n",err);
      
      errs[i,j] += err; #mean
      #errs[i,j] = max(err,errs[i,j]); #\infty norm
    end
    errs[i,j] /= nSamples; #mean
    @printf("%12.8f ",errs[i,j]);
  end
end
@printf("\n");

