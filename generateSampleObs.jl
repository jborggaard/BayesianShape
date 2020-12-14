## Creates the input/output map for an optimization problem

#* 1.) Sample coefficients of a sinusoidal description of the inner boundary
#* 2.) Find coefficients for an approximating B-Spline
#* 3.) Build a finite element mesh using Gmsh routines
#* 4.) Approximate Stokes equations using Taylor-Hood elements and a penalty method
#* 5.) Approximate the advection-diffusion equations using quadratic finite elements
#* 6.) Compute outputs

using Gmsh:gmsh
using LinearAlgebra
#using Makie
using CairoMakie
using AbstractPlotting
using SparseArrays
using SpecialMatrices
#using Plots
using Polynomials
using Printf
using Random
using WriteVTK
using Distributions
using HDF5

include("makeMesh.jl")
include("fitBSpline2Fourier.jl")
include("saveFEMasVTK.jl")
include("sampleInnerGeometry.jl")
include("twodQuadratureRule.jl")
include("twodShape.jl")
include("twodBilinear.jl")
include("twodLinForm.jl")
include("twodStokesRotatingOuter.jl")
include("twodAdvectionDiffusion.jl")
include("twodProjectDerivatives.jl")
include("computeC.jl")
include("computeVorticity.jl")

### Number of samples to draw
nSamples = 10;

### Dimension of unknowns
unkDim = 160;

### BSplines
N = 160;      # number of BSplines used to represent the inner boundary

outFile = "out_F$(unkDim)_B$(N).h5";


###  Define parameters for the simulation
ω = 10.0;    # rotational velocity
κ = 0.10;     # diffusion constant

#observations
#  1. avg temp (subregion)
#  2. avg temp (whole domain)
#  3. avg vort (subregion)
#  4. avg vort (whole domain)
#  5. scalar variance (whole domain)
obs = zeros(nSamples,5);


#prior
#prior = MvNormal(zeros(unkDim),sqrt(sqrt(2)).^-(0:unkDim-1));
prior = MvNormal(zeros(unkDim),2.0.^(-(0:unkDim-1)./4));

aArr = zeros(nSamples,unkDim);
bArr = zeros(nSamples,unkDim);

for i=1:nSamples
  istr = @sprintf("%04d",i);

  #a0,a,b = sampleInnerGeometry()
  a0 = 1.0;
  a  = rand(prior);
  b  = rand(prior);
  aArr[i,:] = a;
  bArr[i,:] = b;
  
  #### First draw 40 parameters from a normal distribution
  #Random.seed!(1);                        # for sane reproducibility in debugging
  #param = randn(Float64,N);
  #param = ones(Float64,N);
  #param = zeros(N);
  #param[1] = 0; param[2] = 1;
  
  ### Generate the finite element mesh using Gmsh (implemented in makeMesh)
  
  r,err = fitBSpline2Fourier(a0,a,b,N)
  #r = ones(N,1); # uncomment for some degubbing
  @printf("B-spline approximation error is: %12.8f\n",err);
  
  ### Then map them to a distribution between 0.5 and 1.5 using the arctan function, small alpha values (e.g. 0.1) cluster the results of the B-spline parameters around 1
  
  α = 1.0;
  r = 1.0 .+ atan.(α*r)/π;
  
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
  
  
  velocity, pressure = twodStokesRotatingOuter(xT,eC,innerNodes,outerNodes,ω)
  
  #added the mass matrix to the items here
  temperature,A = twodAdvectionDiffusion(xT,eC,innerNodes,outerNodes,velocity, κ);
  
  #  Output the solution or visualize
  scalarLabels = ["temperature"]
  vectorLabels = ["velocity"]
  #saveFEMasVTK("mixing",xT,eC,scalarLabels,temperature,vectorLabels,velocity)
  
  #poly(xT, eC[:,1:3], color = velocity[:,2], strokecolor = (:black, 0.6), strokewidth = .2)
  
  vorticity = computeVorticity(xT,eC,velocity)
  saveFEMasVTK("mixing",xT,eC,scalarLabels,vorticity,vectorLabels,velocity)
  p = poly(xT, eC[:,1:3], color = vorticity[:,1], strokecolor = (:black, 0.6), strokewidth = 0.2)
  save("vort_$(istr).png",p);
  
  velMag = sqrt.( velocity[:,1].*velocity[:,1] + velocity[:,2].*velocity[:,2] )
  #poly(xT, eC[:,1:3], color = velMag, strokecolor = (:black, 0.6), strokewidth = .3)
  
  p = poly(xT, eC[:,1:3], color = temperature[:,1], strokecolor = (:black, 0.6), strokewidth = .2)
  save("temp_$(istr).png",p);
  
  # compute the average temperature, x- and y-components of velocity, and vorticity.
  
  V2 = sum(C)
  Volume = sum(Call)
  
  #temperature
  tMeanSub = ((C*temperature)/V2)[1];
  #@printf("the average temperature over the subregion is %g\n",tMeanSub)
  obs[i,1] = tMeanSub;
  tMean = ((Call*temperature)/Volume)[1];
  #@printf("                    and over the entire region is %g\n\n",tMean)
  obs[i,2] = tMean;
  
  #vorticity
  vortMeanSub = ((C*vorticity)/V2)[1];
  #@printf("the average vorticity over the subregion is %g\n",vortMeanSub)
  obs[i,3] = vortMeanSub;
  vortMean = ((Call*vorticity)/Volume)[1];
  #@printf("                  and over the entire region is %g\n\n\n",vortMean)
  obs[i,4] = vortMean;

  #scalar variance
  tdel = temperature .- tMean;
  obs[i,5] = dot(tdel,A*tdel)/Volume;
end

@printf("%3s: ","#");
@printf("%12s ","T (sub)");   #  1. avg temp (subregion)
@printf("%12s ","T (dom)");   #  2. avg temp (whole domain)
@printf("%12s ","vort (sub)");#  3. avg vort (subregion)
@printf("%12s ","vort (dom)");#  4. avg vort (whole domain)
@printf("%12s ","sv (dom)");  #  5. scalar variance (whole domain)
@printf("\n");
for i=1:nSamples
  @printf("%3d: ",i);
  for j=1:size(obs,2)
    @printf("%12.8f ",obs[i,j]);
  end
  @printf("\n");
end


h5write(outFile,"a",aArr);
h5write(outFile,"b",bArr);
h5write(outFile,"obs",obs);
