using Gmsh:gmsh
using LinearAlgebra
##using Makie
#using CairoMakie
#using AbstractPlotting
using SparseArrays
using SpecialMatrices
#using Plots
using Polynomials
using Printf
using Random
using WriteVTK
using Distributions
using HDF5

using FEMfunctions

include("../../src/makeMesh.jl")
include("../../src/fitBSpline2Fourier.jl")
include("../../src/saveFEMasVTK.jl")
include("../../src/twodStokesRotatingOuter.jl")
include("../../src/twodAdvectionDiffusion.jl")
include("../../src/twodProjectDerivatives.jl")
include("../../src/computeC.jl")
include("../../src/computeVorticity.jl")
include("../../src/solutionArray.jl");
include("../../src/twodStokesOnly.jl");

using InfDimMCMC

#ADR_ROOT=ENV["ADR_ROOT"];

#defaults (typically overwritten by arguments to run.jl)
def_datafile  ="dummy";#ADR_ROOT*"/data/point_twohump_012.h5";
def_mcmc  = "pcn|2^-2";
def_ar    = 0.25;
def_nsamp = 10;
def_nburn = 0;

def_regularity = 1.25; #want samples in H_s for s < regularity
def_omega  = 10.0;
def_rmin   = 0.5;
def_rmax   = 1.5;
def_a0     = 1.0;
#def_obsmean = 40.0 .- 10.0*cos.( 0.5*pi*collect(0:7) );
def_obsmean = [30.0; 40.0; 50.0; 40.0; 40.0; 40.0; 30.0; 50.0 ];
def_obsstd  = 1.0;

#parameters
regularity = (@isdefined regularity) ? regularity : def_regularity;
omega  = (@isdefined omega)  ? omega  : def_omega;
rMin    = (@isdefined rmin   ) ? rmin   : def_rmin;
rMax    = (@isdefined rmax   ) ? rmax   : def_rmax;
a0      = (@isdefined a0     ) ? a0     : def_a0;
obsMean = (@isdefined obsmean) ? obsmean : def_obsmean;
obsStd  = (@isdefined obsstd ) ? obsstd  : def_obsstd;

@printf("Using omega=%12.6f\n",omega);
println("Regularity = $(regularity)");

#data#
datafile = (@isdefined datafile) ? datafile : def_datafile;
@printf("Using obsMean:\n");
display(obsMean);
@printf("Using obsStd:\n");
display(obsStd);

#dimension of unknown (number of sines and cosines)
unkDim    = 160;
#number of B splines to approximate interior boundary
nBsplines = 160;

#circle centers (subregions)
nCircles= length(obsMean);
aInterp = collect(0:(nCircles-1)).*pi/4;
rInterp = 1.75;
circleCenters = rInterp .* [ cos.(aInterp) sin.(aInterp) ];


#define squash methodology
squashE = 0.1;
include("../../src/squash/squashPolyinterp.jl");
radiusSquash(r) = squashPolyinterp(r,rMin,rMax;e=squashE);
squashMethod="squashPolyinterp, e=$(squashE)";
println("Squashing with: $(squashMethod)");

## Setup MCMC Problem ##

mcmcP = mcmcProb();

# Number of samples # 
mcmcP.nburn = (@isdefined nburn) ? nburn : def_nburn;
mcmcP.nsamp = (@isdefined nsamp) ? nsamp : def_nsamp;

# Define sampler # 
mcmc = (@isdefined mcmc) ? mcmc : def_mcmc;
mcmcSetSampler(mcmcP,mcmc);
mcmcP.computeGradients = false;
targetAR = (@isdefined ar) ? ar : def_ar;
println("targetAR = $(targetAR)");

## Setup sample space ##

# Prior #
p = 2*regularity + 1; #see Dashti-Stuart Thm 2.12 
sinCosStd = (1:unkDim).^(-0.5*p); #2.0.^(-(0:unkDim-1)./4); 
#sinCosStd = 2.0.^(-(0:unkDim-1)./4); 
prStd = zeros(2*unkDim);
prStd[1:2:end] = sinCosStd; #cos
prStd[2:2:end] = sinCosStd; #sin
mcmcP.prior = MvNormal(zeros(2*unkDim),prStd);

# Likelihood #
llh = MvNormal(obsMean,obsStd);

# # Map from samples to vector components #
# let nk=kdiscV.nk, sampInd=sampInd
#   function padZeros(s)
#     p = zeros(nk);
#     p[sampInd] = s;
#     return p;
#   end
#   InfDimMCMC.mcmcSampToParamMap(s) = padZeros(s.samp);
# end
# let nk=kdiscV.nk
#   gspMat = [ zeros(2,nSampInd); I ];
#   InfDimMCMC.mcmcGradSampToParamMap(s) = gspMat;
# end

# Forward map and observations #
let nBsplines=nBsplines,omega=omega,circleCenters=circleCenters,a0=a0
  function adSolve(ab)
    a = ab[1:2:end]; 
    b = ab[2:2:end];
    return twodStokesOnly(a,b,a0,nBsplines;ω=omega,circleCenters=circleCenters);
  end
  InfDimMCMC.mcmcForwardMap(s) = adSolve(s.param);
end

# Observation map #
let
  function obs(sa)
    #compute vorticity
    vorticity = computeVorticity(sa.x,sa.eConn,sa.velocity);

    #mean vorticity (subregions)
    vortMeanSub = zeros(length(sa.eConn2));
    for i=1:length(sa.eConn2)
      C2 = computeC(sa.x,sa.eConn2[i]);
      V2 = sum(C2);
      vortMeanSub[i] = ((C2*vorticity)/V2)[1];
    end

    return vortMeanSub;
  end
  InfDimMCMC.mcmcObsMap(s) = obs(s.sol);
end

# Potential map #
let llh=llh
  InfDimMCMC.mcmcPotMap(s) = -logpdf(llh,s.obs);
end

# # Gradient of potential map #
# let
#   lolsrc = "../../../lolmc/src";
#   
#   romFile=ADR_ROOT*"/lolmc/rom_gradpot_trig_klsgp.h5";
#   f = h5open(romFile,"r");
#   Ssgprn     = read(f,"Ssgprn");
#   Vsgprn     = read(f,"Vsgprn");
#   sMeansgprn = read(f,"sMeansgprn");
#   kb         = read(f,"kb");
#   alph       = read(f,"alph");
#   close(f);
#   
#   #Use KL built on samples and gradPhi for dimension reduction
#   include(lolsrc*"/kl.jl");
#   include(lolsrc*"/reduceKL.jl");
#   #reduceDim(samp;thresh=1.0,verbose=false) = reduceDim(samp,Ssgprn,Vsgprn,sMeansgprn;thresh=thresh,verbose=verbose);
#   #reduceDim(samp;verbose=false) = reduceDim(samp,Ssgprn,Vsgprn,sMeansgprn;thresh=1.0,verbose=verbose);
#   idataDel = Ssgprn .> 1.0;
#   Sr = Ssgprn[ idataDel ];
#   Vr = Vsgprn[ :, idataDel];
#   global reduceDim(samp;verbose=false) = reduceDim(samp,Sr,Vr,sMeansgprn;verbose=verbose);
# 
#   #Basis functions
#   include(lolsrc*"/basisTrig.jl");
# 
#   include(lolsrc*"/computeGradVandermond.jl");
#   include(lolsrc*"/reducedGradPhi.jl");
#   #gradPotMap(s) = reducedGradPhi(s.samp,kb,alph);
#   InfDimMCMC.mcmcGradPotMap(s) = transpose(s.gradSampToParam)*reducedGradPhi(s.param,kb,alph);
# end
# 
# #TEST
# s = mcmcSample(); s.samp = rand(mcmcP.prior); mcmcFillSample(s,mcmcP);
# #@printf("norm(G(sTrue)-y) = %12.8f\n",norm(s.obs-y));
