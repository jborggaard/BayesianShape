using Gmsh:gmsh
using LinearAlgebra
using Printf
using SparseArrays
using Arpack
using FunctionZeros   # for the validation example
using Distributions
using HDF5

using FEMfunctions

include("../../src/drum/fitBSpline2Fourier.jl")
include("../../src/drum/computeFEMmatrices.jl")
include("../../src/drum/makeDrumMesh.jl")
include("../../src/drum/fitBSpline2Fourier.jl")
#include("../../src/radiusSquash.jl")
include("../../src/drum/inputOutput.jl")
include("../../src/drum/triangleEVs.jl")

using InfDimMCMC

#defaults (typically overwritten by arguments to run.jl)
def_datafile  ="dummy";#ADR_ROOT*"/data/point_twohump_012.h5";
def_mcmc  = "pcn|2^-2";
def_nsamp = 1000;
def_nburn = 100;

def_regularity = 1.25; #want samples in H_s for s < regularity
def_nev    = 20;
def_kappa  = 1.00;
def_rmin   = 0.2;
def_rmax   = 3.0;#9.0;
def_a0     = 1.5;
def_lc     = 7e-3;
    
regularity = (@isdefined regularity) ? regularity : def_regularity;
nEigVals   = (@isdefined nev    )    ? nev        : def_nev;
kappa      = (@isdefined kappa  )    ? kappa      : def_kappa;
rMin       = (@isdefined rmin   )    ? rmin       : def_rmin;
rMax       = (@isdefined rmax   )    ? rmax       : def_rmax;
a0         = (@isdefined a0     )    ? a0         : def_a0;
lc         = (@isdefined lc     )    ? lc         : def_lc;

#def_obsmean = triangleEVs(def_nev); #inputOutput(1.0,zeros(2),zeros(2);nev=def_nev,κ=def_kappa); #zeros(def_nev);
#def_obsstd  = sqrt.(sqrt.(def_obsmean));
#obsMean    = (@isdefined obsmean)    ? obsmean    : def_obsmean;
#obsStd     = (@isdefined obsstd )    ? obsstd     : def_obsstd;
obsMean = triangleEVs(nEigVals); #inputOutput(1.0,zeros(2),zeros(2);nev=def_nev,κ=def_kappa); #zeros(def_nev);
obsMean ./= obsMean[1]; #normalize to eliminate scaling
obsStd  = 0.01*obsMean; #0.08*obsMean; #obsMean.^0.25;

println("Identifying drum shape with triangle A data\n\n");
println("Regularity = $(regularity)");
@printf("Using nev=%12.6f and kappa=%12.6f\n",nEigVals,kappa);
println("Radius constraints: ($(rMin),$(rMax))");
println("Mean radius (a0) = $(a0)");
println("Mesh resolution (lc) = $(lc)");

#data#
datafile = (@isdefined datafile) ? datafile : def_datafile;
#@printf("Using obsMean=%12.6f and obsStd=%12.6f\n",obsMean,obsStd);
println("obsMean is:");
display(obsMean)
println("");

#dimension of unknown (number of sines and cosines)
unkDim    = 160;
#number of B splines to represent boundary
nBsplines = 160;

#define squash methodology
squashE = 0.1;
include("../../src/squash/squashPolyinterp.jl");
radiusSquash(r) = squashPolyinterp(r,rMin,rMax;e=squashE);
squashMethod="squashPolyinterp, e=$(squashE)";
println("Squashing with: $(squashMethod)");


# Setup MCMC Problem ##

mcmcP = mcmcProb();

# Number of samples # 
mcmcP.nburn = (@isdefined nburn) ? nburn : def_nburn;
mcmcP.nsamp = (@isdefined nsamp) ? nsamp : def_nsamp;

# Define sampler # 
mcmc = (@isdefined mcmc) ? mcmc : def_mcmc;
mcmcSetSampler(mcmcP,mcmc);
mcmcP.computeGradients = false;

## Setup sample space ##

sampInd  = 3:(2*unkDim+1);
nSampInd = length(sampInd);

#sampComp = :sincos; #sample sines and cosines
#sampComp = :cos;    #sample only cosines

# Prior #
p = 2*regularity + 1; #see Dashti-Stuart Thm 2.12 
sinCosStd = (1:unkDim).^(-0.5*p); #2.0.^(-(0:unkDim-1)./4); 
prStd = zeros(2*unkDim+1);
prStd[2:2:end] = sinCosStd; #cos
prStd[3:2:end] = sinCosStd; #sin
prStd = prStd[sampInd]; #truncate to sampling indices
mcmcP.prior = MvNormal(zeros(length(prStd)),prStd);

# Likelihood #
llh = MvNormal(obsMean,obsStd);

# Map from samples to vector components #
if nSampInd != 2*unkDim+1
  sampNoInd = setdiff(1:2*unkDim+1,sampInd); #fixed coefficients
  let unkDim=unkDim, nSampInd=nSampInd, sampInd=sampInd, a0=a0
    function spm(s)
      p = zeros(2*unkDim+1);
      p[1] = a0;
      p[sampInd] = s;
      return p;
    end
    InfDimMCMC.mcmcSampToParamMap(s) = spm(s.samp);
  end
  let unkDim=unkDim, nSampInd=nSampInd, sampInd=sampInd
    gspMat = zeros(2*unkDim+1,nSampInd);
    gspMat[CartesianIndex.(sampInd,1:nSampInd)] .=  1.0;
    InfDimMCMC.mcmcGradSampToParamMap(s) = gspMat;
  end
end

# Forward map and observations #
let nBsplines=nBsplines,nEigVals=nEigVals,kappa=kappa,lc=lc
  function drumSolve(ab)
    evs = inputOutput(ab;N=nBsplines,nev=nEigVals,κ=kappa,lc=lc);
    return ( evs ./ evs[1] );
  end
  InfDimMCMC.mcmcForwardMap(s) = drumSolve(s.param);
end

# Observation map #
#let nSectors=nSectors
#  InfDimMCMC.mcmcObsMap(s) = obs(s.sol);
#end

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
