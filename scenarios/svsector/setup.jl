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
include("../../src/sampleInnerGeometry.jl")
#include("../../src/twodQuadratureRule.jl")
#include("../../src/twodShape.jl")
#include("../../src/twodMassMatrix.jl")
#include("../../src/twodBilinear.jl")
#include("../../src/twodLinForm.jl")
include("../../src/twodStokesRotatingOuter.jl")
include("../../src/twodAdvectionDiffusion.jl")
include("../../src/twodProjectDerivatives.jl")
include("../../src/computeC.jl")
include("../../src/computeVorticity.jl")
include("../../src/solutionArray.jl");
include("../../src/twodStokesAD.jl");

#using SpectralDiscrete2D
#using AdvectionDiffusion
using InfDimMCMC
#using AdVecMCMC

#ADR_ROOT=ENV["ADR_ROOT"];

#defaults (typically overwritten by arguments to run.jl)
def_datafile  ="dummy";#ADR_ROOT*"/data/point_twohump_012.h5";
def_mcmc  = "pcn|2^-2";
def_ar    = 0.25;
def_nsamp = 10;
def_nburn = 0;

nSectors = 4;

def_regularity = 1.0; #want samples in H_s for s < regularity
def_omega  = 10.0;
def_kappa  = 1.00;
def_rmin   = 0.5;
def_rmax   = 1.5;
def_a0     = 1.0;
def_svmean = [0.4; 0.0; 0.4; 0.0];
def_svstd  = 0.05;
    
# data=adPointData(ad,kdisc,dataTX,y;computeL2Kern=false);

#parameters
#omega  = 10.0;
#kappa  = 1.00;
#svMean = 0.030;
#svStd  = 0.0005;
#kappa  = 0.10;
#svMean = 0.07; #0.02;
#svStd  = 0.01;
#kappa  = 0.01;
#svMean = 0.70; #0.10;
#svStd  = 0.10;
regularity = (@isdefined regularity) ? regularity : def_regularity;
omega  = (@isdefined omega)  ? omega  : def_omega;
kappa  = (@isdefined kappa)  ? kappa  : def_kappa;
rMin    = (@isdefined rmin   ) ? rmin   : def_rmin;
rMax    = (@isdefined rmax   ) ? rmax   : def_rmax;
a0      = (@isdefined a0     ) ? a0     : def_a0;
svMean = (@isdefined svmean) ? svmean : def_svmean;
svStd  = (@isdefined svstd ) ? svstd  : def_svstd;

sourceXY=[1.5;1.0]; #source location

@printf("Using omega=%12.6f and kappa=%12.6f\n",omega,kappa);
println("Regularity = $(regularity)");

#data#
datafile = (@isdefined datafile) ? datafile : def_datafile;
#@printf("Using svMean=%12.6f and svStd=%12.6f\n",svMean,svStd);

#dimension of unknown (number of sines and cosines)
unkDim    = 160;
#number of B splines to approximate interior boundary
nBsplines = 160;

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
llh = MvNormal(svMean,svStd);

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
let nBsplines=nBsplines,omega=omega,kappa=kappa,a0=a0
  function adSolve(ab)
    a = ab[1:2:end]; 
    b = ab[2:2:end];
    return twodStokesAD(a,b,a0,nBsplines;ω=omega,κ=kappa,sourceXY=sourceXY);
  end
  InfDimMCMC.mcmcForwardMap(s) = adSolve(s.param);
end

# Observation map #
let nSectors=nSectors
  function obs(sa)
    # #find element centroids
    # centroids  = s.sol.xT[s.sol.eC[:,1],:];
    # centroids += s.sol.xT[s.sol.eC[:,2],:];
    # centroids += s.sol.xT[s.sol.eC[:,3],:];
    # centroids ./= 3.0;

    # #figure out element sectors
    # centAngles = atan.(centroids[:,2],centroids[:,1]); 
    # centAngles .+= 2*pi*( centAngles .< 0 ); #(-pi,pi) => (0,2pi)
    # centSector = floor.(Int,nSectors*centAngles./(2*pi));

    #find node sectors
    xAngles   = atan.(sa.x[:,2],sa.x[:,1]); 
    xAngles .+= 2*pi*( xAngles .< 0 ); #(-pi,pi) => (0,2pi)
    xSector   = floor.(Int,nSectors*xAngles./(2*pi)) .+ 1; #ceil produces zero for (2,0)
    
    #mean temp and deviation
    Call     = computeC(sa.x,sa.eConn);
    Volume   = sum(Call);
    tMean    = ((Call*sa.temperature)/Volume)[1];
    tdel = sa.temperature .- tMean;
    
    #compute scalar variance by sector
    massMat = twodMassMatrix(sa.x,sa.eConn);
    sv   = zeros(nSectors+1);
    sv[1] = dot(tdel,massMat*tdel)/Volume;
    for i=1:nSectors
      idx   = (xSector .== i);
      sv[i+1] = dot(tdel[idx],massMat[idx,idx]*tdel[idx])/Volume;
    end
    
    return sv;
  end
  InfDimMCMC.mcmcObsMap(s) = obs(s.sol);
end

# Potential map #
let llh=llh
  InfDimMCMC.mcmcPotMap(s) = -logpdf(llh,s.obs[2:end]);
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
