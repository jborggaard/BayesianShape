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

include("../../src/makeMesh.jl")
include("../../src/fitBSpline2Fourier.jl")
include("../../src/saveFEMasVTK.jl")
include("../../src/twodQuadratureRule.jl")
include("../../src/twodShape.jl")
include("../../src/twodBilinear.jl")
include("../../src/twodLinForm.jl")
include("../../src/twodStokesRotatingOuter.jl")
include("../../src/twodAdvectionDiffusion.jl")
include("../../src/twodProjectDerivatives.jl")
include("../../src/computeC.jl")
include("../../src/computeVorticity.jl")
include("../../src/solutionArray.jl");
include("../../src/twodStokesAD.jl");
include("../../src/TriMesh_ElementAdjacency.jl");
include("../../src/TriMesh_Search.jl");
include("../../src/TriMesh_Interpolate.jl");

#using SpectralDiscrete2D
#using AdvectionDiffusion
using InfDimMCMC
#using AdVecMCMC

#ADR_ROOT=ENV["ADR_ROOT"];

#defaults (typically overwritten by arguments to run.jl)
def_datafile  ="dummy";#ADR_ROOT*"/data/point_twohump_012.h5";
def_mcmc  = "pcn|2^-2";
def_nsamp = 10;
def_nburn = 0;

def_omega  = 10.0;
def_kappa  = 1.00;
def_svmean = 0.030;
def_svstd  = 0.0005;

# ## ADVECTION-DIFFUSION PROBLEM DEFINITION: KAPPA, INITIAL CONDITION, OBSERVATIONS/DATA ##
# ad = adProb();
# 
# datafile = (@isdefined datafile) ? datafile : def_datafile;
# @printf("Reading problem setup from datafile: %s\n",datafile);
# maxnorm   = 8;#h5read(datafile,"maxnorm"); 
# ad.kappa  = h5read(datafile,"kappa"); 
# ad.t0     = h5read(datafile,"t0");
# ad.dt     = h5read(datafile,"dt");
# ad.tf     = h5read(datafile,"tf"); 
# ad.th0    = h5read(datafile,"th0"); 
# dataTX    = h5read(datafile,"dataTX");
# y         = h5read(datafile,"y");
# #nsStd     = h5read(datafile,"nsStd");
# 
# #ad.dt = 0.001;
# ad.t      = collect( ad.t0:ad.dt:ad.tf ); ad.nt = length(ad.t);
# 
# #kdisc = specDiscrete2D(maxnorm;multiplicity=true,computeKdiff=true);
# kdisc = specDiscrete2D(maxnorm;multiplicity=true,kDiffFile=@sprintf("%s/lookup/kdiff_kd%03d.h5",ADR_ROOT,maxnorm));
# 
# datafileMaxNorm = h5read(datafile,"maxnorm");
# if datafileMaxNorm > kdisc.maxnorm
#   @printf("datafile maxnorm > kdisc.maxnorm (%d < %d). Truncating th0 to match...\n", datafileMaxNorm, kdisc.maxnorm);
#   ad.th0    = ad.th0[1:kdisc.nk];
# elseif datafileMaxNorm < kdisc.maxnorm
#   @printf("datafile maxnorm < kdisc.maxnorm (%d < %d). Padding th0 with zeros to match...\n", datafileMaxNorm, kdisc.maxnorm);
#   ad.th0    = [ad.th0; zeros(kdisc.nk - length(ad.th0))];
# end
# 
# 
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
omega  = (@isdefined omega)  ? omega  : def_omega;
kappa  = (@isdefined kappa)  ? kappa  : def_kappa;
svMean = (@isdefined svmean) ? svmean : def_svmean;
svStd  = (@isdefined svstd ) ? svstd  : def_svstd;

sourceXY=[1.5;1.0]; #source location

@printf("Using omega=%12.6f and kappa=%12.6f\n",omega,kappa);

#data#
datafile = (@isdefined datafile) ? datafile : def_datafile;
@printf("Using svMean=%12.6f and svStd=%12.6f\n",svMean,svStd);

#dimension of unknown (number of sines and cosines)
unkDim    = 160;
#number of B splines to approximate interior boundary
nBsplines = 160;

#circle centers (subregions)
nCircles= 8;
aInterp = collect(0:(nCircles-1)).*pi/4;
rInterp = 1.75;
xInterp = rInterp .* [ cos.(aInterp) sin.(aInterp) ];

#dummy
dataY = rand(nCircles);
dataStd = 0.05*ones(nCircles);

## Setup MCMC Problem ##

mcmcP = mcmcProb();

# Number of samples # 
mcmcP.nburn = (@isdefined nburn) ? nburn : def_nburn;
mcmcP.nsamp = (@isdefined nsamp) ? nsamp : def_nsamp;

# Define sampler # 
mcmc = (@isdefined mcmc) ? mcmc : def_mcmc;
mcmcSetSampler(mcmcP,mcmc);
mcmcP.computeGradients = false;

## Setup sample space ##

# Prior #
sinCosStd = 2.0.^(-(0:unkDim-1)./4); 
prStd = zeros(2*unkDim);
prStd[1:2:end] = sinCosStd; #cos
prStd[2:2:end] = sinCosStd; #sin
mcmcP.prior = MvNormal(zeros(2*unkDim),prStd);

# Likelihood #
llh = MvNormal(dataY,dataStd);

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
let nBsplines=nBsplines,omega=omega,kappa=kappa,circleCenters=xInterp
  function adSolve(ab)
    a = ab[1:2:end]; 
    b = ab[2:2:end];
    return twodStokesAD(a,b,1.0,nBsplines;ω=omega,κ=kappa,sourceXY=sourceXY,circleCenters=circleCenters);
  end
  InfDimMCMC.mcmcForwardMap(s) = adSolve(s.param);
end

# Observation map #
let xInterp=xInterp
  function obs(sa)
    #use eConn2 to try to guess observation locations
    ecSize = [ size(ec,1) for ec in sa.eConn2 ];
    ecSzCum = [ sum(ecSize[i+1:end]) + round(Int64,ecSize[i]/2) for i=1:length(ecSize) ];
    elGuess = size(sa.eConn,1) .- ecSzCum;

    #interpolate
    eAdjacency = TriMesh_ElementAdjacency(sa.eConn);
    fValues,elInterp = TriMesh_Interpolate(sa.x,sa.eConn,eAdjacency,sa.temperature,xInterp;verbose=1,startEl=elGuess);
    println("Observations:");
    display(fValues)

    return fValues[:];
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
