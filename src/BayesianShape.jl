module BayesianShape

using ArgParse
using Arpack
using CairoMakie: poly, lines
using Distributions
using FastGaussQuadrature
using FunctionZeros   # for the validation example 
using Gmsh:gmsh
using HDF5
using LinearAlgebra
using Plots
using Plots.Measures
using Polynomials
using Printf
using Random
using SparseArrays
using SpecialMatrices
using Statistics
using WriteVTK

#don't think we need these anymore
#using Makie
#using AbstractPlotting

using FEMfunctions
using InfDimMCMC

#miscellaneous functions
export computeC;
export computeFourier;
export computeRadii;
export computeVorticity;
export fitBSpline2Fourier;
export fourierBasis;
export generateSampleObs;
export getMap;
export makeMesh;
export solutionArray;
include("computeC.jl");
include("computeFourier.jl");
include("computeRadii.jl");
include("computeVorticity.jl");
include("fitBSpline2Fourier.jl");
include("fourierBasis.jl");
include("generateSampleObs.jl");
include("getMap.jl");
include("makeMesh.jl");
include("solutionArray.jl");

#solvers
export twodAdvectionDiffusion;
export twodNavierStokesAD;
export twodNavierStokesRotatingOuter;
export twodNavierStokesRotatingOuterNewton;
export twodStokesAD;
export twodStokesOnly;
export twodStokesRotatingOuter;
include("twodAdvectionDiffusion.jl");
include("twodNavierStokesAD.jl");
include("twodNavierStokesRotatingOuter.jl");
include("twodNavierStokesRotatingOuterNewton.jl");
include("twodStokesAD.jl");
include("twodStokesOnly.jl");
include("twodStokesRotatingOuter.jl");

#drum-specific functions
export circleEVs;
export computeFEMmatrices;
export inputOutput;
export isoEVs;
export makeDrumMesh;
export triangleEVs;
include("drum/circleEVs.jl");
include("drum/computeFEMmatrices.jl");
include("drum/inputOutput.jl");
include("drum/isoEVs.jl");
include("drum/makeDrumMesh.jl");
include("drum/triangleEVs.jl");
export triangleRecenter;
export trianglePolar;
include("triangle/triangleRecenter.jl");
include("triangle/trianglePolar.jl");

#squash functions
export squashArctan
export squashErf
export squashPolyinterp
export squashSigmoid
export squashSmoothstep
include("squash/squashArctan.jl");
include("squash/squashErf.jl");
include("squash/squashPolyinterp.jl");
include("squash/squashSigmoid.jl");
include("squash/squashSmoothstep.jl");

export radiusSquash; #placeholder to be overwritten by scenarios
radiusSquash(r) = error("You must define the radius squash.");

#plotting functions
export histmatrix;
export plotMapIBs;
export plotMap;
export plotQuantiles;
export plotRadiiCorr;
export plotRadiiHist;
export plotRadiiHistMatrix;
export plotRadiiQuantiles;
export plotSampleGrid;
export plotSample;
export plotSamplesIBs;
export plotSamplesLpdfs;
export plotSamplesSV;
export plotSamplesVsTruth;
export plotSave;
include("histmatrix.jl");
include("plotMapIBs.jl");
include("plotMap.jl");
include("plotQuantiles.jl");
include("plotRadiiCorr.jl");
include("plotRadiiHist.jl");
include("plotRadiiHistMatrix.jl");
include("plotRadiiQuantiles.jl");
include("plotSampleGrid.jl");
include("plotSample.jl");
include("plotSamplesIBs.jl");
include("plotSamplesLpdfs.jl");
include("plotSamplesSV.jl");
include("plotSamplesVsTruth.jl");
include("plotSave.jl");
export plotDrumMap;
export plotSampleShapes;
include("drum/plotDrumMap.jl");
include("drum/plotSampleShapes.jl");

#don't think we need these anymore
#export sampleInnerGeometry;
#include("sampleInnerGeometry.jl");

#I believe these are now in FEMfunctions
#export TriMesh_ElementAdjacency;
#export TriMesh_Interpolate;
#export TriMesh_Search;
#export twodBilinear;
#export twodLinForm;
#export twodMassMatrix;
#export twodProjectDerivatives;
#export twodQuadratureRule;
#export twodShape;
#
#include("TriMesh_ElementAdjacency.jl");
#include("TriMesh_Interpolate.jl");
#include("TriMesh_Search.jl");
#include("twodBilinear.jl");
#include("twodLinForm.jl");
#include("twodMassMatrix.jl");
#include("twodProjectDerivatives.jl");
#include("twodQuadratureRule.jl");
#include("twodShape.jl");


end # module
