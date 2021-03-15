
using Gmsh:gmsh
using LinearAlgebra
#using Makie
using CairoMakie
using AbstractPlotting
using SparseArrays
using SpecialMatrices
using Plots
using Polynomials
using Printf
using Random
using WriteVTK
using HDF5

include("makeMesh.jl")
include("fitBSpline2Fourier.jl")
include("saveFEMasVTK.jl")
#include("sampleInnerGeometry.jl")
include("twodQuadratureRule.jl")
include("twodShape.jl")
include("twodBilinear.jl")
include("twodLinForm.jl")
include("twodStokesRotatingOuter.jl")
include("twodAdvectionDiffusion.jl")
include("twodProjectDerivatives.jl")
include("computeC.jl")
include("computeVorticity.jl")
include("computeFourier.jl")

function plotSample(ab,outFile; N = length(ab), a0 = 1.0, ω = -10.0, κ = 1.0, circleCenters=[], verbose=true, quivNpts = 1000, quivScale = 0.05)

  verbose && println("omega = $(ω), kappa = $(κ)");

  a = ab[1:2:end]; 
  b = ab[2:2:end];

  # ### Generate the finite element mesh using Gmsh (implemented in makeMesh)
  # 
  # r,err = fitBSpline2Fourier(a0,a,b,N)
  # #r = ones(N,1); # uncomment for some degubbing
  # @printf("B-spline approximation error (%d B-splines) is: %12.8f\n",N,err);
  # 
  # ### Then map them to a distribution between 0.5 and 1.5 using the arctan function, small alpha values (e.g. 0.1) cluster the results of the B-spline parameters around 1
  # α = 1.0;
  # r = 1.0 .+ atan.(α*r)/π;
  # 
  # x,eConn,eConn2, innerNodes,innerX, outerNodes,outerX = makeMesh(r)
  # 
  # sort!(innerNodes);
  # sort!(outerNodes);
  # 
  # #   Let's look at the mesh...
  # nNodes = size(x,2)
  # nElements = size(eConn,2)
  # nElementsSensor = size(eConn2,2)
  # 
  # xT = zeros(Float64,nNodes,2)   # xT = transpose(x[1:2,:])
  # eC = zeros(Int64,nElements,6)  # eC = transpose(eConn), converted to Int64
  # for i=1:nNodes
  #   xT[i,1] = x[1,i]
  #   xT[i,2] = x[2,i]
  # end
  # for i=1:nElements#-nElementsSensor # the sensor elements are orientated correctly (why?)
  # #  for j=1:6
  # #    eC[i,j] = convert(Int64,eConn[j,i])
  # #  end
  #    eC[i,1] = convert(Int64,eConn[1,i])
  #    eC[i,2] = convert(Int64,eConn[3,i])
  #    eC[i,3] = convert(Int64,eConn[2,i])
  #    eC[i,4] = convert(Int64,eConn[6,i])
  #    eC[i,5] = convert(Int64,eConn[5,i])
  #    eC[i,6] = convert(Int64,eConn[4,i])
  # end
  # 
  # eC2  = zeros(Int64,nElementsSensor,6)
  # for i=1:nElementsSensor
  #   for j=1:6
  #     eC2[i,j] = convert(Int64,eConn2[j,i])
  #   end
  # end  # elements are positively orientated
  # C    = computeC(xT,eC2)
  # 
  # Call = computeC(xT,eC)
  # 
  # #compute and plot vorticity
  # velocity = twodStokesRotatingOuter(x,eConn,innerNodes,outerNodes,ω);
  #  
  # #compute and plot temperature
  # temperature,A = twodAdvectionDiffusion(x,eConn,innerNodes,outerNodes,velocity,κ);
  
  
  #sa = twodStokesAD(a,b,1.0,nBsplines;ω=omega,κ=kappa,circleCenters=circleCenters);
  sa = twodStokesAD(a,b,1.0,nBsplines;ω=omega,κ=kappa);
  
  vorticity = computeVorticity(sa.x,sa.eConn,sa.velocity)
  p1 = poly(sa.x, sa.eConn[:,1:3], color = vorticity[:,1], strokecolor = (:black, 0.6), strokewidth = 0.2, aspect_ratio=:equal)
  plotName = outFile*"_vort.png";
  save(plotName,p1);
  println("Wrote: $plotName");
  
  p2 = poly(sa.x, sa.eConn[:,1:3], color = sa.temperature[:,1], strokecolor = (:black, 0.6), strokewidth = .2, aspect_ratio=:equal)
  plotName = outFile*"_temp.png";
  save(plotName,p2);
  println("Wrote: $plotName");


  #quiver plot
  pidx = rand(1:size(sa.x,1),quivNpts);
  qxy = sa.x[pidx,:]; 
  vf  = sa.velocity[pidx,:];
  qvf = quivScale .* vf;
  
  th = (0.0:360.0)*pi/180;   #angles in radians
  r = computeFourier(ab,th); #inner boundary

  p3 = Plots.quiver(qxy[:,1],qxy[:,2],quiver=(qvf[:,1],qvf[:,2]),marker=(:none),color=:black);
  Plots.plot!(p3,r.*cos.(th),r.*sin.(th),color=:blue,lab=:none);
  Plots.plot!(p3,2.0.*cos.(th),2.0.*sin.(th),color=:blue,lab=:none);
  Plots.plot!(p3,[1.5],[0.75],markershape=:xcross,markercolor=:red,markersize=15,markerstrokewidth=3,line=false,lab=:none);
  plotName = outFile*"_quiver.png";
  Plots.savefig(p3,plotName);
  println("Wrote: $plotName");
  
  return p1,p2,p3;
end
