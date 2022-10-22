
using Gmsh:gmsh
using LinearAlgebra
#using Makie
using CairoMakie
#using AbstractPlotting
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
#include("twodQuadratureRule.jl")
#include("twodShape.jl")
#include("twodBilinear.jl")
#include("twodLinForm.jl")
include("twodStokesRotatingOuter.jl")
include("twodAdvectionDiffusion.jl")
include("twodProjectDerivatives.jl")
include("computeC.jl")
include("computeVorticity.jl")
include("computeFourier.jl")

function plotSampleGrid(samples,outFile; idx=round.(Int,range(size(samples,1)/2, size(samples,1), length=4)), nBsplines = size(samples,2), a0 = 1.0, ω = -10.0, κ = 1.0, sourceXY=[1.5;1.0], circleCenters=[], quivNpts = 1000, quivScale = 0.05, figsize=800, verbose=true)

  verbose && println("omega = $(ω), kappa = $(κ)");

  #p1 = Figure(resolution=(figsize*1.2,figsize));
  p2 = Figure(resolution=(figsize*1.2,figsize));
  #p3 = Plots.quiver(qxy[:,1],qxy[:,2],quiver=(qvf[:,1],qvf[:,2]),marker=(:none),color=:black);

  #limits for shared colorbar
  cmin = 0.0;
  cmax = 0.0;

  for i = 1:length(idx)
    sIdx = idx[i];
    ab = samples[sIdx,:];
    #outFileTmp = @sprintf("%s_s%08d",outFile,sIdx);

    row = (i+1) ÷ 2; col = i+2 - 2*row;

    a = ab[1:2:end]; 
    b = ab[2:2:end];

    sa = twodStokesAD(a,b,1.0,nBsplines;ω=ω,κ=κ,circleCenters=circleCenters);
    
    # vorticity = computeVorticity(sa.x,sa.eConn,sa.velocity)
    # #p1 = poly(sa.x, sa.eConn[:,1:3], color = vorticity[:,1], strokecolor = (:black, 0.6), strokewidth = 0.2, aspect_ratio=:equal)
    # #p1 = Figure(resolution=(figsize*1.2,figsize));
    # #ax1 = p1[1,1] = Axis(p1, aspect=AxisAspect(1), xlabel="x", ylabel="y");
    # #poly!(ax1, sa.x, sa.eConn[:,1:3], color = vorticity[:,1], strokecolor = (:black, 0.6), strokewidth = 0.2);
    # poly(p1[1,1], sa.x, sa.eConn[:,1:3], color = vorticity[:,1], strokecolor = (:black, 0.6), strokewidth = 0.2, axis=(aspect=AxisAspect(1),xlabel="x",ylabel="y"));
    # Colorbar(p1[1,2], width=20, limits = extrema(vorticity[:,1]));
    # plotName = outFileTmp*"_vort.png";
    # save(plotName,p1);
    # println("Wrote: $plotName");
    
    #p2 = poly(sa.x, sa.eConn[:,1:3], color = sa.temperature[:,1], strokecolor = (:black, 0.6), strokewidth = .2, aspect_ratio=:equal)
    #p2 = Figure(resolution=(figsize*1.2,figsize));
    #ax2 = Axis(p2[row,col], xlabel = "x", ylabel = "y", title = "Sample $(sIdx)");
    poly(p2[row,col], sa.x, sa.eConn[:,1:3], color = sa.temperature[:,1], strokecolor = (:black, 0.6), strokewidth = 0.2, axis=(aspect=AxisAspect(1),xlabel="x",ylabel="y",title="Sample $(sIdx)"));
    if i == 1
      #cmin = minimum(sa.temperature[:,1]);
      #cmax = maximum(sa.temperature[:,1]);
      cmin,cmax = extrema(sa.temperature[:,1]);
    else
      cmin = min(cmin,minimum(sa.temperature[:,1]));
      cmax = max(cmax,maximum(sa.temperature[:,1]));
    end
    #Colorbar(p2[1,2], width=20, limits = extrema(sa.temperature[:,1]));
    #plotName = outFileTmp*"_temp.png";
    #save(plotName,p2);
    #println("Wrote: $plotName");


    # #quiver plot
    # pidx = rand(1:size(sa.x,1),quivNpts);
    # qxy = sa.x[pidx,:]; 
    # vf  = sa.velocity[pidx,:];
    # qvf = quivScale .* vf;
    # 
    # th = (0.0:360.0)*pi/180;   #angles in radians
    # #r = computeFourier(ab,th); #inner boundary
    # r = computeRadii(ab,th); #inner boundary

    # #p3 = Plots.quiver(qxy[:,1],qxy[:,2],quiver=(qvf[:,1],qvf[:,2]),marker=(:none),color=:black);
    # Plots.plot!(p3,r.*cos.(th),r.*sin.(th),color=:blue,lab=:none);
    # Plots.plot!(p3,2.0.*cos.(th),2.0.*sin.(th),color=:blue,lab=:none);
    # #Plots.plot!(p3,[1.5],[0.75],markershape=:xcross,markercolor=:red,markersize=15,markerstrokewidth=3,line=false,lab=:none);
    # Plots.plot!(p3,[sourceXY[1]],[sourceXY[2]],markershape=:xcross,markercolor=:red,markersize=15,markerstrokewidth=3,line=false,lab=:none);
    # Plots.plot!(p3,aspect_ratio=:equal,size=(figsize,figsize));
    # plotName = outFileTmp*"_quiver.png";
    # Plots.savefig(p3,plotName);
    # println("Wrote: $plotName");
  end

  Colorbar(p2[:,3], width=20, limits = (cmin,cmax));
  plotName = outFile*"_temp.png";
  save(plotName,p2);
  println("Wrote: $plotName");
  
  #return p1,p2,p3;
  return p2;
end

function plotSampleGrid(inFile::String; kwargs...)
  f = h5open(inFile,"r");
  samples = read(f,"samples");
  close(f);
  outFile = replace(inFile,".h5"=>"_sample_grid");
  plotSampleGrid(samples,outFile; kwargs...);
end
