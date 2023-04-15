
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

function plotSampleGrid(samples,outFile; idx=round.(Int,range(size(samples,1)/2, size(samples,1), length=4)), computeScalar=true, nBsplines = size(samples,2), a0 = 1.0, ω = -10.0, κ = 1.0, sourceXY=[1.5;1.0], circleCenters=[], circleRadius = 0.1, circleColors = [ :black for i=1:size(circleCenters,1) ], quivNpts = 500, quivScale = 0.025, figsize=800, verbose=true)

  verbose && println("omega = $(ω), kappa = $(κ)");

  p1 = Figure(resolution=(figsize*1.2,figsize));
  p2 = Figure(resolution=(figsize*1.2,figsize));
  p3 = Plots.plot(layout=4,aspect_ratio=:equal,size=(figsize,figsize));
  #p3 = Plots.quiver(qxy[:,1],qxy[:,2],quiver=(qvf[:,1],qvf[:,2]),marker=(:none),color=:black);

  #limits for shared colorbar
  cmin = zeros(3);
  cmax = zeros(3);

  #angles
  th = collect(0:360).*pi/180.0;

  for i = 1:length(idx)
    sIdx = idx[i];
    ab = samples[sIdx,:];
    #outFileTmp = @sprintf("%s_s%08d",outFile,sIdx);

    row = (i+1) ÷ 2; col = i+2 - 2*row;

    a = ab[1:2:end]; 
    b = ab[2:2:end];

    if computeScalar
      sa = twodStokesAD(a,b,a0,nBsplines;ω=ω,κ=κ,circleCenters=circleCenters);
    else
      sa = twodStokesOnly(a,b,a0,nBsplines;ω=ω,circleCenters=circleCenters);
    end
    
    vorticity = computeVorticity(sa.x,sa.eConn,sa.velocity)
    # #p1 = poly(sa.x, sa.eConn[:,1:3], color = vorticity[:,1], strokecolor = (:black, 0.6), strokewidth = 0.2, aspect_ratio=:equal)
    # #p1 = Figure(resolution=(figsize*1.2,figsize));
    # #ax1 = p1[1,1] = Axis(p1, aspect=AxisAspect(1), xlabel="x", ylabel="y");
    # #poly!(ax1, sa.x, sa.eConn[:,1:3], color = vorticity[:,1], strokecolor = (:black, 0.6), strokewidth = 0.2);
    # poly(p1[1,1], sa.x, sa.eConn[:,1:3], color = vorticity[:,1], strokecolor = (:black, 0.6), strokewidth = 0.2, axis=(aspect=AxisAspect(1),xlabel="x",ylabel="y"));
    poly(p1[row,col], sa.x, sa.eConn[:,1:3], color = vorticity[:,1], strokecolor = (:black, 0.6), strokewidth = 0.2, axis=(aspect=AxisAspect(1),xlabel="x",ylabel="y",title="Sample $(sIdx)"));
    if i == 1
      cmin[1],cmax[1] = extrema(vorticity[:,1]);
    else
      cmin[1] = min(cmin[1],minimum(vorticity[:,1]));
      cmax[1] = max(cmax[1],maximum(vorticity[:,1]));
    end
    for j=1:size(circleCenters,1)
      circle = circleRadius .* [ cos.(th)  sin.(th) ] .+ circleCenters[j,:]';
      #lines!(p1,circle[:,1],circle[:,2],color=:red,lab=:none);
      lines!(p1[row,col],circle[:,1],circle[:,2],color=circleColors[j],linewidth=2);
    end
    # Colorbar(p1[1,2], width=20, limits = extrema(vorticity[:,1]));
    # plotName = outFileTmp*"_vort.png";
    # save(plotName,p1);
    # println("Wrote: $plotName");
    
    if computeScalar
      #p2 = poly(sa.x, sa.eConn[:,1:3], color = sa.temperature[:,1], strokecolor = (:black, 0.6), strokewidth = .2, aspect_ratio=:equal)
      #p2 = Figure(resolution=(figsize*1.2,figsize));
      #ax2 = Axis(p2[row,col], xlabel = "x", ylabel = "y", title = "Sample $(sIdx)");
      poly(p2[row,col], sa.x, sa.eConn[:,1:3], color = sa.temperature[:,1], strokecolor = (:black, 0.6), strokewidth = 0.2, axis=(aspect=AxisAspect(1),xlabel="x",ylabel="y",title="Sample $(sIdx)"));
      if i == 1
        #cmin = minimum(sa.temperature[:,1]);
        #cmax = maximum(sa.temperature[:,1]);
        cmin[2],cmax[2] = extrema(sa.temperature[:,1]);
      else
        cmin[2] = min(cmin[2],minimum(sa.temperature[:,1]));
        cmax[2] = max(cmax[2],maximum(sa.temperature[:,1]));
      end
      #Colorbar(p2[1,2], width=20, limits = extrema(sa.temperature[:,1]));
      #plotName = outFileTmp*"_temp.png";
      #save(plotName,p2);
      #println("Wrote: $plotName");
    end


    #quiver plot
    pidx = rand(1:size(sa.x,1),quivNpts);
    qxy = sa.x[pidx,:]; 
    vf  = sa.velocity[pidx,:];
    qvf = quivScale .* vf;
    
    r = computeRadii(ab,th); #inner boundary

    Plots.quiver!(p3[i],qxy[:,1],qxy[:,2],quiver=(qvf[:,1],qvf[:,2]),marker=(:none),color=:black);
    Plots.plot!(p3[i],r.*cos.(th),r.*sin.(th),color=:blue,lab=:none);
    Plots.plot!(p3[i],2.0.*cos.(th),2.0.*sin.(th),color=:blue,lab=:none);
    #Plots.plot!(p3,[1.5],[0.75],markershape=:xcross,markercolor=:red,markersize=15,markerstrokewidth=3,line=false,lab=:none);
    if computeScalar
      Plots.plot!(p3[i],[sourceXY[1]],[sourceXY[2]],markershape=:xcross,markercolor=:red,markersize=15,markerstrokewidth=3,line=false,lab=:none);
    end
    for j=1:size(circleCenters,1)
      circle = circleRadius .* [ cos.(th)  sin.(th) ] .+ circleCenters[j,:]';
      Plots.plot!(p3[i],circle[:,1],circle[:,2],lc=:red,lab=:none,linewidth=2);
    end
    Plots.plot!(p3[i],xlabel="x",ylabel="y",title="Sample $(sIdx)");
    #Plots.plot!(p3[i],aspect_ratio=:equal,size=(figsize,figsize));
  end

  Colorbar(p1[:,3], width=20, limits = (cmin[1],cmax[1]));
  plotName = outFile*"_vort.png";
  save(plotName,p1);
  println("Wrote: $plotName");

  if computeScalar
    Colorbar(p2[:,3], width=20, limits = (cmin[2],cmax[2]));
    plotName = outFile*"_temp.png";
    save(plotName,p2);
    println("Wrote: $plotName");
  end

  plotName = outFile*"_quiver.png";
  Plots.savefig(p3,plotName);
  println("Wrote: $plotName");
  plotName = outFile*"_quiver.pdf";
  Plots.savefig(p3,plotName);
  println("Wrote: $plotName");
  
  return p1,p2,p3;
end

function plotSampleGrid(inFile::String; labelwith="", kwargs...)
  f = h5open(inFile,"r");
  samples = read(f,"samples");
  if labelwith == "likelihood"
    labelvalues = read(f,"lpdfs")[:,2];
    labelprefix = "Likelihood = ";
  end
  close(f);
  outFile = replace(inFile,".h5"=>"_sample_grid");
  plotSampleGrid(samples,outFile; kwargs...);
end
