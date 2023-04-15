using Plots
using Plots.Measures
using HDF5
using LinearAlgebra

function plotSamplesVsTruth(samples,obs,trueSamp,trueObs; 
  idx=round.(Int,range(size(samples,1)/2, size(samples,1), length=4)), 
  #computeScalar=true, nBsplines = size(samples,2), a0 = 1.0, ω = -10.0, κ = 1.0, sourceXY=[1.5;1.0], 
  circleCenters=[], circleRadius = 0.1, circleColors = [ :grey for i=1:size(circleCenters,1) ], 
  figsize=800, verbose=true,
  obsXlab="Observation #", obsYlab="Observation Value",
  kwargs...)

  sampLabels = ["Sample $(i)" for i=idx];
  
  #angles
  #th = pi*(0:360)/180;
  th = collect(0:360).*pi/180.0;

  #gr();
  p = plot(size=(2*figsize,figsize),layout=(1,2),kwargs...);

  ## first plot (radii) ##
  for j=1:size(circleCenters,1)
    circle   = circleRadius .* [ cos.(th)  sin.(th) ] .+ circleCenters[j,:]';
    #plot!(p[1],circle[:,1],circle[:,2],color=circleColors[j],lab=nothing);
    circleR  = sqrt.( circle[:,1].^2 + circle[:,2].^2 );
    circleTh = atan.( circle[:,2], circle[:,1] );
    plot!(p[1],circleTh,circleR,color=circleColors[j],lab=nothing);
  end
  plot!(p[1],margin=10mm,proj=:polar,leg=true);
  plot!(p[1], th, 2.0.*ones(length(th)), c=:black, label=nothing); #outer boundary
  plot!(p[2],margin=10mm,xlab=obsXlab, ylab=obsYlab);

  #plot truth
  trueRadii = computeRadii(trueSamp,th);
  plot!(p[1], th, trueRadii, lc=:black, ls=:dash, lab="Truth");
  plot!(p[2], 1:length(trueObs), trueObs, seriestype=:scatter, c=:black, markershape=:x, ms=8, lab="Truth");
  #plot!(p[2], 1:length(trueObs), trueObs, lc=:black, lab="Truth");

  for i=1:length(idx)
      #get sample & observation
      ab = samples[idx[i],:];
      sampObs = obs[idx[i],:];
  
      #compute radii
      r = computeRadii(ab,th);
  
      #plot
      plot!(p[1], th, r, c=i, lab=sampLabels[i]);
      plot!(p[2], 1:length(sampObs), sampObs, seriestype=:scatter, c=i, lab=sampLabels[i]);
      #plot!(p[2], 1:length(sampObs), sampObs, c=i, lab=sampLabels[i]);
  end
  return p;
end

function plotSamplesVsTruth(samples::AbstractArray,obs::AbstractArray,trueSamp::AbstractArray,trueObs::AbstractArray,outFile::String; exts=["png"], kwargs...)
  p = plotSamplesVsTruth(samples, obs, trueSamp, trueObs; kwargs...);
  plotSave(p,outFile,exts); 
end

function plotSamplesVsTruth(inFile::String; kwargs...)
  samples  = h5read(inFile,"samples");
  obs      = h5read(inFile,"obs");
  trueSamp = h5read(inFile,"trueSamp");
  trueObs  = h5read(inFile,"obsMean");
  outFile = replace(inFile,".h5"=>"_samples_vs_truth");
  plotSamplesVsTruth(samples,obs,trueSamp,trueObs,outFile; kwargs...);
end
