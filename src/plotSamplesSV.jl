using Plots
using HDF5

function plotSamplesSV(obs,svMean,svStd;nthin = 10, kwargs...)
  y = svMean; ystd = svStd;
  idx = 1:nthin:size(obs,1);
  p = plot(xlab="Sample #",ylab="Scalar Variance"; kwargs...);
  plot!(p,idx,obs[idx,1],lab="G(u)");
  plot!(p,idx,y*ones(length(idx)),lab="y");
  plot!(p,idx,y*ones(length(idx)).-ystd,color=3,ls=:dash,lab="y +/- std");
  plot!(p,idx,y*ones(length(idx)).+ystd,color=3,ls=:dash,lab=:none);
  plot!(p,leg=:bottomright);
  return p
end

function plotSamplesSV(obs,svMean,svStd,outFile;nthin = 10,ext=["png"], kwargs...)
  p = plotSamplesSV(obs,svMean,svStd; nthin=nthin, kwargs...);
  #for ex in ext
  #  oFl = outFile*"."*ex;
  #  savefig(p,oFl);
  #  println("Wrote: $oFl");
  #end
  plotSave(p,outFile,ext);
end

function plotSamplesSV(inFile,outFile; nthin=10,ext=["png"], kwargs...)
  #read in data
  f = h5open(inFile);
  obs    = read(f,"obs");
  svMean = read(f,"svMean");
  svStd  = read(f,"svStd");
  close(f);
  #plot
  plotSamplesSV(obs,svMean,svStd,outFile;nthin=nthin,ext=ext, kwargs...);
end

function plotSamplesSV(inFile;nthin = 10,ext=["png"], kwargs...)
  outFile = replace(inFile,".h5"=>"_sv_evolve");
  plotSamplesSV(inFile,outFile;nthin=nthin,ext=ext, kwargs...);
end
