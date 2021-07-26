using Plots
using HDF5
using LinearAlgebra

function plotMapIBs(mapSamples::AbstractArray; kwargs...)
  labels=["MPR","MLE","MAP"];
  
  #gr();
  p = plot(proj=:polar,size=(900,300),layout=(1,3),leg=false,kwargs...);
  
  #angles
  th = pi*(0:360)/180;#range(0.0,stop=2.0*pi,length=size(samples,2));


  for i=1:3
      #get sample
      #ab = samples[maxIdx[i][1],:];
      ab = mapSamples[i,:];
  
      #compute radii
      r = computeRadii(ab,th);
  
      #plot
      plot!(p[i], th, r, c=:black);
      plot!(p[i], th, 2.0.*ones(length(th)), c=:black);
      plot!(p[i], title=labels[i]);
  end
  return p;
end

function plotMapIBs(mapSamples::AbstractArray,outFile::String; exts=["png"], kwargs...)
  p = plotMapIBs(mapSamples; kwargs...);
  plotSave(p,outFile,exts); 
end

function plotMapIBs(inFile::String; kwargs...)
  mapSamples = getMap(inFile,0); 
  outFile = replace(inFile,".h5"=>"_map_ibs");
  plotMapIBs(mapSamples,outFile; kwargs...);
end
