using Plots
using HDF5
using LinearAlgebra

function plotSamplesIBs(samples::AbstractArray, svs::AbstractArray; idx=[1,10,50,100,1000,2000,3000,4000,5000], kwargs...)
  
  #truncate indices if we don't have enough samples
  idx = idx[idx .<= size(samples,1)];
  
  #gr();
  sz = sqrt(length(idx))*300;
  p = plot(proj=:polar,size=(sz,sz),layout=(length(idx)),leg=false,kwargs...);
  
  #angles
  th = pi*(0:360)/180;#range(0.0,stop=2.0*pi,length=size(samples,2));

  for i=1:length(idx)
      #get sample
      ab = samples[idx[i],:];
  
      #compute fourier representation
      r = computeFourier(ab,th);
  
      #plot
      plot!(p[i], th, r, c=:black);
      plot!(p[i], th, 2.0.*ones(length(th)), c=:black);
      sv = round(svs[idx[i]]; digits=5);
      plot!(p[i], title="$(idx[i]): $sv");
  end
  return p;
end

function plotSamplesIBs(samples::AbstractArray,svs::AbstractArray,outFile::String; exts=["png"], kwargs...)
  p = plotSamplesIBs(samples,svs; kwargs...);
  plotSave(p,outFile,exts); 
end

function plotSamplesIBs(inFile::String; kwargs...)
  f = h5open(inFile,"r");
  samples = read(f,"samples");
  svs = read(f,"obs")[:,1];
  close(f);
  outFile = replace(inFile,".h5"=>"_sample_ibs");
  plotSamplesIBs(samples,svs,outFile; kwargs...);
end

