#using Plots
#using HDF5
#using LinearAlgebra

function plotSamplesIBs(samples::AbstractArray, svs::AbstractArray; a0=0.0, idx=round.(Int,range(1, size(samples,1), length=9)), kwargs...)
  
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
      #r = computeFourier(ab,th);
      r = computeRadii(ab,th,a0);
  
      #plot
      plot!(p[i], th, r, c=:black);
      plot!(p[i], th, 2.0.*ones(length(th)), c=:black);
      sv = round(svs[idx[i]]; digits=5);
      plot!(p[i], title="$(idx[i]): $sv");
  end
  return p;
end

function plotSamplesIBs(samples::AbstractArray,svs::AbstractArray,a0,outFile::String; exts=["png"], kwargs...)
  p = plotSamplesIBs(samples,svs,a0=a0; kwargs...);
  plotSave(p,outFile,exts); 
end

function plotSamplesIBs(inFile::String; kwargs...)
  f = h5open(inFile,"r");
  samples = read(f,"samples");
  a0 = read(f,"a0");
  svs = read(f,"obs")[:,1];
  close(f);
  outFile = replace(inFile,".h5"=>"_sample_ibs");
  plotSamplesIBs(samples,svs,a0,outFile; kwargs...);
end

