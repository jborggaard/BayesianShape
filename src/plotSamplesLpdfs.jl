using Plots
using HDF5

function plotSamplesLpdfs(lpdfs::AbstractArray;nthin=10,kwargs...)
  idx = 1:nthin:size(lpdfs,1);
  labels=["Log prior","Log LLH","Log (u)posterior"];
  p = plot(layout=(3,1),size=(600,600),leg=false;kwargs...);
  for i=1:length(p)
    plot!(p[i],idx,lpdfs[idx,i],xlab="Sample #",ylab=labels[i]);
  end
  return p
end

function plotSamplesLpdfs(lpdfs::AbstractArray,outFile::String;nthin = 10,ext=["png"],kwargs...)
  p = plotSamplesLpdfs(lpdfs; nthin=nthin,kwargs...);
  plotSave(p,outFile,ext);
end

function plotSamplesLpdfs(inFile::String,outFile::String; nthin=10,ext=["png"],kwargs...)
  #read in data
  f = h5open(inFile);
  lpdfs    = read(f,"lpdfs");
  close(f);
  #plot
  plotSamplesLpdfs(lpdfs,outFile;nthin=nthin,ext=ext,kwargs...);
end

function plotSamplesLpdfs(inFile::String;nthin = 10,ext=["png"],kwargs...)
  outFile = replace(inFile,".h5"=>"_lpdf_evolve");
  plotSamplesLpdfs(inFile,outFile;nthin=nthin,ext=ext,kwargs...);
end
