using Plots
using HDF5
using LinearAlgebra
using Statistics
using Printf

#plot
function plotQuantiles(pltData::AbstractArray; nburn=0, ps=[0.1,0.25,0.5,0.75,0.9], targetData = [], targetLabel="target", kwargs...)
  #compute quantiles
  q = zeros(size(pltData,2),length(ps));
  for j=1:size(q,1)
     q[j,:] = quantile(pltData[:,j],ps);
  end

  #plot
  p = plot(;kwargs...);
  if length(targetData) > 0
      plot!(p, 1:size(pltData,2), targetData, lab=targetLabel, lc=:black, ls=:dash);
  end
  for j=1:length(ps)
      label = @sprintf("%d%%",100*ps[j]);
      plot!(p, 1:size(pltData,2), q[:,j], lab=label);
  end

  return p;
end

#plots and saves
function plotQuantiles(pltData::AbstractArray,outFile::String; exts=["png"], kwargs...)
  p = plotQuantiles(pltData; kwargs...);
  plotSave(p,outFile,exts); 
end

#read in, plot, and save
function plotQuantiles(inFile::String; kwargs...)
  f = h5open(inFile,"r");
  pltData    = read(f,"obs");
  targetData = read(f,"obsMean");
  close(f);
  outFile = replace(inFile,".h5"=>"_obs_quantiles");
  plotQuantiles(pltData,outFile; targetData=targetData, targetLabel="y", kwargs...);
end


