using Plots
using HDF5
using LinearAlgebra
using Statistics
using Printf

#plot
function plotRadiiQuantiles(samples::AbstractArray; th=pi*(0:360)/180, nburn=0, ps=[0.1,0.25,0.5,0.75,0.9], kwargs...)
  #compute radius by angle
  fb = fourierBasis(size(samples,2)÷2,th);
  sampleAngles = computeRadii(samples[nburn+1:end,:],fb);
  
  #compute quantiles
  q = zeros(size(sampleAngles,2),length(ps));
  for j=1:size(q,1)
     q[j,:] = quantile(sampleAngles[:,j],ps);
  end

  #plot
  sz = 500;
  p = plot(proj=:polar,size=(sz,sz);kwargs...);
  for j=1:length(ps)
      label = @sprintf("%d%%",100*ps[j]);
      plot!(p, th, q[:,j], lab=label);
  end

  return p;
end

#plots and saves
function plotRadiiQuantiles(samples::AbstractArray,outFile::String; exts=["png"], kwargs...)
  p = plotRadiiQuantiles(samples; kwargs...);
  plotSave(p,outFile,exts); 
end

#read in, plot, and save
function plotRadiiQuantiles(inFile::String; kwargs...)
  f = h5open(inFile,"r");
  samples = read(f,"samples");
  close(f);
  outFile = replace(inFile,".h5"=>"_radii_quantiles");
  plotRadiiQuantiles(samples,outFile; kwargs...);
end


