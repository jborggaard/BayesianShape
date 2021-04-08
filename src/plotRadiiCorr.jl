using Plots
using HDF5
using LinearAlgebra
using Statistics
using Printf

#plot
function plotRadiiCorr(samples::AbstractArray; th=pi*(0:360)/180, nburn=0, kwargs...)
  nAngles = length(th);

  thDeg = round.(Int,th*180/pi);

  #compute radius by angle
  fb = fourierBasis(size(samples,2)÷2,th);
  sampleAngles = computeRadii(samples[nburn+1:end,:],fb);
  
  #compute correlations
  idx = [1;91;181;271];
  c = zeros(nAngles,length(idx));
  for i=1:nAngles
    lagIdx = mod.(idx .+ i .-1,nAngles).+1;
    for j=1:length(idx)
      c[i,j] = cor(sampleAngles[:,idx[j]],sampleAngles[:,lagIdx[j]]);
    end
  end

  #plot
  #sz = 500;
  #p = plot(proj=:polar,size=(sz,sz),kwargs...);
  p = plot(xlab="Lag",ylab="Correlation",ylim=(-1,1),leg=:bottomright,xticks=0:90:360; kwargs...);
  for j=1:length(idx)
      label = @sprintf("%d",thDeg[idx[j]]);
      plot!(p, thDeg, c[:,j], lab=label);
  end

  return p;
end

#plots and saves
function plotRadiiCorr(samples::AbstractArray,outFile::String; exts=["png"], kwargs...)
  p = plotRadiiCorr(samples; kwargs...);
  plotSave(p,outFile,exts); 
end

#read in, plot, and save
function plotRadiiCorr(inFile::String; kwargs...)
  f = h5open(inFile,"r");
  samples = read(f,"samples");
  close(f);
  outFile = replace(inFile,".h5"=>"_radii_corr");
  plotRadiiCorr(samples,outFile; kwargs...);
end


