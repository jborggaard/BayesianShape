using Plots
using HDF5
using LinearAlgebra
using Statistics
using Printf

#plot
function plotRadiiHistMatrix(samples::AbstractArray; th=pi*(0:45:359)/180, nburn=0, kwargs...)
  
  nAngles=length(th);

  thDeg = round.(Int,th*180/pi);

  #compute radius by angle
  fb = fourierBasis(size(samples,2)÷2,th);
  sampleAngles = computeRadii(samples[nburn+1:end,:],fb);
  
  #plot histograms
  p = histmatrix(sampleAngles; label=[ "$t" for t in thDeg ], kwargs...);

  return p;
end

#plots and saves
function plotRadiiHistMatrix(samples::AbstractArray,outFile::String; exts=["png"], kwargs...)
  p = plotRadiiHistMatrix(samples; kwargs...);
  plotSave(p,outFile,exts); 
end

#read in, plot, and save
function plotRadiiHistMatrix(inFile::String; kwargs...)
  f = h5open(inFile,"r");
  samples = read(f,"samples");
  close(f);
  outFile = replace(inFile,".h5"=>"_radii_histmatrix");
  plotRadiiHistMatrix(samples,outFile; kwargs...);
end


