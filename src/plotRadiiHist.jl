using Plots
using HDF5
using LinearAlgebra
using Statistics
using Printf

#plot
function plotRadiiHist(samples::AbstractArray; th=pi*(0:45:359)/180, nburn=0, plotType=:graphical, kwargs...)
  
  nAngles=length(th);

  if plotType == :graphical && nAngles != 8
    warn("plotType $(plotType) not supported with the number of angles != 8. Changing to plotType=:sequential.");
    plotType = :sequential
  end

  thDeg = round.(Int,th*180/pi);

  #compute radius by angle
  fb = fourierBasis(size(samples,2)÷2,th);
  sampleAngles = computeRadii(samples[nburn+1:end,:],fb);
  
  #plot histograms
  if plotType == :single
    #all histograms in a single plot
    p = stephist(sampleAngles, bins = :scott, labels=[ "$(ang)" for ang=thDeg ] );
  elseif plotType == :sequential
    #separate plots, sequential order
    p = Plots.plot(layout=(nAngles),leg=false,xlab="Radius",xlim=extrema(sampleAngles[:]));
    for i=1:nAngles
      stephist!(p[i],sampleAngles[:,i], bins = :scott, ylab="$(thDeg[i])");
    end
  elseif plotType == :graphical
    #separate plots, try to lay out along with angles
    p = Plots.plot(layout=(9),leg=false,xlab="Radius",xlim=extrema(sampleAngles[:]));
    pltIndex = [ 6; 3; 2; 1; 4; 7; 8; 9 ];
    for i=1:nAngles
      stephist!(p[pltIndex[i]],sampleAngles[:,i], bins = :scott, ylab="$(thDeg[i])");
    end
  else
    error("plotType $(plotType) not recognized.");
  end

  return p;
end

#plots and saves
function plotRadiiHist(samples::AbstractArray,outFile::String; exts=["png"], kwargs...)
  p = plotRadiiHist(samples; kwargs...);
  plotSave(p,outFile,exts); 
end

#read in, plot, and save
function plotRadiiHist(inFile::String; kwargs...)
  f = h5open(inFile,"r");
  samples = read(f,"samples");
  close(f);
  outFile = replace(inFile,".h5"=>"_radii_hist");
  plotRadiiHist(samples,outFile; kwargs...);
end


