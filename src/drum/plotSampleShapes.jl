#single sample
function plotSampleShapes(ab::Array{Float64,1}; rMin=0.5, rMax=1.5, kwargs...)

  sz = 300;
  p = plot(proj=:polar,size=(sz,sz),leg=false);
  
  #angles
  th = pi*(0:360)/180;#range(0.0,stop=2.0*pi,length=size(samples,2));
  
  #compute radii
  r = computeRadii(ab,th;rMin=rMin,rMax=rMax);
  
  #plot
  plot!(p, th, r, c=:black);  

  return p;
end

function plotSampleShapes(samples::AbstractArray; idx=sort(rand(1:size(samples,1),9)), rMin=0.5, rMax=1.5, kwargs...)

  sz = sqrt(length(idx))*300;
  p = plot(proj=:polar,size=(sz,sz),layout=(length(idx)),leg=false);
  
  #angles
  th = pi*(0:360)/180;#range(0.0,stop=2.0*pi,length=size(samples,2));
  
  for i=1:length(idx)
    #get sample
    ab = samples[idx[i],:];
  
    #compute radii
    r = computeRadii(ab,th;rMin=rMin,rMax=rMax);
  
    #plot
    plot!(p[i], th, r, c=:black);  
    plot!(p[i], title="Sample $(idx[i])");
  end
  return p;
end

function plotSampleShapes(samples::AbstractArray, outFile::String; exts=["png"], idx=sort(rand(1:size(samples,1),9)), kwargs...)
  p = plotSampleShapes(samples; kwargs...);
  plotSave(p,outFile,exts);
end

function plotSampleShapes(inFile::String; kwargs...)
  f = h5open(inFile,"r");
  samples = read(f,"samples");
  rMin = read(f,"rMin");
  rMax = read(f,"rMax");
  close(f);
  outFile = replace(inFile,".h5"=>"_sample_shapes");
  plotSampleShapes(samples,outFile; rMin=rMin, rMax=rMax, kwargs...);
end
