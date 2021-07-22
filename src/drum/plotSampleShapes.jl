function plotSampleShapes(samples::AbstractArray; idx=sort(rand(1:size(samples,1),9)), kwargs...)

  sz = sqrt(length(idx))*300;
  p = plot(proj=:polar,size=(sz,sz),layout=(length(idx)),leg=false);
  
  #angles
  th = pi*(0:360)/180;#range(0.0,stop=2.0*pi,length=size(samples,2));
  
  for i=1:length(idx)
    #get sample
    ab = samples[idx[i],:];
  
    #compute fourier representation
    r = computeFourier(ab,th);
  
    #plot
    plot!(p[i], th, r, c=:black);  
    plot!(p[i], title="Sample $(idx[i])");
  end
  return p;
end
