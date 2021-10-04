function polarTriangle(vert::Array{Float64,2}, pts::AbstractArray)
  #slopes
  delv = vert[[2;3;1], :] - vert;
  m = delv[:,2] ./ delv[:,1];
  #intercepts
  b = vert[:,2] - m .* vert[:,1];
  
  #rescale pts to [0,1]
  ptsr = (pts .- pts[1]) ./ (pts[end]-pts[1]);

  #now polar
  thmin = atan.(vert[:,2],vert[:,1]);
  thmax = thmin[[2;3;1]]; #thmax[3] += 2*pi;
  thmax += 2*pi*( thmax .< thmin );
 
  #angles
  th = ptsr * (thmax - thmin)' .+ thmin';

  #radii
  r  = b' ./ ( sin.(th) - m'.*cos.(th) );
 
  return th, r;
end

# original script to build and plot
#
# using Plots, Statistics
#
# #vertices
# #vert = [ 0 0; 1 0; 0.5 1.0 ];
# vert = 5.0.*[ 0 0; 1 0; 0.635 0.275 ];
# #centroid
# cent = mean(vert,dims=1);
# #recenter
# vert0 = vert .- cent;
# #slopes
# delv = vert0[[2;3;1], :] - vert0;
# m = delv[:,2] ./ delv[:,1];
# #intercepts
# b = vert0[:,2] - m .* vert0[:,1];
# 
# #plot to check that we got the slopes right
# xmin,xmax = extrema(vert0[:,1]);
# x = collect(range(xmin,xmax,length=5));
# lines = m' .* x .+ b';
# plot(x,lines,leg=false); #lines
# scatter!(vert0[:,1],vert0[:,2]); #vertices
# 
# #now polar
# thmin = atan.(vert0[:,2],vert0[:,1]);
# thmax = thmin[[2;3;1]]; thmax[3] += 2*pi;
# #if given th, here's how to figure out which interval it lives in:
# # ((th .>= thmin') .& (th .< thmax')) .| ((th .>= (thmin.+2*pi)') .& (th .< (thmax.+2*pi)'))
# npts = 10;
# th = zeros(npts,3);
# r  = zeros(npts,3);
# p = plot(proj=:polar,leg=false);
# for i = 1:3
#   th[:,i] = collect(range(thmin[i],thmax[i],length=npts));
#   r[:,i]  = b[i] ./ ( sin.(th[:,i]) - m[i].*cos.(th[:,i]) );
#   plot!(p,th[:,i],r[:,i]);
# end
# scatter!(p,thmin,sqrt.(vert0[:,1].^2 + vert0[:,2].^2)); #add points
# display(p)
