#this function has been deprecated in favor of computeRadii
#function computeFourier(ab,th,a0=1.0)
#  a = ab[1:2:end];
#  b = ab[2:2:end];
#
#  #compute fourier representation
#  r = a0.*ones(length(th));
#  for j=1:length(a)
#      r += a[j].*cos.(j*th) + b[j].*sin.(j*th);
#  end
#
#  #squash
#  α = 1.0;
#  r = 1.0 .+ atan.(α*r)/π;
#
#  return r;
#end
