#this function has been deprecated in favor of computeRadii
#function computeFourier(ab,th; rMin=0.2, rMax=5.0)
#  a = ab[1:2:end];
#  b = ab[2:2:end];
#
#  #compute fourier representation
#  r = a0.*ones(length(th));
#  for j=1:length(a)
#      r += a[j].*cos.(j*th) + b[j].*sin.(j*th);
#  end
#  fb = fourierBasis(length(ab)÷2,th);
#  sampleAngles = computeRadii(samples[nburn+1:end,:],fb);
#
#  #α = 1.0;
#  #r0 = 0.5*(rMax+rMin);
#  #r = r0 .+ (rMax-rMin)*atan.(α*(r.-r0))/π;
#  ##squash
#  #α = 1.0;
#  #r = 1.0 .+ atan.(α*r)/π;
#  r = radiusSquash(r; rMin=rMin, rMax=rMax);
#
#  return r;
#end
