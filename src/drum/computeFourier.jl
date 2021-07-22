function computeFourier(ab,th;a0=1.0, rMin=0.2, rMax=5.0)
  a = ab[1:2:end];
  b = ab[2:2:end];

  #compute fourier representation
  r = a0.*ones(length(th));
  for j=1:length(a)
      r += a[j].*cos.(j*th) + b[j].*sin.(j*th);
  end

  α = 1.0;
  r0 = 0.5*(rMax+rMin);
  r = r0 .+ (rMax-rMin)*atan.(α*(r.-r0))/π;
  ##squash
  #α = 1.0;
  #r = 1.0 .+ atan.(α*r)/π;

  return r;
end
