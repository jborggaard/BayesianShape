#squashSmoothstep() uses atan() to squeeze radii to a given (rMin,rMax)
#see https://en.wikipedia.org/wiki/Smoothstep
function squashSmoothstep(r, rMin, rMax)
  width = rMax-rMin;
  rMean = 0.5*(rMin+rMax);

  #default derivative is not 1 at the midpoint, this corrects for that
  #d = 1.5; #order 3
  d = 30/16; #order 5   
  rc = (r .- rMean)./d .+ rMean; 
  
  #rescale to [0,1]
  rc = collect( (rc.-rMin)./width );
  
  #clamp
  rc[ rc .< 0.0 ] .= 0.0;
  rc[ rc .> 1.0 ] .= 1.0;

  #hermite polynomial
  #rs = 3.0 * rc.^2 - 2.0*rc.^3; #order 3
  rs = 6*rc.^5 - 15*rc.^4 + 10*rc.^3;  #order 5

  #rescale back to [rMin,rMax]
  rs = rMin .+ width.* rs;
  
  return rs;
end

