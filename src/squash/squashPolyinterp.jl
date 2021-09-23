#squashPolyinterp() uses polynomial interpolation to squash within epsilon of rMin and rMax
function squashPolyinterp(r, rMin, rMax; e=0.1)
  rs = collect(copy(r));

  #squash: soft lower bound
  b = 0.25/e; c = 0.5; d = 0.25*e;
  fl(x) = b*x.^2 + c*x .+ d;
  idx = abs.(r .- rMin) .<= e;
  rs[idx] = fl( r[idx] .- rMin ) .+ rMin;

  #squash: soft upper bound
  b = -0.25/e; c = 0.5; d = -0.25*e;
  fu(x) = b*x.^2 + c*x .+ d;
  idx = abs.(r .- rMax) .<= e;
  rs[idx] = fu( r[idx] .- rMax ) .+ rMax;

  #hard bounds on (-inf,rMin-e) and (rMax+e,inf)
  rs[ rs.<rMin ] .= rMin;
  rs[ rs.>rMax ] .= rMax;

  return rs;
end

