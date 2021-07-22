function circleEVs(nev)
  if nev>13
    error("Only 13 eigenvalues are currently available for the circle problem. Talk to Jeff to get more. ;-)");
  end

  lam = zeros(13);
  lam[1]  = besselj_zero(0,1)^2;
  lam[2]  = besselj_zero(1,1)^2;
  lam[3]  = besselj_zero(1,1)^2;
  lam[4]  = besselj_zero(2,1)^2;
  lam[5]  = besselj_zero(2,1)^2;
  lam[6]  = besselj_zero(0,2)^2;
  lam[7]  = besselj_zero(3,1)^2;
  lam[8]  = besselj_zero(3,1)^2;
  lam[9]  = besselj_zero(1,2)^2;
  lam[10] = besselj_zero(1,2)^2;
  lam[11] = besselj_zero(4,1)^2;
  lam[12] = besselj_zero(4,1)^2;
  lam[13] = besselj_zero(2,2)^2;

  ##this was a draft of semi-hacky code to compute a more general list
  ##but it failed because of the repeated roots
  #lam = zeros(100);
  #cnt = 0;
  #for nu=0:6
  #  for n=1:6
  #    cnt += 1;
  #    lam[cnt]  = besselj_zero(nu,n)^2;
  #  end
  #end
  #lam = lam[1:cnt];
  #sort!(lam);

  return lam[1:nev];
end
