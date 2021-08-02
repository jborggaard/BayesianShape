function isoEVs(nev)
  #  Writes out the first 30 eigenvalues that were computed using runME1
  #
 
  if nev>30
    error("Only 30 eigenvalues are currently available for the isospectral problem. Talk to Jeff to get more. ;-)");
  end

  lam = zeros(30);
  lam[ 1] =  6.540813909 
  lam[ 2] =  6.601734625 
  lam[ 3] = 10.86849532 
  lam[ 4] = 11.06411183 
  lam[ 5] = 11.18288225 
  lam[ 6] = 11.58919436 
  lam[ 7] = 13.60955176 
  lam[ 8] = 13.60959103 
  lam[ 9] = 13.60964353 
  lam[10] = 13.96635365 
  lam[11] = 14.66879231 
  lam[12] = 14.9068866 
  lam[13] = 15.66712442 
  lam[14] = 16.06653466 
  lam[15] = 17.46166006 
  lam[16] = 18.42146297 
  lam[17] = 19.72903042 
  lam[18] = 20.96754459 
  lam[19] = 21.80938303 
  lam[20] = 23.49062664 
  lam[21] = 25.71041391 
  lam[22] = 25.71094091 
  lam[23] = 25.71138743 
  lam[24] = 27.46260917 
  lam[25] = 29.48354845 
  lam[26] = 30.15003172 
  lam[27] = 32.19788286 
  lam[28] = 34.05536264 
  lam[29] = 35.62081149 
  lam[30] = 36.73464678

  return lam[1:nev];
end
