function isoEVs(nev)
  #  Writes out the first 30 eigenvalues that were computed using 
  #  runTriangleExample.  These are eigenvalues associated with the
  #  triangular drum in "Any three eigenvalues do not determine a
  #  triangle," by J. Gomez-Serrano and G. Orriols, arXiv 1911.06758v2.
  #  See case A in Section 3.
  #
 
  if nev>30
    error("Only 30 eigenvalues are currently available for the triangle problem. Talk to Jeff to get more. ;-)");
  end

  lam = zeros(30);
  lam[ 1] =  233.4680566
  lam[ 2] =  391.4452799
  lam[ 3] =  546.1468169
  lam[ 4] =  698.9388722
  lam[ 5] =  790.8456206
  lam[ 6] =  906.2111733
  lam[ 7] = 1096.400444
  lam[ 8] = 1123.579023
  lam[ 9] = 1347.672249
  lam[10] = 1360.963646
  lam[11] = 1554.866078
  lam[12] = 1603.927283
  lam[13] = 1727.497587
  lam[14] = 1881.806161
  lam[15] = 1927.809755
  lam[16] = 2142.927905
  lam[17] = 2182.742933
  lam[18] = 2248.155789
  lam[19] = 2478.75141
  lam[20] = 2536.988674
  lam[21] = 2545.328114
  lam[22] = 2754.529731
  lam[23] = 2845.17537
  lam[24] = 2909.648902
  lam[25] = 3046.57542
  lam[26] = 3192.060662
  lam[27] = 3274.594339
  lam[28] = 3335.977015
  lam[29] = 3525.072572
  lam[30] = 3583.034133

  return lam[1:nev];
end
