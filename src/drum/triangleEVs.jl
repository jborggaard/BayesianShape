function isoEVs(nev)
  #  Writes out the first 30 eigenvalues that were computed using 
  #  runTriangleExample.  These are eigenvalues associated with the
  #  triangular drum in "Any three eigenvalues do not determine a
  #  triangle," by J. Gomez-Serrano and G. Orriols, arXiv 1911.06758v2.
  #  See case A in Section 3.
  #
 
  if nev>50
    error("Only 50 eigenvalues are currently available for the triangle problem. Talk to Jeff to get more. ;-)");
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
  lam[31] = 3680.141151
  lam[32] = 3739.531294
  lam[33] = 3968.362045
  lam[34] = 4014.653006
  lam[35] = 4080.193257
  lam[36] = 4134.924514
  lam[37] = 4354.085871
  lam[38] = 4404.379962
  lam[39] = 4529.760298
  lam[40] = 4599.451684
  lam[41] = 4732.60348
  lam[42] = 4850.174512
  lam[43] = 4950.449534
  lam[44] = 5033.889615
  lam[45] = 5131.098288
  lam[46] = 5254.38567
  lam[47] = 5319.00187
  lam[48] = 5455.787778
  lam[49] = 5542.936332
  lam[50] = 5598.991273

  return lam[1:nev];
end
