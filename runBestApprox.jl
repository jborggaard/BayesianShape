a0 = 1.0; a = [-1/2, 1/4, -1/8, 1/8, 1/8]; b = [1/4, 1/8, 0.0, 1/4, 1/8];
#plotFourier(a0,a,b)
#hold on
N = 40;
coefficients = fitBSpline2Fourier(a0,a,b,N);
#plotBSpline(coefficients,N)