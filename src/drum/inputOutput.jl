# N:   number of BSplines used to represent the drum shape
# nev: how many frequencies we hear
function inputOutput(a0,a,b; N=40, nev = 1, κ=1.0, rMin=0.2, rMax=5.0)

  r,err = fitBSpline2Fourier(a0,a,b,N)
  #@printf("B-spline approximation error (%d B-splines) is: %12.8f\n",N,err)
  #println("Extrema of r (pre-smoothing) are:");
  #display(extrema(r));
  
  ### Then map them to a distribution between rMin and rMax using the arctan function, small alpha values (e.g. 0.1) cluster the results of the B-spline parameters around 1
  α = 1.0;
  r0 = 0.5*(rMax+rMin);
  r = r0 .+ (rMax-rMin)*atan.(α*(r.-r0))/π;
  #println("Extrema of r (post-smoothing) are:");
  #display(extrema(r));

  x,eConn,boundaryNodes = makeDrumMesh(r)

  nBoundary = length(boundaryNodes)
  dBoundary = zeros(Float64,nBoundary,1)
  A,M = computeFEMmatrices(x,eConn,boundaryNodes,dBoundary;κ=κ)

  λ, ϕ = eigs(A,M; tol=1e-6, which=:SM, nev)

  realλ = real(λ)

  return realλ

end
