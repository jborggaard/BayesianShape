function inputOutput(a0,a,b)

  N = 40  # number of BSplines used to represent the drum shape


  r,err = fitBSpline2Fourier(a0,a,b,N)
  #@printf("B-spline approximation error (%d B-splines) is: %12.8f\n",N,err)

  x,eConn,boundaryNodes = makeDrumMesh(r)

  nBoundary = length(boundaryNodes)
  dBoundary = zeros(Float64,nBoundary,1)
  A,M = computeFEMmatrices(x,eConn,boundaryNodes,dBoundary)

  #  how many frequencies do we hear?
  nev = 1
  λ, ϕ = eigs(A,M; tol=1e-6, which=:SM, nev)

  realλ = zeros(Float64,nev)
  for i=1:nev
    realλ[i] = real(λ[i])
  end

  return realλ

end
