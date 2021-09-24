function fitBSpline2Fourier(a0,a,b,N)
# Reads in Fourier coefficients and returns the coefficients of
# the best-L2 BSpline approximation.

  m = zeros(2*N-1);
  m[N-3:N+3] = pi*[1.0/2520.0, 1.0/21.0, 397.0/840.0, 302.0/315.0, 397.0/840.0, 1.0/21.0, 1.0/2520.0]/N;

  m[(2*N-3):(2*N-1)] = pi*[1.0/2520.0, 1.0/21.0, 397.0/840.0]/N;
  m[1:3] = pi*[397.0/840.0, 1.0/21.0, 1.0/2520.0]/N;

  M = Toeplitz(m);

  r = a0 * 2*pi/N *ones(N);

  for n=1:length(a)
    rhs_a = computeRHSa(n,N);
    rhs_b = computeRHSb(n,N);
    r = r + a[n]*rhs_a + b[n]*rhs_b;
  end

  #return coefficients only
  #return M\r

  #return coefficients and the error in the approximation
  #err^2 = ||f||^2 -2 r^T x + x^T M x
  #best approximating coefficients solve Mx=r
  x = M\r; 
  err = sqrt( max(0.0,2*pi*a0^2 + pi*(norm(a)^2+norm(b)^2) - 2*dot(r,x) + dot(x,M*x)) ); #max() because numerical errors can make err^2 < 0
  return x,err;
end

#in this version, a0,a,b are passed in as part of the same vector,
#so we split them and pass them to the original version
function fitBSpline2Fourier(ab,N)
  a0 = ab[1];
  a  = ab[2:2:end];
  b  = ab[3:2:end];
  return fitBSpline2Fourier(a0,a,b,N);
end

function computeRHSa(n,N)

  rhs = zeros(N);

  for k=1:N
    rhs[k] = (N^3*(cos((2*n*pi*(-1 + k))/N) - 4*cos((2*pi*k*n)/N) + 6*cos((2*n*pi*(1 + k))/N) - 4*cos((2*n*pi*(2 + k))/N) + cos((2*n*pi*(3 + k))/N)))/(8*n^4*pi^3);
  end

  return rhs
end

function computeRHSb(n,N)

  rhs = zeros(N);
  for k=1:N
    rhs[k] = (N^3*(sin((2*n*pi*(-1 + k))/N) + 6*sin((2*n*pi*(1 + k))/N) - 4*sin((2*n*pi*(2 + k))/N) + sin((2*n*pi*(3 + k))/N) - 4*sin((2*pi*k*n)/N)))/(8*n^4*pi^3);
  end

  return rhs
end
