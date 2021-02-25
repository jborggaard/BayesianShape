function twodBilinear( kernel, phi, test, wg )
##----------------------------------------------------------------------
#  twodBilinear - routine to compute \int{ kernel*phi*test }
#
#  Copyright (c) 2001, Jeff Borggaard, Virginia Tech
#  Version: 1.0
#
#  Usage:    M = twodBilinear(kernel, phi, test, wg)
#
#  Variables:     kernel
#                        Kernel function in the integral evaluated
#                        at the Gauss points
#
#                 phi
#                        matrix of element test functions evaluated
#                        at the Gauss points (dim: n_gauss, n_dof)
#
#                 test
#                        matrix of test functions evaluated at the
#                        Gauss points (dim: n_gauss, n_dof)
#
#                 wg
#                        Row vector of quadrature weights
## ---------------------------------------------------------------------

#  M = test'*diagm( wg.*kernel )*phi;


  wk = wg.*kernel;

  M = test'*(LinearAlgebra.Diagonal( wk )*phi);

  return M

end
