function twodLinForm( Ff, test, wg )
#-----------------------------------------------------------------------
#  twodLinForm - routine to compute \int{ f*test }
#
#  Copyright (c) 2001, Jeff Borggaard, Virginia Tech
#  Version: 1.0
#
#  Usage:    F = twodLinForm( Ff, test, wg )
#
#  Variables:     Ff    
#                        Function values at the Gauss points
#
#                 test
#                        matrix of test functions evaluated at the
#                        Gauss points (dim: n_gauss, n_dof)
#
#                 wg
#                        Row vector of Gauss weights
#-----------------------------------------------------------------------

# n_dof = size(test,2);

# F = zeros(n_dof,1);
# for j=1:n_dof
#   F(j) = test(:,j)' * ( wg .* Ff );
# end
# F = sum(test'*(wg.*Ff),2);

  F = test'*(wg.*Ff);

  return F

end
