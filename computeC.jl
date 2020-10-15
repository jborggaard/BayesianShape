function computeC(x,eConn)
#  Given a list of nodes and a subset of elements, computes the matrix 
#  representation of the integral of a finite element solution over that
#  subset of elements.
#
#  Usage:    C = computeC(x,eConn)
#
#  The average of the finite element solution with nodal values z over
#  the region described by the elements in eConn is C*z
#
#  Note that the average over the domain and the area of the domain
#  can be readily found using this function as well.
#
#  An internal parameter is the number of quadrature points which needs
#  to be sufficient to exactly integrate a finite element basis function.
##

  nNodes    = size(x    ,1)
  nElements = size(eConn,1)
  nElDOF    = size(eConn,2)

  rule  = 3
  r,s,w = twodQuadratureRule(rule)
  o     = ones(rule)

  C   = zeros(Float64,1,nNodes)

  # preallocate storage for variables defined at quadrature points
  xg  = zeros(Float64,rule,2)
  wg  = zeros(Float64,rule)
  phi = zeros(Float64,rule,nElDOF)
  p_x = zeros(Float64,rule,nElDOF)
  p_y = zeros(Float64,rule,nElDOF)

  index = 0
  for k=1:nElements

    nLocal = eConn[k,:][:]
    xLocal = x[nLocal,:]

    xg,wg,phi,p_x,p_y = twodShape( xLocal, r, s, w )

    c = twodLinForm(o,phi,wg)    # o can be replaced by a weighting function...

    for nt = 1:nElDOF
      nTest = nLocal[nt]
      C[nTest] = C[nTest] + c[nt]
    end
  end

  return C
end
