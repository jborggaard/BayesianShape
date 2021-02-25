function twodMassMatrix(x,eConn)
#  Creates the 2D finite element mass matrix (Natural boundary conditions
#
#  In the future, this should have an optional argument (ide) that allows
#  for various boundary condition types.

  nNodes    = size(x    ,1);
  nElements = size(eConn,1);

  #  Set the index into equation numbers
  ide = zeros(Int64,nNodes);

#  #  If there are no Dirichlet boundary conditions, and no ide array...
#  nUnk = 0;
#  iU = zeros(Int64,nNodes);
#  for i=1:nNodes
#    if ide[i]==0
#      nUnk = nUnk+1;
#      iU[nUnk] = i;
#      ide[i] = nUnk;
#    end
#  end

  ide = 1:nNodes;

  rule  = 3;
  r,s,w = twodQuadratureRule(rule);
  o     = ones(rule);

  nElDOF = size(eConn,2);
  nEntries = nElements*nElDOF^2;
  II  = zeros(Int64,nEntries);
  JJ  = zeros(Int64,nEntries);
  MM  = zeros(Float64,nEntries);

  xg  = zeros(Float64,rule,2);
  wg  = zeros(Float64,rule);
  phi = zeros(Float64,rule,nElDOF);
  p_x = zeros(Float64,rule,nElDOF);
  p_y = zeros(Float64,rule,nElDOF);

  index = 0;
  for k=1:nElements

    nLocal = eConn[k,:][:];
    xLocal = x[nLocal,:];

    xg,wg,phi,p_x,p_y = twodShape( xLocal, r, s, w )

    MLoc = twodBilinear(  o, phi, phi, wg )

    lDOF = ide[nLocal];
    for nt = 1:nElDOF
      for nu = 1:nElDOF
        index = index + 1;
        II[index] = lDOF[nu];
        JJ[index] = lDOF[nt];
        MM[index] = MLoc[nt,nu];
      end
    end
  end

  M = sparse(II,JJ,MM);

  return M
end
