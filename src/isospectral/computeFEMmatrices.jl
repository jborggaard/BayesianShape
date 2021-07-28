function computeFEMmatrices(x,eConn,
                            dirichletNodes,dirichletValues;
                            κ::Float64=1.0)  # optional argument
#
#  Computes the finite element matrices (mass and stiffness) for the problem
#  of hearing the shape of a drum.
#---------------------------------------------------------------------------78

  rule = 7   # integrates a 5th degree polynomial exactly

  nNodes = size(x,1)
  nElements = size(eConn,1)
  nDirichlet = length(dirichletNodes)

  #  Set the index into equation numbers (known and unknowns)
  #---------------------------------------------------------------------------
  ide = zeros(Int64,nNodes,1)
  
  global nUnk = 0
  for i=1:nNodes
    global nUnk = nUnk+1
    ide[i,1] = nUnk
  end

  # ide = 1:nNodes

  #  Integrate system matrices, element-by-element
  #---------------------------------------------------------------------------
  r,s,w = twodQuadratureRule(rule)
  κ_g   = κ*ones(rule)
  o_g   = ones(rule)

  nElDOF   = size(eConn,2)
  nElDOF2  = nElDOF*nElDOF
  nEntries = nElements*nElDOF2
  II  = Array{Int64,1}(undef, nEntries)    #SharedArray(Int32,nEntries);
  JJ  = Array{Int64,1}(undef, nEntries)    #SharedArray(Int32,nEntries);
  AA  = Array{Float64,1}(undef, nEntries)  #SharedArray(Float64,nEntries);
  MM  = Array{Float64,1}(undef, nEntries)  #SharedArray(Float64,nEntries);

#  @sync @parallel for k=1:nElements
  for k=1:nElements
    nLocal = eConn[k,:][:]
    xLocal = x[nLocal,:]

    xg,wg,ϕ,ϕ_x,ϕ_y = twodShape( xLocal, r, s, w )

    ALoc = twodBilinear( κ_g, ϕ_x, ϕ_x, wg ) + twodBilinear( κ_g, ϕ_y, ϕ_y, wg )
    MLoc = twodBilinear( o_g, ϕ  , ϕ  , wg )

    index = (k-1)*nElDOF2    # compute base index (k could be in any order)

    for nt = 1:nElDOF
      nTest = ide[nLocal[nt],1]
      for nu = 1:nElDOF
        nUnkU = ide[nLocal[nu],1]

        index = index + 1
        II[index] = nTest
        JJ[index] = nUnkU
        AA[index] = ALoc[nt,nu]
        MM[index] = MLoc[nt,nu]
      end
    end
  end

  A = sparse(II,JJ,AA)
  M = sparse(II,JJ,MM)

  interiorNodes = convert(Array{Int64,1},1:nNodes)
  setdiff!(interiorNodes,dirichletNodes)

  knownIndex = dirichletNodes

  unknownIndex = interiorNodes

  Amat = A[unknownIndex,unknownIndex]
  Mmat = M[unknownIndex,unknownIndex]

  return Amat,Mmat
end
