function twodNavierStokesRotatingOuterNewton(x,eConn,innerNodes,outerNodes,ω,velocity0,pressure0,μ)


ϵ = 0.001

rule  = 3   # points in quadrature formula

function fx(x::Array{Float64,2})
    return 0.0 #1.0 - 6.0*μ*x[2] #1.0 + 8.0*μ*x[2]
end
#  @everywhere function fy(x::Array{Float64,2})
function fy(x::Array{Float64,2})
    return 0.0 #1.0 + 6.0*μ*x[1] #4.0*μ*(1.0-2.0*x[1])
end

  #  Create equation numbers and set boundary conditions
  #-----------------------------------------------------------------------------
  nNodes = size(x,1)
  nElements = size(eConn,1)

  nInnerNodes = length(innerNodes)
  nOuterNodes = length(outerNodes)
  nDirichlet = nInnerNodes+nOuterNodes

  #  Set the index into equation numbers
  #-----------------------------------------------------------------------------
  ideU = zeros(Int64,nNodes,2)

  global nUnk = 0
  for i=1:nNodes
    if ideU[i,1]==0
      global nUnk = nUnk+1
      ideU[i,1] = nUnk
    end
    if ideU[i,2]==0
      global nUnk = nUnk+1
      ideU[i,2] = nUnk
    end
  end

  ideP = zeros(Int64,nNodes)
  for n_el=1:nElements
    vertex = eConn[n_el,1:3]
    if ( ideP[vertex[1]]==0 )
      global nUnk = nUnk + 1
      ideP[vertex[1]] = nUnk
    end

    if ( ideP[vertex[2]]==0 )
      global nUnk = nUnk + 1
      ideP[vertex[2]] = nUnk
    end

    if ( ideP[vertex[3]]==0 )
      global nUnk = nUnk + 1
      ideP[vertex[3]] = nUnk
    end
  end
  nP = nUnk - 2*nNodes

  #  Integrate system matrices, element-by-element
  #-----------------------------------------------------------------------------
  r,s,w = twodQuadratureRule(rule)
  μ_g   = μ*ones(rule)
  ϵ_g   = ϵ*ones(rule)
  one   = ones(rule)

  nElDOF   = size(eConn,2)
  nElDOF2  = nElDOF*nElDOF
  nEntries = nElements*(2*nElDOF+3)*(2*nElDOF+3)

  II  = Array{Int64,1}(undef, nEntries)
  JJ  = Array{Int64,1}(undef, nEntries)
  AA  = Array{Float64,1}(undef, nEntries)
  b   = zeros(Float64,nUnk,1)

for k=1:nElements


    nLocal = eConn[k,:][:]
    xLocal = x[nLocal,:]

    xg,wg,phi,phi_x,phi_y = twodShape( xLocal, r, s, w )
    xg,wg,psi,psi_x,psi_y = twodShape( xLocal[1:3,:], r, s, w )

    # evaluate forcing function at quadrature points
    fx_g   = fx(xg)
    fy_g   = fy(xg)


    # u_n, v_n and p_n are the previous solution coeffecients at the element level. 

    u_n   = phi*velocity0[nLocal,1]
    u_n_x = phi_x*velocity0[nLocal,1]
    u_n_y = phi_y*velocity0[nLocal,1]

    v_n   = phi*velocity0[nLocal,2]
    v_n_x = phi_x*velocity0[nLocal,2]
    v_n_y = phi_y*velocity0[nLocal,2]

    p_n   = psi*pressure0[nLocal[1:3]]

    A11Loc =-twodBilinear( 2*μ_g, phi_x, phi_x, wg) - twodBilinear(   μ_g, phi_y, phi_y, wg) 
    A12Loc =-twodBilinear(   μ_g, phi_x, phi_y, wg) 
    A21Loc =-twodBilinear(   μ_g, phi_y, phi_x, wg) 
    A22Loc =-twodBilinear(   μ_g, phi_x, phi_x, wg) - twodBilinear( 2*μ_g, phi_y, phi_y, wg) 

    W_11_Loc =-twodBilinear(u_n_x,phi,phi,wg)
    W_12_Loc =-twodBilinear(u_n_y,phi,phi,wg)
    W_21_Loc =-twodBilinear(v_n_x,phi,phi,wg) 
    W_22_Loc =-twodBilinear(v_n_y,phi,phi,wg)

    N_Loc  =-twodBilinear(  u_n, phi_x, phi, wg) - twodBilinear(v_n, phi_y, phi, wg) 

    A11Loc = A11Loc + W_11_Loc +N_Loc
    A12Loc = A12Loc + W_12_Loc
    A21Loc = A21Loc + W_21_Loc
    A22Loc = A22Loc + W_22_Loc +N_Loc

    B1Loc  = twodBilinear(   one, phi_x, psi  , wg) 
    B2Loc  = twodBilinear(   one, phi_y, psi  , wg) 

    MLoc   = twodBilinear(   ϵ_g, psi  , psi  , wg) 
  
    F1Loc  = twodLinForm(   fx_g, phi  ,        wg) 
    F2Loc  = twodLinForm(   fy_g, phi  ,        wg) 

    F_NL_x_Loc = twodLinForm(u_n.*u_n_x, phi, wg)  + twodLinForm(v_n.*u_n_y, phi, wg)
    F_NL_y_Loc = twodLinForm(u_n.*v_n_x, phi, wg)  + twodLinForm(v_n.*v_n_y, phi, wg)

    F_L_x_Loc  = twodLinForm(2*μ*u_n_x, phi_x, wg)     + twodLinForm(μ*u_n_y, phi_y, wg)       + twodLinForm(μ*v_n_x, phi_y, wg)
    F_L_y_Loc  = twodLinForm(2*μ*v_n_y, phi_y, wg)     + twodLinForm(μ*u_n_y, phi_x, wg)       + twodLinForm(μ*v_n_x, phi_x, wg)

    F_p_x_Loc  = twodLinForm(p_n, phi_x, wg)
    F_p_y_Loc  = twodLinForm(p_n, phi_y, wg)

    F1Loc = F1Loc-F_NL_x_Loc-F_L_x_Loc+F_p_x_Loc

    F2Loc = F2Loc-F_NL_y_Loc-F_L_y_Loc+F_p_y_Loc


    g_Loc =     twodLinForm(u_n_x+v_n_y, psi, wg) + ϵ*twodLinForm(p_n, psi, wg)
                

    
    index = (k-1)*(2*nElDOF+3)^2  # compute base index (since k is distributed)
    lDOF = ideU[nLocal,1][:]
    # u-momentum equations
    for nt = 1:nElDOF
      nTest = ideU[nLocal[nt],1]
      for nu = 1:nElDOF
        nUnkU = ideU[nLocal[nu],1]
        nUnkV = ideU[nLocal[nu],2]

        index = index + 1
        II[index] = nTest
        JJ[index] = nUnkU
        AA[index] = A11Loc[nt,nu]

        index = index + 1
        II[index] = nTest
        JJ[index] = nUnkV
        AA[index] = A12Loc[nt,nu]
      end
      for np = 1:3  # special to this linear element
        nUnkP = ideP[nLocal[np]]

        index = index + 1
        II[index] = nTest
        JJ[index] = nUnkP
        AA[index] = B1Loc[np,nt]
      end

      b[nTest] = b[nTest] + F1Loc[nt]
    end

    # v-momentum equations
    for nt = 1:nElDOF
      nTest = ideU[nLocal[nt],2]
      for nu = 1:nElDOF
        nUnkU = ideU[nLocal[nu],1]
        nUnkV = ideU[nLocal[nu],2]

        index = index + 1
        II[index] = nTest
        JJ[index] = nUnkU
        AA[index] = A21Loc[nt,nu]

        index = index + 1
        II[index] = nTest
        JJ[index] = nUnkV
        AA[index] = A22Loc[nt,nu]
      end
      for np = 1:3  # special to this linear element
        nUnkP = ideP[nLocal[np]]

        index = index + 1
        II[index] = nTest
        JJ[index] = nUnkP
        AA[index] = B2Loc[np,nt]
      end

      b[nTest] = b[nTest] + F2Loc[nt]
    end

    # continuity equations
    for nt = 1:3
      nTest = ideP[nLocal[nt]]
      for nu = 1:nElDOF
        nUnkU = ideU[nLocal[nu],1]
        nUnkV = ideU[nLocal[nu],2]

        index = index + 1
        II[index] = nTest
        JJ[index] = nUnkU
        AA[index] = B1Loc[nt,nu]

        index = index + 1
        II[index] = nTest
        JJ[index] = nUnkV
        AA[index] = B2Loc[nt,nu]
      end
      for np = 1:3  # special to this linear element
        nUnkP = ideP[nLocal[np]]

        index = index + 1
        II[index] = nTest
        JJ[index] = nUnkP
        AA[index] = MLoc[np,nt]
      end
      b[nTest] = b[nTest] + g_Loc[nt]  #This is new 
    end

  end
  A = sparse(II,JJ,AA)

  #  Apply the boundary conditions and solve
  #-------------------------------------------------------------------------78--

  interiorNodes = 1:nNodes
  interiorNodes = setdiff(interiorNodes,innerNodes)
  interiorNodes = setdiff(interiorNodes,outerNodes)

  knownIndexU = [2*innerNodes-ones(UInt64,nInnerNodes,1) 
                 2*outerNodes-ones(UInt64,nOuterNodes,1)]
  knownIndexV = [2*innerNodes 
                 2*outerNodes]

  dirichletU = Array{Float64,1}(undef,nDirichlet)
  dirichletV = Array{Float64,1}(undef,nDirichlet)
  for i=1:nInnerNodes
    dirichletU[i] = 0.0
    dirichletV[i] = 0.0
  end
  for i=(nInnerNodes+1):(nInnerNodes+nOuterNodes)
    dirichletU[i] = 0.0
    dirichletV[i] = 0.0
  end

  nUnknowns = 2*length(interiorNodes) + nP
  unknownIndex = Array{UInt64,1}(undef,nUnknowns)
  uIndex = Array{UInt64,1}(undef,2*length(interiorNodes))
  vIndex = Array{UInt64,1}(undef,2*length(interiorNodes))
  uInverseIndex = Array{UInt64,1}(undef,nNodes)
  vInverseIndex = Array{UInt64,1}(undef,nNodes)

  global nUnknownCounter = 0
  for i ∈ interiorNodes
    global nUnknownCounter = nUnknownCounter+1
    unknownIndex[nUnknownCounter] = 2*i-1
    uIndex[nUnknownCounter] = i
    uInverseIndex[i]=nUnknownCounter

    global nUnknownCounter = nUnknownCounter+1
    unknownIndex[nUnknownCounter] = 2*i
    vIndex[nUnknownCounter] = i
    vInverseIndex[i]=nUnknownCounter
  end

  getP = zeros(Int64,3*nNodes,1)
  for nNode=1:nNodes
    if ideP[nNode]>0
      global nUnknownCounter = nUnknownCounter + 1
      unknownIndex[nUnknownCounter] = ideP[nNode]
      getP[ideP[nNode]] = nUnknownCounter
    end
  end

  rhs = -b[unknownIndex]-A[unknownIndex,knownIndexU[:]]*dirichletU-A[unknownIndex,knownIndexV[:]]*dirichletV

  uvp = A[unknownIndex,unknownIndex]\rhs

  dvelocity = zeros(Float64,nNodes,2)

  # here is where we will eventually apply the "inverse" index mappings
  global nUnknownCounter = 0
  for i=interiorNodes
    global nUnknownCounter = nUnknownCounter+1
    dvelocity[i,1] = uvp[nUnknownCounter]

    global nUnknownCounter = nUnknownCounter + 1
    dvelocity[i,2] = uvp[nUnknownCounter]
  end

  for i=(nInnerNodes+1):(nInnerNodes+nOuterNodes)
    dvelocity[outerNodes[i-nInnerNodes],1] = dirichletU[i]
    dvelocity[outerNodes[i-nInnerNodes],2] = dirichletV[i]
  end

  # to be computed in the future if needed
  dpressure = zeros(Float64,nNodes,1)

  vv = eConn[:,1:3]
  vvv = vv[:]
  vertices = unique(sort(vvv))
  for i ∈ vertices
    dpressure[i] = uvp[getP[ideP[i]]]
  end
  dpressure = TriMesh_PromoteL2Q(dpressure,eConn)

  



  return dvelocity,dpressure
  
end
