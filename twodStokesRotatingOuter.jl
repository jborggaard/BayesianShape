function twodStokesRotatingOuter(x,eConn,innerNodes,outerNodes,ω)
#
#  Solves the Stokes equation in 2D with Dirichlet boundary conditions
#     - ∇⋅(∇z+∇z') + ∇p = f,  Ω is the domain between a Bspline disk and an
#  outer circle.  We prescribe a rotation to the outer surface.
#---------------------------------------------------------------------------78--

#  include("twodQuadratureRule.jl")
#  @everywhere include("twodShape.jl")
#  @everywhere include("twodBilinear.jl")
#  @everywhere include("twodLinForm.jl")
#  include("twodShape.jl")
#  include("twodBilinear.jl")
#  include("twodLinForm.jl")

  #  Define problem parameters
#  @everywhere μ = 0.001;
#  @everywhere ϵ = 0.001;
  μ = 0.001;
  ϵ = 0.001;

  rule  = 3;   # points in quadrature formula

  function fx(x::Array{Float64,2})
    return 0.0; #1.0 - 6.0*μ*x[2]; #1.0 + 8.0*μ*x[2];
  end
#  @everywhere function fy(x::Array{Float64,2})
  function fy(x::Array{Float64,2})
    return 0.0 #1.0 + 6.0*μ*x[1]; #4.0*μ*(1.0-2.0*x[1]);
  end

  #  Create equation numbers and set boundary conditions
  #-----------------------------------------------------------------------------
  nNodes = size(x,1);
  nElements = size(eConn,1);

  nInnerNodes = length(innerNodes);
  nOuterNodes = length(outerNodes);
  nDirichlet = nInnerNodes+nOuterNodes;

  #  Set the index into equation numbers
  #-----------------------------------------------------------------------------
  ideU = zeros(Int64,nNodes,2);

  global nUnk = 0;
  for i=1:nNodes
    if ideU[i,1]==0
      global nUnk = nUnk+1;
      ideU[i,1] = nUnk;
    end
    if ideU[i,2]==0
      global nUnk = nUnk+1;
      ideU[i,2] = nUnk;
    end
  end

  ideP = zeros(Int64,nNodes);
  for n_el=1:nElements
    vertex = eConn[n_el,1:3];
    if ( ideP[vertex[1]]==0 )
      global nUnk = nUnk + 1;
      ideP[vertex[1]] = nUnk;
    end

    if ( ideP[vertex[2]]==0 )
      global nUnk = nUnk + 1;
      ideP[vertex[2]] = nUnk;
    end

    if ( ideP[vertex[3]]==0 )
      global nUnk = nUnk + 1;
      ideP[vertex[3]] = nUnk;
    end
  end
  nP = nUnk - 2*nNodes;

  #  Integrate system matrices, element-by-element
  #-----------------------------------------------------------------------------
  r,s,w = twodQuadratureRule(rule);
  μ_g   = μ*ones(rule);
  ϵ_g   = ϵ*ones(rule);
  one   = ones(rule);

  nElDOF   = size(eConn,2);
  nElDOF2  = nElDOF*nElDOF;
  nEntries = nElements*(2*nElDOF+3)*(2*nElDOF+3);
#  II  = SharedArray{Int64}(nEntries);
#  JJ  = SharedArray{Int64}(nEntries);
#  AA  = SharedArray{Float64}(nEntries);
#  b   = SharedArray{Float64}(2*nNodes);
  II  = Array{Int64,1}(undef, nEntries);
  JJ  = Array{Int64,1}(undef, nEntries);
  AA  = Array{Float64,1}(undef, nEntries);
  b   = zeros(2*nNodes,1);

#  @sync @distributed for k=1:nElements
  for k=1:nElements
    xg  = zeros(Float64,rule,2);
    wg  = zeros(Float64,rule);
    phi = zeros(Float64,rule,nElDOF);
    p_x = zeros(Float64,rule,nElDOF);
    p_y = zeros(Float64,rule,nElDOF);

    nLocal = eConn[k,:][:];
    xLocal = x[nLocal,:];

    xg,wg,phi,phi_x,phi_y = twodShape( xLocal, r, s, w );
    xg,wg,psi,psi_x,psi_y = twodShape( xLocal[1:3,:], r, s, w );

    fx_g   = fx(xg);           # forcing function evaluated at quadrature points
    fy_g   = fy(xg);           # forcing function evaluated at quadrature points

    A11Loc =-twodBilinear( 2*μ_g, phi_x, phi_x, wg) -
             twodBilinear(   μ_g, phi_y, phi_y, wg) ;
    A12Loc =-twodBilinear(   μ_g, phi_x, phi_y, wg) ;
    A21Loc =-twodBilinear(   μ_g, phi_y, phi_x, wg) ;
    A22Loc =-twodBilinear(   μ_g, phi_x, phi_x, wg) -
             twodBilinear( 2*μ_g, phi_y, phi_y, wg) ;
    B1Loc  = twodBilinear(   one, phi_x, psi  , wg) ;
    B2Loc  = twodBilinear(   one, phi_y, psi  , wg) ;
    MLoc   =-twodBilinear(   ϵ_g, psi  , psi  , wg) ;
    F1Loc  = twodLinForm(  fx_g, phi  ,        wg) ;
    F2Loc  = twodLinForm(  fy_g, phi  ,        wg) ;

    index = (k-1)*(2*nElDOF+3)^2;  # compute base index (since k could be in any order)
    #@printf("%g:\n",index)
    lDOF = ideU[nLocal,1][:];
  #  @printf("%g %g %g %g %g %g\n",lDOF[1],lDOF[2],lDOF[3],lDOF[4],lDOF[5],lDOF[6])
    # u-momentum equations
    for nt = 1:nElDOF
      nTest = ideU[nLocal[nt],1];
#      @printf("%g\n",nTest)
      for nu = 1:nElDOF
        nUnkU = ideU[nLocal[nu],1];
        nUnkV = ideU[nLocal[nu],2];

        index = index + 1;
        II[index] = nTest;
        JJ[index] = nUnkU;
        AA[index] = A11Loc[nt,nu];

        index = index + 1;
        II[index] = nTest;
        JJ[index] = nUnkV;
        AA[index] = A12Loc[nt,nu];
      end
      for np = 1:3  # special to this linear element
        nUnkP = ideP[nLocal[np]];

        index = index + 1;
        II[index] = nTest;
        JJ[index] = nUnkP;
        AA[index] = B1Loc[np,nt];
      end

      b[nTest] = b[nTest] + F1Loc[nt];
    end

    # v-momentum equations
    for nt = 1:nElDOF
      nTest = ideU[nLocal[nt],2];
      for nu = 1:nElDOF
        nUnkU = ideU[nLocal[nu],1];
        nUnkV = ideU[nLocal[nu],2];

        index = index + 1;
        II[index] = nTest;
        JJ[index] = nUnkU;
        AA[index] = A21Loc[nt,nu];

        index = index + 1;
        II[index] = nTest;
        JJ[index] = nUnkV;
        AA[index] = A22Loc[nt,nu];
      end
      for np = 1:3  # special to this linear element
        nUnkP = ideP[nLocal[np]];

        index = index + 1;
        II[index] = nTest;
        JJ[index] = nUnkP;
        AA[index] = B2Loc[np,nt];
      end

      b[nTest] = b[nTest] + F2Loc[nt];
    end

    # continuity equations
    for nt = 1:3
      nTest = ideP[nLocal[nt]];
      for nu = 1:nElDOF
        nUnkU = ideU[nLocal[nu],1];
        nUnkV = ideU[nLocal[nu],2];

        index = index + 1;
        II[index] = nTest;
        JJ[index] = nUnkU;
        AA[index] = B1Loc[nt,nu];

        index = index + 1;
        II[index] = nTest;
        JJ[index] = nUnkV;
        AA[index] = B2Loc[nt,nu];
      end
      for np = 1:3  # special to this linear element
        nUnkP = ideP[nLocal[np]];

        index = index + 1;
        II[index] = nTest;
        JJ[index] = nUnkP;
        AA[index] = MLoc[np,nt];
      end

    end
    #@printf("%g\n",index)
  end
  #print("$(II)")
  A = sparse(II,JJ,AA);

  interiorNodes = 1:nNodes;
  interiorNodes = setdiff(interiorNodes,innerNodes);
  interiorNodes = setdiff(interiorNodes,outerNodes);

  knownIndexU = [2*innerNodes-ones(UInt64,nInnerNodes,1); 2*outerNodes-ones(UInt64,nOuterNodes,1)];
  knownIndexV = [2*innerNodes; 2*outerNodes];

  dirichletU = Array{Float64,1}(undef,nDirichlet);
  dirichletV = Array{Float64,1}(undef,nDirichlet);
  for i=1:nInnerNodes
    dirichletU[i] = 0.0;#-ω*x[innerNodes[i],2];
    dirichletV[i] = 0.0;#ω*x[innerNodes[i],1];
  end
  for i=(nInnerNodes+1):(nInnerNodes+nOuterNodes)
    dirichletU[i] = ω*x[outerNodes[i-nInnerNodes],2];
    dirichletV[i] = -ω*x[outerNodes[i-nInnerNodes],1];
  end

  nUnknowns = length(interiorNodes);
  unKnownIndex = Array{UInt64,1}(undef,2*nUnknowns)
  uIndex = Array{UInt64,1}(undef,2*nUnknowns)
  vIndex = Array{UInt64,1}(undef,2*nUnknowns)
  uInverseIndex = Array{UInt64,1}(undef,nNodes)
  vInverseIndex = Array{UInt64,1}(undef,nNodes)

  global nUnknownCounter = 0;
  for i=interiorNodes
    global nUnknownCounter = nUnknownCounter+1;
    unKnownIndex[nUnknownCounter] = 2*i-1;
    uIndex[nUnknownCounter] = i;
    uInverseIndex[i]=nUnknownCounter;

    global nUnknownCounter = nUnknownCounter+1;
    unKnownIndex[nUnknownCounter] = 2*i;
    vIndex[nUnknownCounter] = i;
    vInverseIndex[i]=nUnknownCounter;
  end

  uvp = A[unKnownIndex,unKnownIndex]\(b[unKnownIndex]-A[unKnownIndex,knownIndexU[:]]*dirichletU-A[unKnownIndex,knownIndexV[:]]*dirichletV)

  velocity = zeros(nNodes,2)

  nUnknownCounter = 0;
  for i=interiorNodes
    global nUnknownCounter = nUnknownCounter+1;
    velocity[i,1] = uvp[nUnknownCounter];

    nUnknownCounter = nUnknownCounter + 1;
    velocity[i,2] = uvp[nUnknownCounter];
  end

  for i=(nInnerNodes+1):(nInnerNodes+nOuterNodes)
    velocity[outerNodes[i-nInnerNodes],1] = dirichletU[i-nInnerNodes];
    velocity[outerNodes[i-nInnerNodes],2] = dirichletV[i-nInnerNodes];
  end

  pressure = zeros(Float64,nNodes,1);
  #for i=(nInnerNodes+1):(nInnerNodes+nOuterNodes)
  #  velocity[outerNodes[i-nInnerNodes],1] = 0;
  #  velocity[outerNodes[i-nInnerNodes],2] = 0;
  #end
#  nUKnowns = length(nInner)+length(nOuter)
#  nVKnowns = length(nInner)+length(nOuter)

#  ideU = zeros(Int64,nNodes,2);
#  ideP = zeros(Int64,nNodes,1);
#
#  nBC = nUKnowns + nVKnowns;
#  dirU = zeros(Float64,nBC);
#
#  for n=nInner
#      ideU(n,1) = -1;
#      ideU(n,2) = -1;
#  end
#  for n=nOuter
#      ideU(n,1) = -2;
#      ideU(n,2) = -2;
#  end
#
#  idUnk = 0;
#  idKnn = 0;
#  for n=1:nNodes
#    if ( ideU(n,1)==0 )
#      idx = idx+1;
#      ideU(n,1) = idx;
#    else
#      idKnn = idKnn+1;
#      dirU(idKnn) = -omega*x[i,2]
#    end
#
#    if ( ideU(n,2)==0 )
#      idx = idx+1;
#      ideU(n,2) = idx;
#    else
#      idKnn = idKnn+1;
#      dirU(idKnn) = omega*x[i,1]
#    end
#  end



  #print("$(A)")
  #print("$(b)")

#=  u = zeros(Float64,nNodes);

  uExact = zeros(Float64,nNodes);
  for n=1:nNodes
    uExact[n] = zExact(x[n,:])
  end

  #print("$(u)")
=#
#  return x,eConn,u,uExact

  return velocity, pressure
end
