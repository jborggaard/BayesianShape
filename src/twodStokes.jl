function twodStokes(N::Int)
#  Solves Stokes equation in 2D with Dirichlet boundary conditions
#     - ∇⋅(∇z+∇z') + ∇p = f,  Ω = (0,1)×(0,1),  homogeneous Dirichlet b.c.
#---------------------------------------------------------------------------78--

  include("twodMesh.jl")
  include("twodQuadratureRule.jl")
#  @everywhere include("twodShape.jl")
#  @everywhere include("twodBilinear.jl")
#  @everywhere include("twodLinForm.jl")
  include("twodShape.jl")
  include("twodBilinear.jl")
  include("twodLinForm.jl")

  #  Define problem parameters
#  @everywhere μ = 0.001;
#  @everywhere ϵ = 0.0;
  μ = 0.001;
  ϵ = 0.0;

  xMin =-1.0;
  xMax = 1.0;
  yMin =-1.0;
  yMax = 1.0;

  rule  = 3;   # points in quadrature formula


  #  Specify manufactured solutions for code verification
  #-----------------------------------------------------------------------------
  function zExact(x)
    return sin(π*x[1])*sin(2.0*π*x[2]);
  end

#  @everywhere function fx(x::Array{Float64,2})
  function fx(x::Array{Float64,2})
    return 1.0 - 6.0*μ*x[2]; #1.0 + 8.0*μ*x[2];
  end
#  @everywhere function fy(x::Array{Float64,2})
  function fy(x::Array{Float64,2})
    return 1.0 + 6.0*μ*x[1]; #4.0*μ*(1.0-2.0*x[1]);
  end

  #  Create mesh and set equation numbers and boundary conditions
  #-----------------------------------------------------------------------------
  nNodesX = N; nNodesY = N;
  x,eConn,iB = twodMesh( xMin,xMax, yMin,yMax, "quadratic", nNodesX,nNodesY );
  nNodes = size(x,1);
  nElements = size(eConn,1);

  #  Set the index into equation numbers
  #-----------------------------------------------------------------------------
  ideU = zeros(Int64,nNodes,2);
  ideP = zeros(Int64,nNodes,1);

  nBC = 2*size(iB,1);
  dirU = zeros(Float64,nBC);

  global n_diru = 0;
#=  for i=iB
    if ( x[i,2]>yMax-1e-8 )
      n_diru = n_diru + 1;
      ideU[i,1] = -n_diru;
      dirU[n_diru] = 4.0*(x[i,1]-xMin)*(xMax-x[i,1])/(xMax-xMin)^2;

      n_diru = n_diru + 1;
      ideU[i,2] = -n_diru;
      dirU[n_diru] = -2*(x[i,2]-yMin)^2/(yMax-yMin) * (xMax+xMin-2*x[i,1])/(xMax-xMin);
    else
      n_diru = n_diru + 1;
      ideU[i,1] = -n_diru;
      dirU[n_diru] = 0.0;

      n_diru = n_diru + 1;
      ideU[i,2] = -n_diru;
      dirU[n_diru] = -2*(x[i,2]-yMin)^2/(yMax-yMin) * (xMax+xMin-2*x[i,1])/(xMax-xMin);
    end
end =#

  ideU = zeros(Int64,nNodes,2);  # rezero (will find a way to identify DBC later)

  global nUnk = 0;
  iU = zeros(Int64,2*nNodes); #iU = zeros(Int64,2*nNodes-nBC);
  for i=1:nNodes
    if ideU[i,1]==0
      nUnk = nUnk+1;
      iU[nUnk] = i;
      ideU[i,1] = nUnk;
    end
    if ideU[i,2]==0
      nUnk = nUnk+1;
      iU[nUnk] = i;
      ideU[i,2] = nUnk;
    end
  end

  ideP = zeros(Int64,nNodes);
  for n_el=1:nElements
    vertex = eConn[n_el,1:3];
    if ( ideP[vertex[1]]==0 )
      nUnk = nUnk + 1;
      ideP[vertex[1]] = nUnk;
    end

    if ( ideP[vertex[2]]==0 )
      nUnk = nUnk + 1;
      ideP[vertex[2]] = nUnk;
    end

    if ( ideP[vertex[3]]==0 )
      nUnk = nUnk + 1;
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
  b   = Array{Float64,1}(undef, 2*nNodes);

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

  #print("$(A)")
  #print("$(b)")

#=  u = zeros(Float64,nNodes);
  u[iU] = A[iU,iU]\b[iU];

  uExact = zeros(Float64,nNodes);
  for n=1:nNodes
    uExact[n] = zExact(x[n,:])
  end

  #print("$(u)")
=#
#  return x,eConn,u,uExact
  return A,b
end
