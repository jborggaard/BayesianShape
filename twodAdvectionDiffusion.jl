function twodAdvectionDiffusion(x,eConn,innerNodes,outerNodes,velocity)
#
#  Solves the advection-diffusion equation in 2D
#     - ∇⋅(κ∇θ) + v⋅∇θ = q,
#  Ω is the domain between a Bspline disk and an
#  outer circle.  We prescribe zero Dirichlet conditions to the outer boundary.
#---------------------------------------------------------------------------78--

#  include("twodQuadratureRule.jl")
#  include("twodShape.jl")
#  include("twodBilinear.jl")
#  include("twodLinForm.jl")

  #  Define problem parameters
  κ = 1.0;

  rule  = 7;   # points in quadrature formula

  function q(x::Array{Float64,2})
    dist2 = (x[:,1].-1.5).^2+(x[:,2].-0.75).^2;
    return exp.(-dist2);
    #d = size(x,1);
    #return zeros(Float64,d,1)
  end

  #  Get problem dimensions
  #-----------------------------------------------------------------------------
  nNodes = size(x,1);
  nElements = size(eConn,1);

  nInnerNodes = length(innerNodes);
  nOuterNodes = length(outerNodes);
#  nDirichlet = nInnerNodes+nOuterNodes;   # use for Dirichlet BCs
  nDirichlet = nOuterNodes;

  #  Set the index into equation numbers
  #-----------------------------------------------------------------------------
  ideθ = zeros(Int64,nNodes,1);

  global nUnk = 0;
  for i=1:nNodes
    global nUnk = nUnk+1;
    ideθ[i,1] = nUnk;
  end

  #  Integrate system matrices, element-by-element
  #-----------------------------------------------------------------------------
  r,s,w = twodQuadratureRule(rule);
  κ_g   = κ*ones(rule);
  one   = ones(rule);

  nElDOF   = size(eConn,2);
  nElDOFsq = nElDOF*nElDOF;
  nEntries = nElements*nElDOFsq;
#  II  = SharedArray{Int64}(nEntries);
#  JJ  = SharedArray{Int64}(nEntries);
#  AA  = SharedArray{Float64}(nEntries);
#  b   = SharedArray{Float64}(nNodes);
  II  = Array{Int64,1}(undef, nEntries);
  JJ  = Array{Int64,1}(undef, nEntries);
  AA  = Array{Float64,1}(undef, nEntries);
  b   = zeros(Float64,nNodes,1);

#  @sync @distributed for k=1:nElements
  for k=1:nElements
    xg  = zeros(Float64,rule,2);
    wg  = zeros(Float64,rule);
    ϕ = zeros(Float64,rule,nElDOF);
    ϕ_x = zeros(Float64,rule,nElDOF);
    ϕ_y = zeros(Float64,rule,nElDOF);

    nLocal = eConn[k,:][:];
    xLocal = x[nLocal,:];

    xg,wg,ϕ,ϕ_x,ϕ_y = twodShape( xLocal, r, s, w );

    q_g   = q(xg);           # forcing function evaluated at quadrature points
    u_g   = ϕ*velocity[nLocal,1];
    v_g   = ϕ*velocity[nLocal,2];

    ALoc = twodBilinear( u_g, ϕ_x, ϕ  , wg) + twodBilinear( v_g, ϕ_y, ϕ  , wg) +
           twodBilinear( κ_g, ϕ_x, ϕ_x, wg) + twodBilinear( κ_g, ϕ_y, ϕ_y, wg) ;
    FLoc = twodLinForm(  q_g, ϕ  ,      wg) ;

    index = (k-1)*nElDOFsq;  # compute base index (since k could be in any order)
    lDOF = ideθ[nLocal,1][:];

    for nt = 1:nElDOF
      nTest = ideθ[nLocal[nt],1];
      for nu = 1:nElDOF
        nUnkT = ideθ[nLocal[nu],1];

        index = index + 1;
        II[index] = nTest;
        JJ[index] = nUnkT;
        AA[index] = ALoc[nt,nu];
      end

      b[nTest] = b[nTest] + FLoc[nt];
    end
  end

  A = sparse(II,JJ,AA);

  interiorNodes = 1:nNodes;
#  interiorNodes = setdiff(interiorNodes,innerNodes);
  interiorNodes = setdiff(interiorNodes,outerNodes);

#  knownIndexθ = [innerNodes;outerNodes];
  knownIndexθ = outerNodes;

  dirichletθ = Array{Float64,1}(undef,nDirichlet);
  for i=1:nDirichlet
    nO = knownIndexθ[i];
    dirichletθ[i] = 0.0;
#    dirichletθ[i] = x[nO,1];
  end

  nUnknowns = length(interiorNodes);
  unKnownIndex = Array{UInt64,1}(undef,nUnknowns)

  global nUnknownCounter = 0;
  for i ∈ interiorNodes
    global nUnknownCounter = nUnknownCounter+1;
    unKnownIndex[nUnknownCounter] = i;
  end

  θunk = A[unKnownIndex,unKnownIndex]\(b[unKnownIndex]-A[unKnownIndex,knownIndexθ[:]]*dirichletθ)

  θ = zeros(nNodes,1)

  global nUnknownCounter = 0;
  for i ∈ interiorNodes
    global nUnknownCounter = nUnknownCounter+1;
    θ[i] = θunk[nUnknownCounter];
  end

#  for i=1:nInnerNodes
#    θ[innerNodes[i]] = dirichletθ[i];
#  end

#  for i=1:nOuterNodes
#    θ[outerNodes[i]] = dirichletθ[nInnerNodes+i];
#  end
  for i=1:nOuterNodes
    θ[outerNodes[i]] = dirichletθ[i];
  end

  return θ
end
