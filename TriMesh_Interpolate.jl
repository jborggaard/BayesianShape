function TriMesh_Interpolate(x,eConn,eAdjacency,fNodes,xInterp; verbose=0, startEl=rand(1:size(eConn,1),size(xInterp,1)))
#= """
  TriMesh_Interpolate - interpolates finite element functions at given points.

  This routine exploits the fact that the incoming mesh is Delauney.  If it 
  isn't, it is possible to land in an infinite loop.  See the functions
  TriMesh_TestDelauney and TriMesh_ElementAdjacency.

  Author: Jeff Borggaard, Virginia Tech
  Version: 1.0

  Usage:
  ```julia
    fValues = TriMesh_Interpolate(x,eConn,eAdjacency,fNodes,xInterp)
  ```

  Arguments:
  - `x`: nodal coordinates of the triangular mesh
  - `eConn`: element connectivity
  - `eAdjacency`: element adjacency
  - `fNodes`: an array of nodal values for m different functions
     (size(fNodes) = nNodes x nFunctions)
  - `xInterp`: a set of interpolation points (must be in the interior)
  - `startEl`: initial guess (optional, default: random)
  - `verbose`: verbosity level (optional, default: 0)

  Output argument:
  - `fValues`: an array of interpolated values
""" =#

  #  Get problem dimensions (we don't check the inputs)
  nNodes = size(x,1)
  nElements = size(eConn,1)
  nPoints = size(xInterp,1)
  nFunctions = size(fNodes,2)

  elementList = Array{Int64,1}(undef,nPoints)
  fValues = zeros(Float64,nPoints,nFunctions)

  for i=1:nPoints
    (verbose > 0) && println("TriMesh_Interpolate: Doing point $(i)...");
    elementList[i], baryCoord1, baryCoord2 = 
            TriMesh_Search(x,eConn[:,1:3],eAdjacency,xInterp[i,:];verbose=verbose,startEl=startEl[i]);
    
    nodesLocal = eConn[ elementList[i], : ]
    xLocal = x[nodesLocal,:]
    xg,wg,ϕ,ϕ_x,ϕ_y = twodShape(xLocal,[baryCoord1],[baryCoord2],[1.0])


    for j=1:nFunctions
      tmp = ϕ*fNodes[nodesLocal,j]
      fValues[i,j] = tmp[1]
    end
  end

  return fValues, elementList
end
