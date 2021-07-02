function makeDrumMesh(r;lc=7e-3,outDir=".",verbose=0)
  #  Reads in Fourier coefficients r and generates a finite element
  #  mesh for the "hearing the shape of a fuzzy drum" problem.

  #  lc (mesh density at the outer boundary)
  gmsh.initialize()
  gmsh.option.setNumber("General.Terminal", verbose)

  gmsh.model.add("drum")

  #  Build points for the outer boundary (using the given Fourier coefficients)
  n = size(r,1)
  theta = range(0.0,stop=2.0*pi,length=n+1)

  for i=1:n
    gmsh.model.geo.addPoint(r[i]*cos(theta[i]), 
                            r[i]*sin(theta[i]), 
                            0.0, 
                            lc, i)
  end

  #  Assemble the boundary points into curves following the right-hand-rule
  gmsh.model.geo.addBSpline([1;collect(n:-1:1)],1)         

  #  Assemble boundary points into curves following the right-hand-rule
  #  (traverse the boundary with your right hand toward the "outside")

  gmsh.model.geo.addCurveLoop([1],1)
 
  gmsh.model.geo.addPlaneSurface([-1],1)

  gmsh.model.addPhysicalGroup(1, [1], 1)
  gmsh.model.setPhysicalName(1, 1, "Boundary")

  gmsh.model.addPhysicalGroup(2, [1], 1)
  gmsh.model.setPhysicalName(2, 1, "Mesh")

  gmsh.model.geo.synchronize()

  gmsh.model.mesh.generate(2)
  gmsh.model.mesh.setOrder(2)  # quadratic elements

  #
  #  Extract the nodes and elements of the entire domain
  #
  ##nodeIds,nodeCoords,_ = gmsh.model.mesh.getNodes()
  ##x = reshape(nodeCoords,3,:)
  nodeTags1,xx1 = gmsh.model.mesh.getNodesForPhysicalGroup(2,1)
  xx1 = reshape(xx1,3,:)
  nn = convert(Int64,maximum(nodeTags1))
  x = zeros(3,nn)
  for i=1:length(nodeTags1)
    x[:,nodeTags1[i]] = xx1[:,i]
  end
  x = copy(transpose(x))

  elemTypes, elemTags, elemConnectivity = gmsh.model.mesh.getElements(2,1)
  eConn = reshape(elemConnectivity[1],6,:)
  eC = copy(transpose(eConn))
  eC = convert(Array{Int64,2},eC)   # for debugging
  eC = eC[:,[1,3,2,6,5,4]]

  #
  #  Extract the nodes of the boundary
  n1,xInner = gmsh.model.mesh.getNodesForPhysicalGroup(1,1) # 1D groups, #1 is on the inner circle

  nBoundary = unique(n1)
  nBoundary = convert(Array{Int64,1},nBoundary)
  sort!(nBoundary)

#  gmsh.write(outDir*"concentric.msh")
  gmsh.finalize()

  return x,eC, nBoundary
end
