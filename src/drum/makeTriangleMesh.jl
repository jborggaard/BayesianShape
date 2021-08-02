function makeTriangleMesh(a; lc=5e-3, outDir=".")
#
#  Uses Gmsh to generate a mesh for a triangular domain described
#  by the points (0,0), (1,0), and (a[1],a[2])
#  

  gmsh.initialize()
  gmsh.option.setNumber("General.Terminal", 0)

  gmsh.model.add("square")

  #add points for the boundary
  gmsh.model.geo.addPoint( 0.0, 0.0, 0.0, lc,  1)
  gmsh.model.geo.addPoint( 1.0, 0.0, 0.0, lc,  2)
  gmsh.model.geo.addPoint(a[1],a[2], 0.0, lc,  3)

  #  Assemble boundary points into curves following the right-hand-rule
  #  (traverse the boundary with your right hand toward the "outside")
  gmsh.model.geo.addLine(  1,  2,  1)
  gmsh.model.geo.addLine(  2,  3,  2)
  gmsh.model.geo.addLine(  3,  1,  3)

  gmsh.model.geo.addCurveLoop([ 1, 2, 3],1)

  gmsh.model.geo.addPlaneSurface([1],1)

  # 1-dimensional; comprised of lines 1, 2, and 3; labeled group 1
  gmsh.model.addPhysicalGroup(1, [ 1, 2, 3], 1)
  # 1-dimensional; comprised of group 1; labeled "Boundary"
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
  #eConn = eConn[[1,3,2,6,5,4],:]
  eC = copy(transpose(eConn))
  eC = convert(Array{Int64,2},eC)   # for debugging

  #
  #  Extract the information for the boundary
  #  Nodes
  n1,xInner = gmsh.model.mesh.getNodesForPhysicalGroup(1,1) # 1D groups, 
  nBoundary = unique(n1)
  nBoundary = convert(Array{Int64,1},nBoundary)

  #  Faces       # 1 is for 1D groups, and 1 is the boundary 
  faceTypes, faceTags, faceConnectivity = gmsh.model.mesh.getElements(1,1)
  fConn1 = reshape(faceConnectivity[1],3,:)
  fConn1 = fConn1[[1,3,2],:]
  fConn1 = copy(transpose(fConn1))
  fC1    = convert(Array{Int64,2},fConn1)

  gmsh.write(outDir*"triangle.msh")
  gmsh.finalize()

  return x,eC, nBoundary, fC1
end
