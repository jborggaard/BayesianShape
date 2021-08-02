function makeIsospectralMesh1(; lc=1e-2, outDir=".")
#
#  Uses Gmsh to generate a mesh on a rectangular domain
#  

  gmsh.initialize()
  gmsh.option.setNumber("General.Terminal", 0)

  gmsh.model.add("square")

  #add points for the boundary
  gmsh.model.geo.addPoint( 3.0, 0.0, 0.0, lc,  1)
  gmsh.model.geo.addPoint( 4.0, 0.0, 0.0, lc,  2)
  gmsh.model.geo.addPoint( 4.0, 1.0, 0.0, lc,  3)
  gmsh.model.geo.addPoint( 5.0, 1.0, 0.0, lc,  4)
  gmsh.model.geo.addPoint( 5.0, 2.0, 0.0, lc,  5)
  gmsh.model.geo.addPoint( 4.0, 2.0, 0.0, lc,  6)
  gmsh.model.geo.addPoint( 4.0, 4.0, 0.0, lc,  7)
  gmsh.model.geo.addPoint( 5.0, 4.0, 0.0, lc,  8)
  gmsh.model.geo.addPoint( 5.0, 5.0, 0.0, lc,  9)
  gmsh.model.geo.addPoint( 4.0, 5.0, 0.0, lc, 10)
  gmsh.model.geo.addPoint( 4.0, 7.0, 0.0, lc, 11)
  gmsh.model.geo.addPoint( 6.0, 7.0, 0.0, lc, 12)
  gmsh.model.geo.addPoint( 6.0, 6.0, 0.0, lc, 13)
  gmsh.model.geo.addPoint( 7.0, 6.0, 0.0, lc, 14)
  gmsh.model.geo.addPoint( 7.0, 7.0, 0.0, lc, 15)
  gmsh.model.geo.addPoint( 6.0, 8.0, 0.0, lc, 16)
  gmsh.model.geo.addPoint( 4.0, 8.0, 0.0, lc, 17)
  gmsh.model.geo.addPoint( 3.0, 7.0, 0.0, lc, 18)
  gmsh.model.geo.addPoint( 3.0, 5.0, 0.0, lc, 19)
  gmsh.model.geo.addPoint( 2.0, 5.0, 0.0, lc, 20)
  gmsh.model.geo.addPoint( 2.0, 4.0, 0.0, lc, 21)
  gmsh.model.geo.addPoint( 3.0, 4.0, 0.0, lc, 22)
  gmsh.model.geo.addPoint( 3.0, 2.0, 0.0, lc, 23)
  gmsh.model.geo.addPoint( 1.0, 2.0, 0.0, lc, 24)
  gmsh.model.geo.addPoint( 1.0, 3.0, 0.0, lc, 25)
  gmsh.model.geo.addPoint( 0.0, 3.0, 0.0, lc, 26)
  gmsh.model.geo.addPoint( 0.0, 2.0, 0.0, lc, 27)
  gmsh.model.geo.addPoint( 1.0, 1.0, 0.0, lc, 28)
  gmsh.model.geo.addPoint( 3.0, 1.0, 0.0, lc, 29)

  #  Assemble boundary points into curves following the right-hand-rule
  #  (traverse the boundary with your right hand toward the "outside")
  gmsh.model.geo.addLine(  1,  2,  1)
  gmsh.model.geo.addLine(  2,  3,  2)
  gmsh.model.geo.addLine(  3,  4,  3)
  gmsh.model.geo.addLine(  4,  5,  4)
  gmsh.model.geo.addLine(  5,  6,  5)
  gmsh.model.geo.addLine(  6,  7,  6)
  gmsh.model.geo.addLine(  7,  8,  7)
  gmsh.model.geo.addLine(  8,  9,  8)
  gmsh.model.geo.addLine(  9, 10,  9)
  gmsh.model.geo.addLine( 10, 11, 10)
  gmsh.model.geo.addLine( 11, 12, 11)
  gmsh.model.geo.addLine( 12, 13, 12)
  gmsh.model.geo.addLine( 13, 14, 13)
  gmsh.model.geo.addLine( 14, 15, 14)
  gmsh.model.geo.addLine( 15, 16, 15)
  gmsh.model.geo.addLine( 16, 17, 16)
  gmsh.model.geo.addLine( 17, 18, 17)
  gmsh.model.geo.addLine( 18, 19, 18)
  gmsh.model.geo.addLine( 19, 20, 19)
  gmsh.model.geo.addLine( 20, 21, 20)
  gmsh.model.geo.addLine( 21, 22, 21)
  gmsh.model.geo.addLine( 22, 23, 22)
  gmsh.model.geo.addLine( 23, 24, 23)
  gmsh.model.geo.addLine( 24, 25, 24)
  gmsh.model.geo.addLine( 25, 26, 25)
  gmsh.model.geo.addLine( 26, 27, 26)
  gmsh.model.geo.addLine( 27, 28, 27)
  gmsh.model.geo.addLine( 28, 29, 28)
  gmsh.model.geo.addLine( 29,  1, 29)

  gmsh.model.geo.addCurveLoop([ 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
                               11,12,13,14,15,16,17,18,19,20,
                               21,22,23,24,25,26,27,28,29],1)

  gmsh.model.geo.addPlaneSurface([1],1)

  gmsh.model.addPhysicalGroup(1, [ 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
                                  11,12,13,14,15,16,17,18,19,20,
                                  21,22,23,24,25,26,27,28,29], 1)
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
  #  Extract the information for the bottom boundary
  #  Nodes
  n1,xInner = gmsh.model.mesh.getNodesForPhysicalGroup(1,1) # 1D groups, 
                                                            #1 is bottom 
  nBoundary = unique(n1)
  nBoundary = convert(Array{Int64,1},nBoundary)

  #  Faces       # 1 is for 1D groups, and 1 is the bottom boundary 
  faceTypes, faceTags, faceConnectivity = gmsh.model.mesh.getElements(1,1)
  fConn1 = reshape(faceConnectivity[1],3,:)
  fConn1 = fConn1[[1,3,2],:]
  fConn1 = copy(transpose(fConn1))
  fC1    = convert(Array{Int64,2},fConn1)

  gmsh.write(outDir*"isomesh1.msh")
  gmsh.finalize()

  return x,eC, nBoundary, fC1
end
