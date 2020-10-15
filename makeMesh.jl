function makeMesh(r)

  n = size(r,1)

  gmsh.initialize()
  gmsh.option.setNumber("General.Terminal", 1)

  gmsh.model.add("mixing")

  lc = 3e-2
  global theta = range(0.0,stop=2.0*pi,length=n+1)

  for i=1:n
    gmsh.model.geo.addPoint(r[i]*cos(theta[i]),
                            r[i]*sin(theta[i]),
                            0.0, lc, i)
  end

  for i=1:n
    gmsh.model.geo.addPoint(2*cos(theta[i]),2*sin(theta[i]),0.0, 2*lc,n+i)
  end

  nCircle = 15
  global theta1 = range(0.0,stop=2.0*pi,length=nCircle+1)
  for i=1:nCircle
    gmsh.model.geo.addPoint(1.25+0.1*cos(theta1[i]),1.25+0.1*sin(theta1[i]),0.0, lc,2n+i)
#    gmsh.model.geo.addPoint(1.75+0.1*cos(theta1[i]),0.0+0.1*sin(theta1[i]),0.0, lc,2n+i)
  end

  #  Assemble boundary points into curves following the right-hand-rule
  #  (traverse the boundary with your right hand toward the "outside")
  gmsh.model.geo.addBSpline([1;collect(n:-1:1)],1)        
  gmsh.model.geo.addSpline([collect((n+1):(2*n));n+1],2)
  gmsh.model.geo.addSpline([2*n+1;collect(2*n+nCircle:-1:2*n+1)],3)

  gmsh.model.geo.addCurveLoop([1],1)
  gmsh.model.geo.addCurveLoop([2],2)
  gmsh.model.geo.addCurveLoop([3],3)
  gmsh.model.geo.addPlaneSurface([1,2,3],1)
  gmsh.model.geo.addPlaneSurface([3],2)

  gmsh.model.addPhysicalGroup(1, [1], 1)
  gmsh.model.addPhysicalGroup(1, [2], 2)
 # gmsh.model.addPhysicalGroup(1, [3], 3)
  gmsh.model.setPhysicalName(1, 1, "Inner")
  gmsh.model.setPhysicalName(1, 2, "Outer")

  gmsh.model.addPhysicalGroup(2, [1], 1)
  gmsh.model.addPhysicalGroup(2, [2], 2)
  gmsh.model.setPhysicalName(2, 1, "Mesh")
  # gmsh.model.setPhysicalName(2, 2, "Sensor")

  gmsh.model.geo.synchronize()

  gmsh.model.mesh.generate(2)
  gmsh.model.mesh.setOrder(2)  # quadratic elements

  ##nodeIds,nodeCoords,_ = gmsh.model.mesh.getNodes()
  ##x = reshape(nodeCoords,3,:)
  nodeTags1,xx1 = gmsh.model.mesh.getNodesForPhysicalGroup(2,1)
  nodeTags2,xx2 = gmsh.model.mesh.getNodesForPhysicalGroup(2,2)
  xx1 = reshape(xx1,3,:)
  xx2 = reshape(xx2,3,:)
  nn = convert(Int64,maximum(nodeTags2))
  x = zeros(3,nn)
  for i=1:length(nodeTags1)
    x[:,nodeTags1[i]] = xx1[:,i]
  end
  for i=1:length(nodeTags2)
    x[:,nodeTags2[i]] = xx2[:,i]
  end

  nInner,xInner = gmsh.model.mesh.getNodesForPhysicalGroup(1,1) # 1D groups, #1 is the inner circle
  nOuter,xOuter = gmsh.model.mesh.getNodesForPhysicalGroup(1,2) # 1D groups, #2 is the outer circle

  elemTypes1, elemTags1, elemConnectivity1 = gmsh.model.mesh.getElements(2,1)
  elemTypes2, elemTags2, elemConnectivity2 = gmsh.model.mesh.getElements(2,2)

  eConn1 = reshape(elemConnectivity1[1],6,:)
  eConn2 = reshape(elemConnectivity2[1],6,:)

  eConn = [eConn1 eConn2]
#  eC = convert(Array{Int64,2},eConn)   # for debugging

  gmsh.write("mixing.msh")
  gmsh.finalize()

 return x,eConn,eConn2, nInner,xInner, nOuter,xOuter
end
