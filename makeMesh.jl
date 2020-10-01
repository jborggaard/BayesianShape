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

  gmsh.model.geo.addBSpline([1;collect(n:-1:1)],1)
  gmsh.model.geo.addSpline([collect((n+1):(2*n));n+1],2)

  gmsh.model.geo.addCurveLoop([1],1)
  gmsh.model.geo.addCurveLoop([2],2)
  gmsh.model.geo.addPlaneSurface([1,2],1)

  gmsh.model.addPhysicalGroup(1, [1], 1)
  gmsh.model.addPhysicalGroup(1, [2], 2)
  gmsh.model.addPhysicalGroup(2, [1], 1)
  gmsh.model.setPhysicalName(1, 1, "Inner")
  gmsh.model.setPhysicalName(1, 2, "Outer")

  gmsh.model.geo.synchronize()

  gmsh.model.mesh.generate(2)
  gmsh.model.mesh.setOrder(2)  # quadratic elements

  ##nodeIds,nodeCoords,_ = gmsh.model.mesh.getNodes()
  ##x = reshape(nodeCoords,3,:)
  nodeTags,xx = gmsh.model.mesh.getNodesForPhysicalGroup(2,1)
  xx = reshape(xx,3,:)
  x = zeros(size(xx))
  for i=1:length(nodeTags)
      x[:,nodeTags[i]] = xx[:,i]
  end

  nInner,xInner = gmsh.model.mesh.getNodesForPhysicalGroup(1,1) # 1D groups, #1 is the inner circle
  nOuter,xOuter = gmsh.model.mesh.getNodesForPhysicalGroup(1,2) # 1D groups, #2 is the outer circle

  elemTypes, elemTags, elemConnectivity = gmsh.model.mesh.getElements(2,1)

  eConn = reshape(elemConnectivity[1],6,:)
#  eC = convert(Array{Int64,2},eConn)   # for debugging

  #gmsh.write("mixing.msh")
  gmsh.finalize()

  return x,eConn, nInner,xInner, nOuter,xOuter
end
