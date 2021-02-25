function makeMesh(r; outDir=".", circleCenters=[], nCirclePts = 15, circleRadius = 0.1)

  n = size(r,1)

  gmsh.initialize()
  gmsh.option.setNumber("General.Terminal", 1)

  gmsh.model.add("mixing")

  lc = 3e-2
  theta = range(0.0,stop=2.0*pi,length=n+1)
  # #counterclockwise -> clockwise
  # theta = range(2.0*pi,stop=0.0,length=n+1)

  cnt=0;

  #add points for inner boundary
  #for i=1:n
  #  gmsh.model.geo.addPoint(r[i]*cos(theta[i]),
  #                          r[i]*sin(theta[i]),
  #                          0.0, lc, i)
  #counterclockwise -> clockwise
  for i=n:-1:1
    gmsh.model.geo.addPoint(r[i]*cos(theta[i]),
                            r[i]*sin(theta[i]),
                            0.0, lc, n-i+1)
    cnt+=1;
  end

  #add points for outer boundary
  #for i=1:n
  #  gmsh.model.geo.addPoint(2*cos(theta[i]),2*sin(theta[i]),0.0, 2*lc,n+i)
  #counterclockwise -> clockwise
  for i=n:-1:1
    gmsh.model.geo.addPoint(2*cos(theta[i]),2*sin(theta[i]),0.0, 2*lc,2*n-i+1)
    cnt+=1;
  end

  #add points for sensors (circles)
  #cAngles = range(0.0,stop=2.0*pi,length=nCirclePts+1)
  cAngles = range(2.0*pi,stop=0.0,length=nCirclePts+1) #reverse: counterclockwise -> clockwise
  circleOffset = circleRadius .* [ cos.(cAngles) sin.(cAngles) ];
  nCircles = size(circleCenters,1);
  #println("nCircles = $(nCircles)");
  for j=1:nCircles
    circlePts = circleCenters[j,:]' .+ circleOffset; #broadcast?
    for i=1:nCirclePts
      gmsh.model.geo.addPoint(circlePts[i,1],circlePts[i,2],0.0, lc,2*n+nCirclePts*(j-1)+i)
      cnt+=1;
      #gmsh.model.geo.addPoint(1.25+0.1*cos(cAngles[i]),1.25+0.1*sin(cAngles[i]),0.0, lc,2n+i)
#      gmsh.model.geo.addPoint(1.75+0.1*cos(cAngles[i]),0.0+0.1*sin(cAngles[i]),0.0, lc,2n+i)
    end
  end
  println("Added $(cnt) points");

  #  Assemble boundary points into curves following the right-hand-rule
  #  (traverse the boundary with your right hand toward the "outside")
  gmsh.model.geo.addBSpline([1;collect(n:-1:1)],1)        
  gmsh.model.geo.addSpline([collect((n+1):(2*n));n+1],2)
  #gmsh.model.geo.addSpline([2*n+1;collect(2*n+nCirclePts:-1:2*n+1)],3)
  for j=1:nCircles
    gmsh.model.geo.addSpline([2*n+nCirclePts*(j-1)+1;collect(2*n+nCirclePts*j:-1:2*n+nCirclePts*(j-1)+1)],j+2)
  end

  gmsh.model.geo.addCurveLoop([1],1)
  gmsh.model.geo.addCurveLoop([2],2)
  for j=1:nCircles
    gmsh.model.geo.addCurveLoop([j+2],j+2)
  end
  #gmsh.model.geo.addPlaneSurface([1,2,3],1)
  #gmsh.model.geo.addPlaneSurface([3],2)
  gmsh.model.geo.addPlaneSurface(collect(1:nCircles+2),1)
  for j=1:nCircles
    gmsh.model.geo.addPlaneSurface([j+2],j+1)
  end

  gmsh.model.addPhysicalGroup(1, [1], 1)
  gmsh.model.addPhysicalGroup(1, [2], 2)
 # gmsh.model.addPhysicalGroup(1, [3], 3)
  gmsh.model.setPhysicalName(1, 1, "Inner")
  gmsh.model.setPhysicalName(1, 2, "Outer")

  gmsh.model.addPhysicalGroup(2, [1], 1)
  gmsh.model.setPhysicalName(2, 1, "Mesh")
  #gmsh.model.addPhysicalGroup(2, [2], 2)
  ## gmsh.model.setPhysicalName(2, 2, "Sensor")
  for j=1:nCircles
    gmsh.model.addPhysicalGroup(2, [j+1], j+1)
    gmsh.model.setPhysicalName(2, j+1, "Sensor$(j)")
  end

  gmsh.model.geo.synchronize()

  gmsh.model.mesh.generate(2)
  gmsh.model.mesh.setOrder(2)  # quadratic elements

  ##nodeIds,nodeCoords,_ = gmsh.model.mesh.getNodes()
  ##x = reshape(nodeCoords,3,:)
  nodeTags1,xx1 = gmsh.model.mesh.getNodesForPhysicalGroup(2,1)
  nodeTags2,xx2 = gmsh.model.mesh.getNodesForPhysicalGroup(2,2)
  xx1 = reshape(xx1,3,:)
  xx2 = reshape(xx2,3,:)
  nn = convert(Int64,maximum([nodeTags1;nodeTags2]))
  x = zeros(3,nn)
  for i=1:length(nodeTags1)
    x[:,nodeTags1[i]] = xx1[:,i]
  end
  for i=1:length(nodeTags2)
    x[:,nodeTags2[i]] = xx2[:,i]
  end

  #circleNFPG = [ gmsh.model.mesh.getNodesForPhysicalGroup(2,j+1) for j=1:nCircles ];
  #nCircleNodes = sum( [ size(cnfpg[2],1) for cnfpg in circleNFPG ] );

  #figure out how long x needs to be
  #nn = 0
  #for j=1:nCircles+1
  #  nodeTags,xx = gmsh.model.mesh.getNodesForPhysicalGroup(2,j)
  #  nn = max(nn,convert(Int64,maximum(nodeTags)))
  #end
  nodeTags,_ = gmsh.model.mesh.getNodesForPhysicalGroup(2,nCircles+1)
  nn = convert(Int64,maximum(nodeTags))
  #now assemble the list
  x = zeros(3,nn)
  for j=1:nCircles+1
    nodeTags,xx = gmsh.model.mesh.getNodesForPhysicalGroup(2,j)
    xx = reshape(xx,3,:)
    for i=1:length(nodeTags)
      x[:,nodeTags[i]] = xx[:,i]
    end
  end

  nInner,xInner = gmsh.model.mesh.getNodesForPhysicalGroup(1,1) # 1D groups, #1 is the inner circle
  nOuter,xOuter = gmsh.model.mesh.getNodesForPhysicalGroup(1,2) # 1D groups, #2 is the outer circle

  #  elemTypes1, elemTags1, elemConnectivity1 = gmsh.model.mesh.getElements(2,1)
  #  elemTypes2, elemTags2, elemConnectivity2 = gmsh.model.mesh.getElements(2,2)
  #
  #  eConn1 = reshape(elemConnectivity1[1],6,:)
  #  eConn2 = reshape(elemConnectivity2[1],6,:)
  #
  #  eConn = [eConn1 eConn2]
  ##  eC = convert(Array{Int64,2},eConn)   # for debugging
  
  #make an array of arrays containing the element connectivities
  eConnArray = [ reshape(gmsh.model.mesh.getElements(2,j)[3][1],6,:) for j=1:nCircles+1 ];

  #ATTEMPT 1: (eConn2 no longer an array of arrays, which we want)
  # #eConn2 = [];
  # #for j=1:nCircles
  # #  eConn2 = [ eConn2 eConnArray[j+1]];
  # #end
  # if nCircles>0
  #   eConn2 = eConnArray[2];
  #   for j=3:nCircles
  #     eConn2 = [ eConn2 eConnArray[j+1]];
  #   end
  # end
  # eConn = [ eConnArray[1] eConn2 ];

  #ATTEMPT 2:
  #eConn = zeros( 6, sum( [ size(eC,2) for eC in eConnArray ] ) );
  #cnt=0;
  #for eC in eConnArray
  #  eConn[:,cnt+1:cnt:size(eC,2)] = eC;
  #  cnt+=size(eC,2);
  #end
  #eConn2 = eConnArray[2:end];
  
  #ATTEMPT 3:
  #concatenate them into a single array
  eConn  = hcat(eConnArray...);
  #keep the arrays for the small circles
  eConn2 = eConnArray[2:end];

  gmsh.write(outDir*"mixing.msh")
  gmsh.finalize()
  
  #Sort
  sort!(nInner);
  sort!(nOuter);

  #Transpose (short and wide -> tall and thin)
  #Also:
  #      x: Leave behind third dimension of zeros
  #  eConn: Convert to Int64
  xT = Matrix(transpose(x[1:2,:])); 
  eConnT  = Matrix(transpose(Int64.(eConn)));
  eConn2T = [ Matrix(transpose(Int64.(eC2))) for eC2 in eConn2 ];

  #return x,eConn,eConn2, nInner,xInner, nOuter,xOuter
  return xT,eConnT,eConn2T, nInner,xInner, nOuter,xOuter
end
