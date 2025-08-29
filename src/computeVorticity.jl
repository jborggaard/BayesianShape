function computeVorticity(x,eConn,velocity)

  dVelocitydx,dVelocitydy,_ = TriMesh_ProjectDerivatives(x,eConn,velocity,[])

  vorticity = dVelocitydx[:,2] - dVelocitydy[:,1]  # dvdx - dudy

  return vorticity
end
