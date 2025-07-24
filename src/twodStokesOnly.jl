## Given a Fourier parametrization of the inner boundary, compute the solutions to the Stokes and Advection-Diffusion problems
#  Inputs:
#    a        coefficients of cosines
#    b        coefficients of sines
#    a0       mean inner boundary
#    N        number of B splines used to approximate Fourier representation
#    ω        angular velocity
#    verbose  verbosity
#
#  Outputs:
#    sa       structure describing the mesh and the velocity, pressure, and temperature fields

function twodStokesOnly(a,b,a0,N; ω = 10.0, verbose=true, circleCenters=[], μ = 0.001)
  
  r,err = fitBSpline2Fourier(a0,a,b,N);
  #r = ones(N,1); # uncomment for some degubbing
  verbose && @printf("B-spline approximation error (%d B-splines) is: %12.8f\n",N,err);
  
  #Squash to the min/max radius
  r = radiusSquash(r);
  
  #Generate the finite element mesh using Gmsh (implemented in makeMesh)
  x,eConn,eConn2, innerNodes,innerX, outerNodes,outerX = makeMesh(r;circleCenters=circleCenters);
  
  #compute Stokes flow
  velocity,pressure = twodStokesRotatingOuter(x,eConn,innerNodes,outerNodes,ω, μ);
  # #solve steady Advection-Diffusion equation
  # temperature = twodAdvectionDiffusion(x,eConn,innerNodes,outerNodes,velocity,κ,sourceXY)
  
  #squish everything we might need into a structure
  sa = solutionArray();
  sa.x      = x;
  sa.eConn  = eConn;
  sa.eConn2 = eConn2;
  sa.velocity = velocity;
  sa.pressure = pressure;
  #sa.temperature = temperature;

  return sa;
end
