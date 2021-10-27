using FastGaussQuadrature

include("../triangle/triangleRecenter.jl");
include("../triangle/trianglePolar.jl");
include("../triangle/triangleProject.jl");

function triangleCoefficients(vert = 5.0.*[ 0 0; 1 0; 0.635 0.275 ])
  #quadrature points
  nodes, weights = gausslegendre( 1000 );
  
  #recenter
  vert0 = triangleRecenter(vert);
  
  #get projection
  tproj = triangleProject( vert0, nodes, weights);
  #println("First five components are:");
  #display(tproj[1:5]);
  
  #rotate to zero out the second (first non-constant) component
  psi = atan(tproj[2]/tproj[3]); 
  R = [ cos(psi) sin(psi); -sin(psi) cos(psi) ]; 
  vert0 = vert0 * R;
  #println("Rotated.");
  
  #get projection
  tproj2 = triangleProject( vert0, nodes, weights);
  #println("First five components are:");
  #display(tproj2[1:5]);

  #convert the first basis from 1/sqrt(2) => 1
  tproj2[1] /= sqrt(2);
  
  return tproj2;
end
