using Plots
using FastGaussQuadrature

include("triangleRecenter.jl");
include("polarTriangle.jl");

function projectTriangle( vert, nodes, weights; nUnk=160)
  #build triangle (quadrature points used on each face)
  th, r = polarTriangle(vert,nodes);
  
  #rescale weights
  w = weights * ( th[end,:] - th[1,:] )' ./ ( nodes[end] - nodes[1] );
  
  #collect into vectors
  thall = th[:]; rall = r[:]; wall = w[:];
  
  #basis (rescaled to be norm 1 on the circle)
  fb = zeros( length(thall), 2*nUnk+1 );
  fb[:,1] .= 1.0/sqrt(2*pi);
  fb[:,2:2:end] = cos.(thall*(1:nUnk)')/sqrt(pi);
  fb[:,3:2:end] = sin.(thall*(1:nUnk)')/sqrt(pi);
  
  #project
  proj = fb' * (wall .* rall);
  
  ##reconstruct
  #rproj = fb * proj;

  return proj;
end

#triangle vertices
vert = 5.0.*[ 0 0; 1 0; 0.635 0.275 ];
#vert = 1.0 .* [ -1.0 0.0; 1.0 0.0; 0.0 sqrt(3) ];

#quadrature points
nodes, weights = gausslegendre( 1000 );

#recenter
vert0 = triangleRecenter(vert);

#get projection
tproj = projectTriangle( vert0, nodes, weights);
println("First five components are:");
display(tproj[1:5]);

#rotate to zero out the second (first non-constant) component
psi = atan(tproj[2]/tproj[3]); 
R = [ cos(psi) sin(psi); -sin(psi) cos(psi) ]; 
vert0 = vert0 * R;
println("Rotated.");

#get projection
tproj2 = projectTriangle( vert0, nodes, weights);
println("First five components are:");
display(tproj2[1:5]);



using Plots
using Plots.Measures

rMin = 0.2;
rMax = 2*tproj[1]/sqrt(2) - rMin; #try to match first component with mean(rMin,rMax)

#fix errors for headless plotting
#GKS: can't connect to GKS socket application
ENV["GKSwstype"] = "100"

include("../../src/plotSave.jl");
include("../../src/getMap.jl");
include("../../src/radiusSquash.jl");
include("../../src/fourierBasis.jl");
include("../../src/computeRadii.jl");
include("../../src/plotRadiiQuantiles.jl");
include("../../src/plotSamplesLpdfs.jl");
include("../../src/plotQuantiles.jl");
include("../../src/drum/plotMap.jl");
include("../../src/drum/plotSampleShapes.jl");

plotSampleShapes(tproj[2:end],"src/triangle/triangleProjectRotate";rMin=rMin,rMax=rMax);
plotSampleShapes(tproj2[2:end],"src/triangle/triangleProjectRotate2";rMin=rMin,rMax=rMax);
