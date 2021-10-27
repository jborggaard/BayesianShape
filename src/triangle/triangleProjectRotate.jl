using Plots
using FastGaussQuadrature

include("triangleRecenter.jl");
include("trianglePolar.jl");
include("triangleProject.jl");

#triangle vertices
vert = 5.0.*[ 0 0; 1 0; 0.635 0.275 ];
#vert = 1.0 .* [ -1.0 0.0; 1.0 0.0; 0.0 sqrt(3) ];

#quadrature points
nodes, weights = gausslegendre( 1000 );

#recenter
vert0 = triangleRecenter(vert);

#get projection
tproj = triangleProject( vert0, nodes, weights);
println("First five components are:");
display(tproj[1:5]);

#rotate to zero out the second (first non-constant) component
psi = atan(tproj[2]/tproj[3]); 
R = [ cos(psi) sin(psi); -sin(psi) cos(psi) ]; 
vert0 = vert0 * R;
println("Rotated.");

#get projection
tproj2 = triangleProject( vert0, nodes, weights);
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
