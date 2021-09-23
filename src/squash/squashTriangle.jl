using Plots, SpecialFunctions

include("squashArctan.jl");
include("squashSigmoid.jl");
include("squashErf.jl");
include("squashSmoothstep.jl");
include("squashPolyinterp.jl");

include("../triangle/triangleRecenter.jl");
include("../triangle/polarTriangle.jl");


rMin = 0.2; rMax = 8.0;

#triangle vertices
vert = 12.0.*[ 0 0; 1 0; 0.635 0.275 ];
#vert = 3.0 .* [ -1.0 0.0; 1.0 0.0; 0.0 sqrt(3) ];

#recenter
vert0 = triangleRecenter(vert);

#polar coordinates of sides
x = 0:0.001:1.0;
th, r = polarTriangle(vert0,x);
r = r[:]; th = th[:];

#squash
rat  = squashArctan(r,rMin,rMax);
rsig = squashSigmoid(r,rMin,rMax);
rerf = squashErf(r,rMin,rMax);
rss  = squashSmoothstep(r,rMin,rMax);
rpi  = squashPolyinterp(r,rMin,rMax;e=0.2);

#plot
p = plot(proj=:polar,size=(600,600),leg=true);
plot!(th,r,ls=:dash,lab="true");
plot!(th,rat,lab="arctan");
plot!(th,rsig,lab="sigmoid");
plot!(th,rerf,lab="erf");
plot!(th,rss,lab="smoothstep");
plot!(th,rpi,lab="poly interp");
plot!(th,repeat([rMin],length(r)),ls=:dash,lab="rMin");
plot!(th,repeat([rMax],length(r)),ls=:dash,lab="rMax");
plot!(th,repeat([0.5*(rMin+rMax)],length(r)),ls=:dash,lab="(rMax+rMin)/2");
#plot!(leg=:topright);

savefig("squashTriangle.png");
savefig("squashTriangle.pdf");

