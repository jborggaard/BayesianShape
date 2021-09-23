using Plots, SpecialFunctions

include("squashArctan.jl");
include("squashSigmoid.jl");
include("squashErf.jl");
include("squashSmoothstep.jl");
include("squashPolyinterp.jl");

r = 0:0.01:12;
rMin = 0.5; rMax = 8.0;

rat  = squashArctan(r,rMin,rMax);
rsig = squashSigmoid(r,rMin,rMax);
rerf = squashErf(r,rMin,rMax);
rss  = squashSmoothstep(r,rMin,rMax);
rpi  = squashPolyinterp(r,rMin,rMax;e=0.5);

plot(r,rat,lab="arctan");
plot!(r,rsig,lab="sigmoid");
plot!(r,rerf,lab="erf");
plot!(r,rss,lab="smoothstep");
plot!(r,rpi,lab="poly interp");
plot!(r,r,ls=:dash,lab="y=x");
plot!(r,repeat([rMin],length(r)),ls=:dash,lab="rMin");
plot!(r,repeat([rMax],length(r)),ls=:dash,lab="rMax");
plot!(r,repeat([0.5*(rMin+rMax)],length(r)),ls=:dash,lab="(rMax+rMin)/2");
plot!(leg=:bottomright);

savefig("squashCompare.png");

