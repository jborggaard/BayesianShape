#  Solves the Navier-Stokes equation in 2D with Dirichlet boundary conditions
#     - ∇⋅(∇z+∇z') + z.Δz + ∇p = f,  Ω is the domain between a Bspline disk and an
#  outer circle.  We prescribe counter-clockwise rotation to the outer circle.
#---------------------------------------------------------------------------78--

using Gmsh:gmsh
using LinearAlgebra
using SparseArrays
using SpecialMatrices
using Polynomials
using Printf
using Random
using WriteVTK
using Distributions
using HDF5
using InfDimMCMC
using FEMfunctions
using Plots
using CairoMakie

include("makeMesh.jl")
include("fitBSpline2Fourier.jl")
include("saveFEMasVTK.jl")
include("twodStokesRotatingOuter.jl")
include("twodAdvectionDiffusion.jl")
include("twodProjectDerivatives.jl")
include("computeC.jl")
include("computeVorticity.jl")
include("solutionArray.jl");
include("twodStokesAD.jl");
include("twodNavierStokesRotatingOuterNewton.jl")
include("twodNavierStokesRotatingOuter.jl")
include("twodNavierStokesAD.jl")
include("twodNavierStokesOnly.jl")




def_rmin   = 0.5;
def_rmax   = 1.5;
rMin    = (@isdefined rmin   ) ? rmin   : def_rmin;
rMax    = (@isdefined rmax   ) ? rmax   : def_rmax;
squashE = 0.1;
include("../../BayesianShape/src/squash/squashPolyinterp.jl");
radiusSquash(r) = squashPolyinterp(r,rMin,rMax;e=squashE);
squashMethod="squashPolyinterp, e=$(squashE)";
println("Squashing with: $(squashMethod)");
#=
#a=[0.1,0.3,0.1]
#b=[0.3,0.1,0.3]

a=[0.0,0.1,0.3]
b=[0.0,0.0,0.0]


regularity = 0.1
unkDim = 30
p = 2*regularity + 1; #see Dashti-Stuart Thm 2.12 
sinCosStd = (1:unkDim).^(-0.5*p); #2.0.^(-(0:unkDim-1)./4); 
#sinCosStd = 2.0.^(-(0:unkDim-1)./4); 
prStd = zeros(2*unkDim);
prStd[1:2:end] = sinCosStd; #cos
prStd[2:2:end] = sinCosStd; #sin
prior = MvNormal(zeros(2*unkDim),prStd);

sample = rand(prior)

#a = sample[1:2:end]; 
#b = sample[2:2:end];
=#

a0=1/2

Re = 200

N = 10
ω = 10
circleCenters=[]

r,err = fitBSpline2Fourier(a0,a,b,N);
#Squash to the min/max radius
r = radiusSquash(r);
  
#Generate the finite element mesh using Gmsh (implemented in makeMesh)
x,eConn,eConn2, innerNodes,innerX, outerNodes,outerX = makeMesh(r;circleCenters=circleCenters);

max_iter = 100  
min_residual = 10^-10

μ = 1/Re


#sa  =  twodNavierStokesOnly(a,b,a0,N; ω = 10.0,          verbose=true, circleCenters=[],                   max_iter,min_residual,μ)
sa1 =  twodNavierStokesAD(a,b,a0,N; ω = 10.0, κ = 1.0, verbose=true, circleCenters=[],sourceXY=[1.5;1.0],max_iter,min_residual,μ)

velocity_Nstokes=sa1.velocity
pressure_Nstokes=sa1.pressure
temperature_Nstokes=sa1.temperature
vorticity = computeVorticity(x, eConn, velocity_Nstokes)


pidx = rand(1:size(x,1),1000);
quivScale=0.01;
qxy = x[pidx,:]; 
vf  = velocity_Nstokes[pidx,:];
qvf = quivScale .* vf;
Plots.quiver(qxy[:,1],qxy[:,2],quiver=(qvf[:,1],qvf[:,2]),marker=(:none),color=:black)






figsize = 800
# Create a figure
fig = Figure(size=(figsize*1.2, figsize))

# Create an axis inside the figure
ax = Axis(fig[1, 1], aspect=1, xlabel="x", ylabel="y",title="Temperature Field for Re=$Re")

# Plot your colored triangles (poly returns a Poly instance you can attach a colorbar to)
plt = poly!(ax, x, eConn[:,1:3],
            color = temperature_Nstokes[:,1],
            colormap=:turbo,
            colorrange=[minimum(temperature_Nstokes),maximum(temperature_Nstokes)],
            strokecolor = (:black, 0.6),
            strokewidth = 0.2)

# Add a colorbar — tell it which plot object to link to
Colorbar(fig[1, 2], plt, label = "Temperature")

# Display figure
fig





figsize = 800
# Create a figure
fig = Figure(size=(figsize*1.2, figsize))

# Create an axis inside the figure
ax = Axis(fig[1, 1], aspect=1, xlabel="x", ylabel="y",title="Vorticity for Re=$Re")

# Plot your colored triangles (poly returns a Poly instance you can attach a colorbar to)
plt = poly!(ax, x, eConn[:,1:3],
            color = vorticity[:,1],
            colormap=:turbo,
            colorrange=[-5,0],
            strokecolor = (:black, 0.6),
            strokewidth = 0.2)

# Add a colorbar — tell it which plot object to link to
Colorbar(fig[1, 2], plt, label = "Vorticity")

# Display figure
fig






figsize = 800
fig = Figure(size=(figsize, figsize))

# Temperature
ax1 = Axis(fig[1, 1], aspect=1, xlabel="x", ylabel="y",title="Temperature Field for Re=$Re")
plt1 = poly!(ax1, x, eConn[:,1:3],
             color = temperature_Nstokes[:,1],
             colormap=:turbo,
             colorrange=[0,3],
             strokecolor = (:black, 0.6),
             strokewidth = 0.2)
Colorbar(fig[1, 2], plt1, label = "Temperature")

# Vorticity
ax2 = Axis(fig[2, 1], aspect=1, xlabel="x", ylabel="y",title="Vorticity for Re=$Re")
plt2 = poly!(ax2, x, eConn[:,1:3],
             color = vorticity[:,1],
             colormap=:turbo,
             colorrange=[-5,110],
             strokecolor = (:black, 0.6),
             strokewidth = 0.2)
Colorbar(fig[2, 2], plt2, label = "Vorticity")

# Display combined figure
fig
