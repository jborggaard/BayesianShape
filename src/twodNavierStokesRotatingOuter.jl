function twodNavierStokesRotatingOuter(x,eConn,innerNodes,outerNodes,ω,max_iter,min_residual,μ)
#
#  Solves the Navier-Stokes equation in 2D with Dirichlet boundary conditions
#  - ∇⋅(∇z+∇z')+ z.Δz + ∇p = f,  Ω is the domain between a Bspline disk and an
#  outer circle.  We prescribe counter-clockwise rotation to the outer circle. 
#  The solve is done using a Newton scheme with the following parameters
#         
#
#
# max_iter:                        Maximum number of iterations in the newton method
#
# min_residual:                    Residual value at convergence
#
# velocity_stokes:                 Warm start for velocity
#
# pressure_stokes:                 Warm start for pressure
#
# velocity_Nstokes:                Resulting velocity from Newton method
#
# pressure_Nstokes:                Resulting pressure from Newton method
#
# μ:                               μ=1/Re
#
#---------------------------------------------------------------------------78--


# We first compute the solution of the stokes (to be using as a warm start for the Newton method) 
velocity_Nstokes,pressure_Nstokes = twodStokesRotatingOuter(x,eConn,innerNodes,outerNodes,ω,μ);

dvelocity=velocity_Nstokes
dpressure=pressure_Nstokes

Residual = norm([dvelocity dpressure])
n=0;

while Residual > min_residual && n < max_iter


local dvelocity
local dpressure

dvelocity,dpressure = twodNavierStokesRotatingOuterNewton(x,eConn,innerNodes,outerNodes,ω,velocity_Nstokes,pressure_Nstokes,μ)

velocity_Nstokes = velocity_Nstokes + dvelocity
pressure_Nstokes = pressure_Nstokes + dpressure

n=n+1

Residual=norm([dvelocity[:];dpressure])



println("Iteration number $n and residual is $Residual")


end

return velocity_Nstokes,pressure_Nstokes

end
