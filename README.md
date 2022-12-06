# BayesianShape
Software used to solve shape estimation of Stokes flows in the Bayesian framework.  An overview of the
software implemented here can be found in the reference

- *A Statistical Framework for Domain Shape Estimation in Stokes Flows*

by Jeff Borggaard, Nathan Glatt-Holtz, and Justin Krometis.

## Installation Notes
```
  git clone https://www.github.com/jborggaard/BayesianShape.git
```
will download this software.  It depends on the InfDimMCMC package that can be downloaded in julia using the commands
```
  julia> using Pkg
  julia> pkg"add https://www.github.com/krometis/InfDimMCMC.jl.git"
```

A sample job submission script using Slurm is in the file submit.sh

The three sets of experiments in the paper are found in the scenarios directory. 

- the first experiment measuring vorticity at sensors near the outer radius is setup in the *vortsensor* folder.

- the second experiment designing a mixer to minimize the scalar variance is setup in the *svglobal* folder.

- the third experiment to specify different values of the scalar variance in sectors is setup in the *svsector* folder.

## Solver Details

- *twodStokesRotatingOuter.jl* :
The main solver used to generate the Stokes flow for a given geometry.

- *twodAdvectionDiffusion.jl* :
Solves the advection diffusion equation for a given advection velocity (Stokes flow)

- *run_twodStokesRotatingOuter.jl* : 
A function that runs one sampled inner boundary.  

  - Generates a sample of Fourier coefficients describing the inner boundary *sampleInnerGeometry*, 
  
  - fits a B-Spline to the inner radius defined by the above coefficients *fitBSpline2Fourier*,
  
  - generates a mesh using Gmsh *makeMesh*,
  
  - internally builds output matrices to define the sensor,
  
  - computes the velocity field *twodStokesRotatingOuter*,
  
  - computes the scalar concentration *twodAdvectionDiffusion*,
  
  - computes the vorticity field *computeVorticity*,
  
  - internally plots solutions and computes averaged scalar quantities: passive scalar and vorticity.
  
