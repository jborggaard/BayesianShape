#Define data type for mesh/solutions to stokes and a-d equations
mutable struct solutionArray
  x::Array{Float64,2}               #mesh
  eConn::Array{Int64,2}             #mesh connectivity
  eConn2::Array{Array{Int64,2},1}   #mesh connectivity (sensors)
  velocity::Array{Float64,2}        #velocity (stokes)
#  pressure::Array{Float64,2}        #pressure (stokes)
  temperature::Array{Float64,2}     #temperature (a-d)
  #massMat::SparseMatrixCSC{Float64,Int64}      #mass matrix

  solutionArray() = new()
end

import Base.copy
function copy(s1::solutionArray)
  s2 = solutionArray();
  s2.x     = copy(s1.x);
  s2.eConn = copy(s1.eConn);
  s2.velocity     = copy(s1.velocity);
#  s2.pressure     = copy(s1.pressure);
  s2.temperature = copy(s1.temperature);
  #s2.massMat = copy(s1.massMat);
  return s2;
end

