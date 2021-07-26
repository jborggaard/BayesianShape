#computeRadii() Computes radii from Fourier coefficients

#this version assumes we are given a Fourier basis
function computeRadii(samples::Array{Float64,2},fb::Array{Float64,2};rMin=0.5,rMax=1.5, α=π/(rMax-rMin))
    #mean
    a0 = 0.5*(rMax+rMin);
    
    #compute unsquashed
    r = a0 .+ samples*fb;

    #squash
    r = radiusSquash(r; rMin=rMin, rMax=rMax, α=α);

    return r;
end
function computeRadii(ab::Array{Float64,1},fb::Array{Float64,2};rMin=0.5,rMax=1.5, α=π/(rMax-rMin))
    #mean
    a0 = 0.5*(rMax+rMin);
    
    #compute unsquashed
    r = a0 .+ fb'*ab;

    #squash
    r = radiusSquash(r; rMin=rMin, rMax=rMax, α=α);

    return r;
end

#this version assumes we are given array of angles
#function computeRadii(samples,th::Array{Float64,1};rMin=0.5,rMax=1.5, α=π/(rMax-rMin))
#  fb = fourierBasis(size(samples,2)÷2,th);
#  return computeRadii(samples,fb;rMin=rMin,rMax=rMax, α=α);
#end
#function computeRadii(samples,th::AbstractRange;rMin=0.5,rMax=1.5, α=π/(rMax-rMin))
#  fb = fourierBasis(size(samples,2)÷2,th);
#  return computeRadii(samples,fb;rMin=rMin,rMax=rMax, α=α);
#end
function computeRadii(samples::Array{Float64,2},th::Union{Array{Float64,1},AbstractRange};rMin=0.5,rMax=1.5, α=π/(rMax-rMin))
  fb = fourierBasis(size(samples,2)÷2,th);
  return computeRadii(samples,fb;rMin=rMin,rMax=rMax, α=α);
end
function computeRadii(ab::Array{Float64,1},th::Union{Array{Float64,1},AbstractRange};rMin=0.5,rMax=1.5, α=π/(rMax-rMin))
  fb = fourierBasis(length(ab)÷2,th);
  return computeRadii(ab,fb;rMin=rMin,rMax=rMax, α=α);
end
