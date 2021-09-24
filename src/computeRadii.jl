#computeRadii() Computes radii from Fourier coefficients

#this version assumes we are given a Fourier basis
function computeRadii(samples::Array{Float64,2},fb::Array{Float64,2})
    
    #compute unsquashed
    if isodd(size(samples,2))
      r = samples[:,1] + samples[:,2:end]*fb;
    else
      r = a0 .+ samples*fb;
    end

    #squash
    r = radiusSquash(r);

    return r;
end
function computeRadii(ab::Array{Float64,1},fb::Array{Float64,2})
    
    #compute unsquashed
    if isodd(length(ab))
      r = ab[1] .+ fb'*ab;
    else
      r = a0 .+ fb'*ab;
    end

    #squash
    r = radiusSquash(r);

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
function computeRadii(samples::Array{Float64,2},th::Union{Array{Float64,1},AbstractRange})
  n = size(samples,2);
  nf = isodd(n) ? (n-1)÷2 : n÷2;
  fb = fourierBasis(nf,th);
  return computeRadii(samples,fb;rMin=rMin,rMax=rMax);
end
function computeRadii(ab::Array{Float64,1},th::Union{Array{Float64,1},AbstractRange})
  n = length(ab);
  nf = isodd(n) ? (n-1)÷2 : n÷2;
  fb = fourierBasis(nf,th);
  return computeRadii(ab,fb;rMin=rMin,rMax=rMax);
end
