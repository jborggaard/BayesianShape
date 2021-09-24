#radiusSquash() uses atan() to squeeze radii to a given (rMin,rMax)
#Default α matches derivative at r=rMean
#Default rMean assumes that r is centered on (rMin,rMax) and need only be squashed
#  (i.e., that r is not already mean zero)
function radiusSquash(r; rMin=0.5, rMax=1.5, α=π/(rMax-rMin), rMean=0.5*(rMin+rMax))
    # #mean
    # a0 = 0.5*(rMax+rMin);
    
    #squash
    #rs = a0 .+ (rMax-rMin)*atan.(α*(r.-rMean))/π;
    rs = rMean .+ (rMax-rMin)*atan.(α*(r.-rMean))/π;

    return rs;
end

