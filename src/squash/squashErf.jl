#squashErf() uses erf() to squeeze radii to a given (rMin,rMax)
using SpecialFunctions
function squashErf(r, rMin, rMax)
    width = rMax-rMin;
    rMean = 0.5*(rMin+rMax);
    rc = (r.-rMean)./width; #recentered with mean at zero

    #squash
    rs = rMean .+ width*erf.(sqrt(π)*rc)/2;

    return rs;
end

