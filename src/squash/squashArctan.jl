#squashArctan() uses atan() to squeeze radii to a given (rMin,rMax)
function squashArctan(r, rMin, rMax)
    width = rMax-rMin;
    rMean = 0.5*(rMin+rMax);
    rc = (r.-rMean)./width; #recentered with mean at zero

    #squash
    rs = rMean .+ width*atan.(π*rc)/π;

    return rs;
end

