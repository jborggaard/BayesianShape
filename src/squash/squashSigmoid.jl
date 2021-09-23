#squashSigmoid() uses a sigmoid to squeeze radii to a given (rMin,rMax)
function squashSigmoid(r, rMin, rMax)
    width = rMax-rMin;
    rMean = 0.5*(rMin+rMax);
    rc = (r.-rMean)./width; #recentered with mean at zero

    #squash
    rs = rMin .+ width./(1 .+ exp.(-4*rc));

    return rs;
end

