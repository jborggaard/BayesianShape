#computeRadii() Computes radii from Fourier samples given a Fourier basis
function computeRadii(samples,fb,a0=1.0)
    #compute unsquashed
    r = a0 .+ samples*fb;

    #squash
    α = 1.0;
    r = 1.0 .+ atan.(α*r)/π;

    return r;
end
