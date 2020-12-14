#fourierBasis() computes a Fourier basis given a number of basis functions and a set of angles
function fourierBasis(nf,th)
    fb = zeros(2*nf,length(th));
    
    fb[1:2:end,:] = cos.((1:nf)*th');
    fb[2:2:end,:] = sin.((1:nf)*th');
        
    return fb;
end
