## For a given scenario, draw random samples from the prior and compute the associated observations
# Example usage:
#   include("scenarios/svsector/setup.jl");
#   samples, obs = generateSampleObs(mcmcP);

using Printf
using InfDimMCMC

function generateSampleObs(mcmcP::mcmcProb; nSamples=10, printObs=true, outFile="none")

  ##don't need to compute gradients
  #mcmc = "pcn|0.25";

  ##load the scenario
  #scen="scenarios/$(scen)/setup.jl";
  #println("Loading $(scen)");
  #include(scen);

  #don't want to compute gradients
  mcmcP.computeGradients = false;

  #initialize a sample
  s = mcmcSample();

  #draw sample and compute observations to get dimensions
  s.samp = rand(mcmcP.prior);
  mcmcFillSample(s,mcmcP);

  #initialize arrays
  samples = zeros(nSamples,length(s.samp));
  obs     = zeros(nSamples,length(s.obs));

  #waste not
  samples[1,:] = s.samp;
  obs[1,:]     = s.obs;

  #compute the rest
  for i=2:nSamples
    s.samp = rand(mcmcP.prior);
    mcmcFillSample(s,mcmcP);
    samples[i,:] = s.samp;
    obs[i,:]     = s.obs;
  end

  #print
  if printObs
    @printf("observations:\n");
    for i=1:nSamples
      @printf("%3d: ",i);
      for j=1:size(obs,2)
        @printf("%12.8f ",obs[i,j]);
      end
      @printf("\n");
    end
  end

  if outFile != "none"
    h5write(outFile,"samples",samples);
    h5write(outFile,"obs",obs);
    println("Wrote: $(outFile)");
  end

  return samples, obs;
end

