include("../vortsensor/setup.jl");

#datafile = "/projects/SIAllocation/stokes/vortsensor/vortsensor_017.h5";
datafile = "rand";

outfile = "data.h5";

#create sample
s = mcmcSample(); 
#sample value
if datafile != "rand"
  f = h5open(datafile);
  s.samp = f["samples"][end,:];
  close(f);
else
  s.samp = rand(mcmcP.prior);
end
#run forward to get observations
#compute the parameter associated with the sample
s.param = mcmcSampToParamMap(s);
#solution
s.sol   = mcmcForwardMap(s);
#observation
s.obs   = mcmcObsMap(s);

#write output
f = h5open(outfile,"w");
write(f,"datafile",datafile);
write(f,"sample",s.samp);
write(f,"obs",s.obs);
write(f,"omega",omega);
write(f,"rMin",rMin);
write(f,"rMax",rMax);
write(f,"a0",a0);
write(f,"regularity",regularity);
write(f,"squashE",squashE);
close(f); 
println("Wrote: $(outfile)");
