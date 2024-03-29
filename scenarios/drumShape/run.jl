using Pkg;
pkg"activate .";
#PKG_ROOT=ENV["PKG_ROOT"];
#pkg"activate $PKG_ROOT";
using HDF5

#Parse command line arguments
using ArgParse
apSettings = ArgParseSettings();
@add_arg_table apSettings  begin
  "--scen"
    help = "scenario to run"
    required = false
    default = "drumShape"
  "--restartfile"
    help = "file to restart from"
    required = false
  "--mcmc"
    help = "mcmc method and parameters (comma-delimited)"
    required = false
  "--nsamp"
    help = "number of samples to run"
    arg_type = Int
    required = false
  "--nburn"
    help = "number of burnin samples to run"
    arg_type = Int
    required = false
  "--kappa"
    help = "XXXXXXX coefficient"
    arg_type = Float64
    required = false
  "--nev"
    help = "number of eigenvalues to observe"
    arg_type = Int
    required = false
  "--regularity"
    help = "regularity implied by the prior"
    arg_type = Float64
    required = false
  "--rmin"
    help = "minimum radius"
    arg_type = Float64
    required = false
  "--rmax"
    help = "maximum radius"
    arg_type = Float64
    required = false
  "--lc"
    help = "mesh resolution"
    arg_type = Float64
    required = false
end
args = parse_args(apSettings,as_symbols=true);

#Print command line arguments and convert them to global variables
#(This allows the setup scripts to be run without having to pass in a dictionary)
println("Parsed args:")
for (arg,val) in args
  if val != nothing
    println("  $arg  =>  $val"); #print
    @eval (($arg) = ($val));     #global variable
  end
end
if (@isdefined restartfile)
  datafileTmp = h5read(restartfile,"datafile");
  if (@isdefined datafile) && (datafileTmp != datafile)
    error("datafile is specified ($(datafile)) but does not match datafile from restartfile ($(datafileTmp))!");
  end
  datafile = datafileTmp;

  regularityTmp = h5read(restartfile,"regularity");
  if (@isdefined regularity) && (regularityTmp != regularity)
    error("regularity is specified ($(regularity)) but does not match regularity from restartfile ($(regularityTmp))!");
  end
  regularity = regularityTmp;

  nevTmp = h5read(restartfile,"nev");
  if (@isdefined nev) && (nevTmp != nev)
    error("nev is specified ($(nev)) but does not match nev from restartfile ($(nevTmp))!");
  end
  nev = nevTmp;

  kappaTmp = h5read(restartfile,"kappa");
  if (@isdefined kappa) && (kappaTmp != kappa)
    error("kappa is specified ($(kappa)) but does not match kappa from restartfile ($(kappaTmp))!");
  end
  kappa = kappaTmp;

  rminTmp = h5read(restartfile,"rMin");
  if (@isdefined rmin) && (rminTmp != rmin)
    error("rMin is specified ($(rmin)) but does not match rMin from restartfile ($(rminTmp))!");
  end
  rmin = rminTmp;

  rmaxTmp = h5read(restartfile,"rMax");
  if (@isdefined rmax) && (rmaxTmp != rmax)
    error("rMax is specified ($(rmax)) but does not match rMax from restartfile ($(rmaxTmp))!");
  end
  rmax = rmaxTmp;

  lcTmp = h5read(restartfile,"lc");
  if (@isdefined lc) && (lcTmp != lc)
    error("lc is specified ($(lc)) but does not match lc from restartfile ($(lcTmp))!");
  end
  lc = lcTmp;

  if (!@isdefined mcmc)
    mcmc = h5read(restartfile,"mcmc");
  end
end


outDir="/projects/SIAllocation/stokes/$(scen)";

if ( ! isdir(outDir) ) 
  println("Output directory $(outDir) does not exist. Creating...");
  mkdir(outDir);
end

println("Starting setup...");
include("setup.jl");
println("Finished setup.");


## Setup output file ##
mkfilename(cnt;outDir=outDir,scen=scen) = @sprintf("%s/%s_%03d.h5",outDir,scen,cnt);
cnt=1;
outFile=mkfilename(cnt);
while (isfile(outFile))
  global cnt, outFile;
  outFile=mkfilename(cnt);
  cnt+=1;
end
println("Writing output to $(outFile)...");
h5write(outFile,"datafile",datafile);
h5write(outFile,"regularity",regularity);
h5write(outFile,"kappa",kappa);
h5write(outFile,"rMin",rMin);
h5write(outFile,"rMax",rMax);
h5write(outFile,"lc",lc);
h5write(outFile,"nEigVals",nEigVals);
h5write(outFile,"obsMean",obsMean);
h5write(outFile,"obsStd",obsStd);
(@isdefined nburn) && h5write(outFile,"nburn",nburn);
(@isdefined nsamp) && h5write(outFile,"nsamp",nsamp);
(@isdefined mcmc ) && h5write(outFile,"mcmc",mcmc);


#Initial sample
function restartSample(filename)
  #return h5read(filename,"samples")[end,:];
  f = h5open(filename);
  dset = f["samples"];
  s = dset[end,:][:];
  close(f);
  return s;
end
if (@isdefined restartfile)
  s0 = restartSample(restartfile);
  h5write(outFile,"restartfile",restartfile);
  h5write(outFile,"init","restart");
else
  #s0 = sTrue[sampInd]; h5write(outFile,"init","true");
  s0 = rand(mcmcP.prior); h5write(outFile,"init","rand");
  #s0 = zeros(nSampInd); h5write(outFile,"init","zeros"); @printf("init = zeros");
end



## Run ##
mcmcTime = @elapsed mcmcRun(mcmcP, s0; verbose=3, outFile=outFile);



## Post process ##

ar      = h5read(outFile,"ar");
lpdfs   = h5read(outFile,"lpdfs");

acceptRatio = sum(ar[1:mcmcP.nsamp])/mcmcP.nsamp;
sampPerSec = mcmcP.nsamp / mcmcTime;
secPerSamp = mcmcTime / mcmcP.nsamp;

h5write(outFile,"acceptRatio",acceptRatio);
h5write(outFile,"mcmcTime",mcmcTime);
h5write(outFile,"sampPerSec",sampPerSec);
h5write(outFile,"secPerSamp",secPerSamp);

println("Wrote $(outFile).");

@printf("MCMC Time for %d samples: %.2f seconds (%.4f seconds/sample)\n",mcmcP.nsamp,mcmcTime,secPerSamp);
@printf("Acceptance ratio: %.4f\n\n",acceptRatio);

stp=max(1,round(Int,mcmcP.nsamp/20));

@printf("\nSummary of run:\n");
@printf("%9s  %9s: %8s ", "start", "end", "ar (avg)"); 
@printf("%16s %16s %16s\n", "logprior", "loglikelihood", "logposterior");
for i=0:stp:mcmcP.nsamp-1
  @printf("%9d -%9d:",i+1,i+stp);
  @printf(" %8.4f",mean( ar[i+1:i+stp] ));
  @printf(" %16.6f %16.6f %16.6f\n", lpdfs[i+stp,1], lpdfs[i+stp,2], lpdfs[i+stp,3]);
end


## Plots ##
include("plot.jl");
