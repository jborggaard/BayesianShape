using HDF5, Statistics, Printf

# #idx = 94:109;
# #scen = "isoShape";
# idx = 1:16;
# scen = "triangleShape";
# 
# dir = "/projects/SIAllocation/stokes";
# 
# #list of files
# files = [ @sprintf("%s/%s/%s_%03d.h5",dir,scen,scen,i) for i=idx ]; 

#extraVars = [ "regularity", "nEigVals", "rMax"];#, "final_mcmc/beta" ];
extraVars = [ "regularity", "kappa", "omega"];#, "final_mcmc/beta" ];

scen = "svglobal";
dir = "/projects/SIAllocation/stokes/$(scen)";

files = [ m.match for m in match.(r".*\.h5",readdir(dir;join=true)) if m!=nothing ];
files = files[end-5:end];

#header
@printf("%22s  %8s  %5s", "file", "# samp", "a/r"); 
for var in extraVars
  @printf("  %8s",var[1:min(8,length(var))]);
end
@printf("\n");

#go file by file
for outFile in files
  sc=h5read(outFile,"sampComplete"); 
  ar=mean(h5read(outFile,"ar")); 
  #nev=h5read(outFile,"nEigVals"); 
  @printf("%22s  %8d  %5.3f", basename(outFile), sc, ar); 
  #@printf("%15s: %8d %5d %5.3f %6.4f %8.6f\n", basename(outFile), sc, nev, r, ar, beta); 

  for var in extraVars
    d = h5read(outFile,var); 
    @printf("  %8.6f",d);
  end

  @printf("\n");
end
