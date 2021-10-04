#generate plots for the scenario
#need to define outFile before running
#may need to have some other files loaded (e.g., twodStokesAD for plotting the MAP points)

using Plots
using Plots.Measures

#fix errors for headless plotting
#GKS: can't connect to GKS socket application
ENV["GKSwstype"] = "100"

include("../../src/plotSave.jl");
include("../../src/getMap.jl");
#include("../../src/radiusSquash.jl");
include("../../src/fourierBasis.jl");
include("../../src/computeRadii.jl");
include("../../src/plotRadiiQuantiles.jl");
include("../../src/plotSamplesLpdfs.jl");
include("../../src/plotQuantiles.jl");
include("../../src/drum/plotMap.jl");
include("../../src/drum/plotSampleShapes.jl");

plotQuantiles(outFile, margin=10mm);
plotSamplesLpdfs(outFile, margin=10mm);

# #assumes sampling of whole space (samples = parameters)
# plotRadiiQuantiles(outFile, margin=10mm);
# plotSampleShapes(outFile);
# plotMap(outFile;lpdfIdx=3);
# plotMap(outFile;lpdfIdx=2);

#samp to parameter map is not the identity
s = mcmcSample();
gsp = mcmcGradSampToParamMap(s); #assume constant
rMin  = h5read(outFile,"rMin");
rMax  = h5read(outFile,"rMax");
lpdfs = h5read(outFile,"lpdfs");
samples = h5read(outFile,"samples");
params = samples * gsp';

plotFile = replace(outFile,".h5"=>"_radii_quantiles");
plotRadiiQuantiles(params,plotFile; margin=10mm);
plotFile = replace(outFile,".h5"=>"_sample_shapes");
plotSampleShapes(params,plotFile);

plotFile = replace(outFile,".h5"=>"_mle");
param_mle = getMap(params,lpdfs,3);
plotSampleShapes(param_mle,plotFile);
plotFile = replace(outFile,".h5"=>"_map");
param_map = getMap(params,lpdfs,2);
plotSampleShapes(param_map,plotFile);


#include("../../src/plotMapIBs.jl");
#include("../../src/plotRadiiCorr.jl");
#include("../../src/plotRadiiHist.jl");
#include("../../src/histmatrix.jl");
#include("../../src/plotRadiiHistMatrix.jl");
#plotMapIBs(outFile);
#plotRadiiCorr(outFile, left_margin=10mm, bottom_margin=10mm);
#plotRadiiHist(outFile);
#plotRadiiHistMatrix(outFile, left_margin=10mm, size=(1600,1600));

