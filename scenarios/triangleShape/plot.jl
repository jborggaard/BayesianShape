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
include("../../src/radiusSquash.jl");
include("../../src/fourierBasis.jl");
include("../../src/computeRadii.jl");
include("../../src/plotRadiiQuantiles.jl");
include("../../src/plotSamplesLpdfs.jl");
include("../../src/plotQuantiles.jl");
include("../../src/drum/plotMap.jl");
include("../../src/drum/plotSampleShapes.jl");
plotRadiiQuantiles(outFile, margin=10mm);
plotQuantiles(outFile, margin=10mm);
plotSampleShapes(outFile);
plotSamplesLpdfs(outFile, margin=10mm);
plotMap(outFile;lpdfIdx=3);
plotMap(outFile;lpdfIdx=2);

#include("../../src/plotMapIBs.jl");
#include("../../src/plotRadiiCorr.jl");
#include("../../src/plotRadiiHist.jl");
#include("../../src/histmatrix.jl");
#include("../../src/plotRadiiHistMatrix.jl");
#plotMapIBs(outFile);
#plotRadiiCorr(outFile, left_margin=10mm, bottom_margin=10mm);
#plotRadiiHist(outFile);
#plotRadiiHistMatrix(outFile, left_margin=10mm, size=(1600,1600));

