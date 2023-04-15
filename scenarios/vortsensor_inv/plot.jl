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
include("../../src/computeRadii.jl");
include("../../src/computeRadii.jl");
include("../../src/fourierBasis.jl");
include("../../src/plotMapIBs.jl");
include("../../src/plotRadiiQuantiles.jl");
include("../../src/plotRadiiCorr.jl");
include("../../src/plotRadiiHist.jl");
include("../../src/histmatrix.jl");
include("../../src/plotRadiiHistMatrix.jl");
include("../../src/plotQuantiles.jl");
include("../../src/plotSamplesIBs.jl");
include("../../src/plotSamplesLpdfs.jl");
#include("../../src/plotSamplesSV.jl");
plotMapIBs(outFile,trueSamp=trueSamp);
plotRadiiQuantiles(outFile, margin=10mm, trueSamp=trueSamp);
plotRadiiCorr(outFile, left_margin=10mm, bottom_margin=10mm);
#plotQuantiles(outFile, targetData="obsMean", margin=10mm);
#plotQuantiles(outFile, targetData="obsMean", targetLabel="Data", margin=10mm, markershape=:circle, xlabel="Sensor", ylabel="Observation",exts=["png","pdf"]);
plotQuantiles(outFile, targetData="obsMean", targetLabel="Data", margin=10mm, markershape=:circle, xlabel="Sensor", ylabel="Observation",exts=["png","pdf"],leg=:outerright);
plotSamplesIBs(outFile);
plotSamplesLpdfs(outFile, margin=10mm);
#plotSamplesSV(outFile, margin=10mm);

plotRadiiHist(outFile);
plotRadiiHistMatrix(outFile, left_margin=10mm, size=(1600,1600));

include("plotMap.jl");
include("plotSample.jl");
plotMap(outFile;lpdfIdx=3,circleCenters=circleCenters);
plotMap(outFile;lpdfIdx=2,circleCenters=circleCenters);

include("../../src/plotSampleGrid.jl");
#obsMean = h5read(outFile,"obsMean");
#circleColors(val) = (val==30) ? :red : ( (val==40) ? :yellow : :green );
#circleColors(val) = (val==30) ? :black : ( (val==40) ? :gray : :white );
function circleColors(val) 
  val_max = maximum(obsMean);
  val_min = minimum(obsMean);
  val_rescale = (val-val_min)/(val_max-val_min);# + val_min;
  return RGB{Float64}(val_rescale, val_rescale, val_rescale);
end
plotSampleGrid(outFile;computeScalar=false,ω=omega,circleCenters=circleCenters,circleColors=circleColors.(obsMean),sourceXY=[]);

include("../../src/plotSamplesVsTruth.jl");
plotSamplesVsTruth(outFile;circleCenters=circleCenters,exts=["png","pdf"],obsXlab="Sensor ID",obsYlab="Average Vorticity");
