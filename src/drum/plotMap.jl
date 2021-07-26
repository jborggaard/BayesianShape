#include("getMap.jl");
#include("plotSampleShapes.jl");

#plotMap() find a map/mle point from an output file and plot it
# inFile    input file
# lpdfIdx   3 (default) for map, 2 for mle, 1 for max prior
# outFile   root of output file name
#
function plotMap(inFile;lpdfIdx=3,outFile="default",kwargs...)
  
  if outFile == "default"
    mapStr = ["mpr","mle","map"][lpdfIdx];
    outFile = replace(inFile,".h5"=>"_"*mapStr);
  end
  
  #find map/mle point
  sample = getMap(inFile,lpdfIdx);

  rMin = h5read(inFile,"rMin");
  rMax = h5read(inFile,"rMax");

  #plot
  ps = plotSampleShapes(sample,outFile;rMin=rMin,rMax=rMax,kwargs...);

  return ps;
end
