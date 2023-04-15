include("../../src/getMap.jl");
include("plotSample.jl");

#plotMap() find a map/mle point from an output file and plot it
# inFile    input file
# lpdfIdx   3 (default) for map, 2 for mle, 1 for max prior
# outFile   root of output file name
#
function plotMap(inFile;lpdfIdx=3,outFile="default",circleCenters=[])
#function plotMap(inFile;lpdfIdx=3,outFile=replace(inFile,".h5"=>"_map"))
  #f = h5open(inFile);
  #lpdfs   = read(f,"lpdfs");
  #samples = read(f,"samples");
  #close(f);
  #
  ##find map/mle point
  #idx = argmax(lpdfs,dims=1)[lpdfIdx][1];
  #
  ##plot
  #ps = plotSample(samples[idx,:],outFile);
  
  if outFile == "default"
    mapStr = ["mpr","mle","map"][lpdfIdx];
    outFile = replace(inFile,".h5"=>"_"*mapStr);
  end
  
  #find map/mle point
  sample = getMap(inFile,lpdfIdx);

  #read some problem parameters
  #kappa = h5read(inFile,"kappa");
  omega = h5read(inFile,"omega");
  #sourceXY = h5read(inFile,"sourceXY");
  
  #plot
  ps = plotSample(sample,outFile;ω=omega,circleCenters=circleCenters);

  return ps;
end
