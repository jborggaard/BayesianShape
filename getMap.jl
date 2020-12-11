#getMap() find a map/mle point from an output file and plot it
# inFile    input file
# lpdfIdx   3 (default) for map, 2 for mle, 1 for max prior
#
function getMap(inFile,lpdfIdx=3)
  f = h5open(inFile);
  lpdfs   = read(f,"lpdfs");
  samples = read(f,"samples");
  close(f);

  #find map/mle point
  #idx = argmax(lpdfs,dims=1)[lpdfIdx][1];
  idx = argmax(lpdfs[:,lpdfIdx])[1];

  return samples[idx,:];
end

