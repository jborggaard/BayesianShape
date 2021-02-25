#getMap() find a map/mle point from an output file and plot it
# inFile    input file
# lpdfIdx   3 (default) for map, 2 for mle, 1 for max prior, 0 for array of the three
#
function getMap(samples,lpdfs,lpdfIdx=3)
  #find map/mle point
  #idx = argmax(lpdfs,dims=1)[lpdfIdx][1];
  if lpdfIdx == 0
    idx = [ x[1] for x=argmax(lpdfs,dims=1) ][:];
  else
    idx = argmax(lpdfs[:,lpdfIdx])[1];
  end

  return samples[idx,:];
end

function getMap(inFile,lpdfIdx=3)
  f = h5open(inFile);
  lpdfs   = read(f,"lpdfs");
  samples = read(f,"samples");
  close(f);

  return getMap(samples,lpdfs,lpdfIdx);
end

