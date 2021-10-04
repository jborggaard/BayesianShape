function triangleRecenter( vert )
  #centroid
  cent = sum(vert,dims=1)/size(vert,1); #mean(vert,dims=1);
  #recenter
  return vert .- cent;
end
