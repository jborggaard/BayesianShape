#project a triangle with given vertices onto fourier components
#note: here the first (constant) fourier component is assumed to be r=1/sqrt(2)
function triangleProject( vert, nodes, weights; nUnk=160)
  #build triangle (quadrature points used on each face)
  th, r = polarTriangle(vert,nodes);
  
  #rescale weights
  w = weights * ( th[end,:] - th[1,:] )' ./ ( nodes[end] - nodes[1] );
  
  #collect into vectors
  thall = th[:]; rall = r[:]; wall = w[:];
  
  #basis (rescaled to be norm 1 on the circle)
  fb = zeros( length(thall), 2*nUnk+1 );
  fb[:,1] .= 1.0/sqrt(2*pi);
  fb[:,2:2:end] = cos.(thall*(1:nUnk)')/sqrt(pi);
  fb[:,3:2:end] = sin.(thall*(1:nUnk)')/sqrt(pi);
  
  #project
  proj = fb' * (wall .* rall);
  
  ##reconstruct
  #rproj = fb * proj;

  return proj;
end
