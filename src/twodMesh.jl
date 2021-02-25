function twodMesh(x_l::Float64, x_r::Float64,
                  y_l::Float64, y_r::Float64,
                  etype::String,
                  n_nodesx::Int64, n_nodesy::Int64)
##------------------------------------------------------------------------------
#  twodMesh.m - Generate a rectangular mesh with a prescribed density.
#               This routine returns nodal coordinates, element connectivity,
#               and the nodal indices of exterior nodes to assist in the
#               specification of boundary conditions..
#
#  Copyright (c) 2001, Jeff Borggaard, Virginia Tech
#  Version: 1.0
#
#  Usage:  x, e_conn, index_b = twodMesh(x_l,x_r,y_l,y_r,etype,n_nodex,n_nodesy)
#
#  Variables:     (x_l,y_l)
#                        coordinates of lower left corner
#                 (x_r,y_r)
#                        coordinates of upper right corner
#
#                 etype
#                        element type
#                          'linear', 'quadratic', 'cubic'
#
#                 n_nodesx, n_nodesy
#                        number of nodes in x and y directions
#                        (must be compatible with element type)
#
#    Outputs:
#                 x
#                        node coordinates of mesh
#                 e_conn
#                        element connectivity of mesh
#                 index_b
#                        node numbers of boundary nodes
#-------------------------------------------------------------------------------

  #  Generate node coordinates
  nNodes = n_nodesx*n_nodesy;
  x = zeros(Float64,nNodes,2);
  dx = (x_r-x_l)/(n_nodesx-1);
  dy = (y_r-y_l)/(n_nodesy-1);
  for i=1:n_nodesx
    for j=1:n_nodesy
      k = i+(j-1)*n_nodesx;
      x[k,1] = x_l + dx*(i-1);
      x[k,2] = y_l + dy*(j-1);
    end
  end

  #  Generate element connectivity
  if etype == "linear"
    n_elements = 2*(n_nodesx-1)*(n_nodesy-1);

    global e_conn = zeros(Int64, n_elements, 3);
    ie     = 0;
    for j=1:n_nodesy-1
      for i=1:n_nodesx-1
        ie = ie + 1;
        k = i + (j-1)*n_nodesx;
        e_conn[ie,:] = [k, k+1+n_nodesx, k+n_nodesx];
        ie = ie + 1;
        e_conn[ie,:] = [k, k+1, k+1+n_nodesx];
      end
    end

  elseif etype == "quadratic"
    n_elements = convert(Int,2*(n_nodesx-1)*(n_nodesy-1)/4);

    e_conn = zeros(Int64, n_elements, 6);
    ie     = 0;
    for j=1:2:n_nodesy-1
      for i=1:2:n_nodesx-1
        ie = ie + 1;
        k = i + (j-1)*n_nodesx;
        e_conn[ie,:] = [k, k+2+2*n_nodesx, k+2*n_nodesx,
                        k + 1 + n_nodesx, k+1+2*n_nodesx, k+n_nodesx ];
        ie = ie + 1;
        e_conn[ie,:] = [k, k+2, k+2+2*n_nodesx,
                        k+1, k+2+n_nodesx, k+1+n_nodesx ];
      end
    end

  elseif etype == "cubic"
    @printf("twod_mesh: cubic elements are not currently implemented\n")
    return
  else
    @printf("twod_mesh: etype is not a valid string\n")
    return
  end

  index_b = [ 2:n_nodesx-1
              1:n_nodesx:(n_nodesy-1)*n_nodesx
              n_nodesx:n_nodesx:n_nodesy*n_nodesx
              (n_nodesy-1)*n_nodesx+1:nNodes-1];

  return x, e_conn, index_b;

end
