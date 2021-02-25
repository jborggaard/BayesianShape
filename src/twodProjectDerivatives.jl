function twodProjectDerivatives(x,eConn,u,node)
#---------------------------------------------------------------------------78--
#  twodProjectDerivatives.m - routine to project the derivative of scalar fields
#                   onto the continuous finite element space (a ZZ projection).
#
#  Copyright (c) 2020, Jeff Borggaard, Virginia Tech
#  Version: 1.0
#
#  Usage:
#     [dudx1, dudx2, e_error, node] = twodProjectDerivatives(x, e_conn, u, node)
#
#                 Multiple scalar fields can be treated at once
#                     u = [u1  u2  ... ]
#
#                  in which case
#                     dudx1 = [d(s1)/d(x1)  d(s2)/d(x1)  ... ]
#                     dudx2 = [d(s1)/d(x1)  d(s2)/d(x1)  ... ]
#
#  Output Variables:
#                 dudx1   - Projection of the x-derivatives evaluated at nodes
#                 dudx2   - Projection of the y-derivatives evaluated at nodes
#
#                 e_error - element error (H1-seminorm for each field)
#                 node    - a structure that contains an element list for
#                           each node
#  Input Variables:
#                 x       - nodal coordinates
#                 e_conn  - element connectivity
#                 u       - nodal values of scalar quantity
#
#                 node    - (optional, future development)
#                           the node structure can be reused for multiple
#                            projections involving the same mesh, e.g. for
#                            time-dependent problems
## -----------------------------------------------------------------------------

  nNodes = size(x,1)
  nElements = size(eConn,1)

  nFields = size(u,2)

  #-----------------------------------------------------------------------------
  #  For every vertex node, construct a list of elements that share it
  #-----------------------------------------------------------------------------
  nodeList = zeros(UInt64,nNodes,15)   # the first entry is reserved as a counter
                                       # we assume <= 14 elements per vertex

  for n=1:nElements  # local nodes 1,2, and 3 are vertices
    nodeList[eConn[n,1],1] = nodeList[eConn[n,1],1] + 1
    nodeList[eConn[n,1],nodeList[eConn[n,1],1]+1] = n

    nodeList[eConn[n,2],1] = nodeList[eConn[n,2],1] + 1
    nodeList[eConn[n,2],nodeList[eConn[n,2],1]+1] = n

    nodeList[eConn[n,3],1] = nodeList[eConn[n,3],1] + 1
    nodeList[eConn[n,3],nodeList[eConn[n,3],1]+1] = n
  end

  #-----------------------------------------------------------------------------
  #  Compute the least-squares projection over each patch,
  #     least-squares of gradients are performed at element quadrature points
  #-----------------------------------------------------------------------------
  nQuadrature = 7;
  r,s,w = twodQuadratureRule(nQuadrature);

  d1_p = fill(Inf,nNodes,nFields)
  d2_p = fill(Inf,nNodes,nFields)  # initialize to ridiculous values

  for n=1:nNodes
    if ( nodeList[n,1]>0 ) # this is a vertex node
      nel = nodeList[n,1]
      P   = zeros(nel*nQuadrature,6);
      d1  = zeros(nel*nQuadrature,nFields);
      d2  = zeros(nel*nQuadrature,nFields);

      for j=1:nodeList[n,1] # loop over elements in the patch
        el = nodeList[n,j+1];
        nodesLocal = eConn[el,:];
        xLocal = x[nodesLocal,:];
        uLocal = u[nodesLocal,:];

        x_g, w_g, ϕ, ϕ_x, ϕ_y = twodShape(xLocal, r, s, w);

        d_x = ϕ_x*uLocal;   # calculate the finite element derivatives
        d_y = ϕ_y*uLocal;   # on this element

        idx = 1+(j-1)*nQuadrature:j*nQuadrature;
        P[idx,:] = [ ones(nQuadrature,1)  x_g[:,1] x_g[:,2]  x_g[:,1].^2 x_g[:,1].*x_g[:,2] x_g[:,2].^2 ];
        d1[idx,:] = d_x;
        d2[idx,:] = d_y;

#       for i=1:n_quadrature
#         P  = [ P ; ...
#                1 x_g(i,1) x_g(i,2) x_g(i,1)^2 x_g(i,1)*x_g(i,2) x_g(i,2)^2 ];
#         d1 = [ d1; d_x(i) ];
#         d2 = [ d2; d_y(i) ];
#       end
      end # finished filling least-squares system

      # calculate the polynomial coefficients for each field
      for i=1:nFields
        a_x = P\d1[:,i];
        a_y = P\d2[:,i];

        #-------------------------------------------------------------------78--
        #  Compute projected derivatives at vertex node
        #-----------------------------------------------------------------------
        xp  = x[n,1];
        yp  = x[n,2];
        d1_p[n,i] = a_x[1]      + a_x[2]*xp    + a_x[3]*yp +
                    a_x[4]*xp^2 + a_x[5]*xp*yp + a_x[6]*yp^2;
        d2_p[n,i] = a_y[1]      + a_y[2]*xp    + a_y[3]*yp +
                    a_y[4]*xp^2 + a_y[5]*xp*yp + a_y[6]*yp^2;

        for j=1:nodeList[n,1]
          el = nodeList[n,j+1];
          nodesLocal = eConn[el,:];
          xLocal = x[nodesLocal,:];

          #---------------------------------------------------------------------
          #  Compute contribution to edge nodes.  These are either on the
          #  boundary or are averaged over two element faces.
          #---------------------------------------------------------------------
          if ( n==nodesLocal[1] )
            x4 = xLocal[4,1];  y4 = xLocal[4,2];
            x6 = xLocal[6,1];  y6 = xLocal[6,2];

            if ( d1_p[nodesLocal[4],i] == Inf ) # contribution not rec'd
              d1_p[nodesLocal[4],i] = 
                  a_x[1]      + a_x[2]*x4    + a_x[3]*y4 + 
                  a_x[4]*x4^2 + a_x[5]*x4*y4 + a_x[6]*y4^2;
              d2_p[nodesLocal[4],i] = 
                  a_y[1]      + a_y[2]*x4    + a_y[3]*y4 + 
                  a_y[4]*x4^2 + a_y[5]*x4*y4 + a_y[6]*y4^2;
            else
              d1_p[nodesLocal[4],i] = .5*( d1_p[nodesLocal[4],i] + 
                  a_x[1]      + a_x[2]*x4    + a_x[3]*y4 + 
                  a_x[4]*x4^2 + a_x[5]*x4*y4 + a_x[6]*y4^2 );
              d2_p[nodesLocal[4],i] = .5*( d2_p[nodesLocal[4],i] + 
                  a_y[1]      + a_y[2]*x4    + a_y[3]*y4 + 
                  a_y[4]*x4^2 + a_y[5]*x4*y4 + a_y[6]*y4^2 );
            end
            if ( d1_p[nodesLocal[6],i] == Inf )
              d1_p[nodesLocal[6],i] = 
                  a_x[1]      + a_x[2]*x6    + a_x[3]*y6 + 
                  a_x[4]*x6^2 + a_x[5]*x6*y6 + a_x[6]*y6^2;
              d2_p[nodesLocal[6],i] = 
                  a_y[1]      + a_y[2]*x6    + a_y[3]*y6 + 
                  a_y[4]*x6^2 + a_y[5]*x6*y6 + a_y[6]*y6^2;
            else
              d1_p[nodesLocal[6],i] = .5*( d1_p[nodesLocal[6],i] + 
                  a_x[1]      + a_x[2]*x6    + a_x[3]*y6 + 
                  a_x[4]*x6^2 + a_x[5]*x6*y6 + a_x[6]*y6^2 );
              d2_p[nodesLocal[6],i] = .5*( d2_p[nodesLocal[6],i] + 
                  a_y[1]      + a_y[2]*x6    + a_y[3]*y6 + 
                  a_y[4]*x6^2 + a_y[5]*x6*y6 + a_y[6]*y6^2 );
            end

          elseif ( n==nodesLocal[2] )
            x4 = xLocal[4,1];  y4 = xLocal[4,2];
            x5 = xLocal[5,1];  y5 = xLocal[5,2];

            if ( d1_p[nodesLocal[4],i] == Inf )
              d1_p[nodesLocal[4],i] = 
                  a_x[1]      + a_x[2]*x4    + a_x[3]*y4 + 
                  a_x[4]*x4^2 + a_x[5]*x4*y4 + a_x[6]*y4^2;
              d2_p[nodesLocal[4],i] = 
                  a_y[1]      + a_y[2]*x4    + a_y[3]*y4 + 
                  a_y[4]*x4^2 + a_y[5]*x4*y4 + a_y[6]*y4^2;
            else
              d1_p[nodesLocal[4],i] = .5*( d1_p[nodesLocal[4],i] + 
                  a_x[1]      + a_x[2]*x4    + a_x[3]*y4 + 
                  a_x[4]*x4^2 + a_x[5]*x4*y4 + a_x[6]*y4^2 );
              d2_p[nodesLocal[4],i] = .5*( d2_p[nodesLocal[4],i] + 
                  a_y[1]      + a_y[2]*x4    + a_y[3]*y4 + 
                  a_y[4]*x4^2 + a_y[5]*x4*y4 + a_y[6]*y4^2 );
            end
            if ( d1_p[nodesLocal[5],i] == Inf )
              d1_p[nodesLocal[5],i] = 
                  a_x[1]      + a_x[2]*x5    + a_x[3]*y5 + 
                  a_x[4]*x5^2 + a_x[5]*x5*y5 + a_x[6]*y5^2;
              d2_p[nodesLocal[5],i] = 
                  a_y[1]      + a_y[2]*x5    + a_y[3]*y5 + 
                  a_y[4]*x5^2 + a_y[5]*x5*y5 + a_y[6]*y5^2;
            else
              d1_p[nodesLocal[5],i] = .5*( d1_p[nodesLocal[5],i] + 
                  a_x[1]      + a_x[2]*x5    + a_x[3]*y5 + 
                  a_x[4]*x5^2 + a_x[5]*x5*y5 + a_x[6]*y5^2 );
              d2_p[nodesLocal[5],i] = .5*( d2_p[nodesLocal[5],i] + 
                  a_y[1]      + a_y[2]*x5    + a_y[3]*y5 + 
                  a_y[4]*x5^2 + a_y[5]*x5*y5 + a_y[6]*y5^2 );
            end

          elseif ( n==nodesLocal[3] )
            x5 = xLocal[5,1];  y5 = xLocal[5,2];
            x6 = xLocal[6,1];  y6 = xLocal[6,2];

            if ( d1_p[nodesLocal[5],i] == Inf )
              d1_p[nodesLocal[5],i] = 
                  a_x[1]      + a_x[2]*x5    + a_x[3]*y5 + 
                  a_x[4]*x5^2 + a_x[5]*x5*y5 + a_x[6]*y5^2;
              d2_p[nodesLocal[5],i] = 
                  a_y[1]      + a_y[2]*x5    + a_y[3]*y5 + 
                  a_y[4]*x5^2 + a_y[5]*x5*y5 + a_y[6]*y5^2;
            else
              d1_p[nodesLocal[5],i] = .5*( d1_p[nodesLocal[5],i] + 
                  a_x[1]      + a_x[2]*x5    + a_x[3]*y5 + 
                  a_x[4]*x5^2 + a_x[5]*x5*y5 + a_x[6]*y5^2 );
              d2_p[nodesLocal[5],i] = .5*( d2_p[nodesLocal[5],i] + 
                  a_y[1]      + a_y[2]*x5    + a_y[3]*y5 + 
                  a_y[4]*x5^2 + a_y[5]*x5*y5 + a_y[6]*y5^2 );
            end
            if ( d1_p[nodesLocal[6],i] == Inf )
              d1_p[nodesLocal[6],i] = 
                  a_x[1]      + a_x[2]*x6    + a_x[3]*y6 + 
                  a_x[4]*x6^2 + a_x[5]*x6*y6 + a_x[6]*y6^2;
              d2_p[nodesLocal[6],i] = 
                  a_y[1]      + a_y[2]*x6    + a_y[3]*y6 + 
                  a_y[4]*x6^2 + a_y[5]*x6*y6 + a_y[6]*y6^2;
            else
              d1_p[nodesLocal[6],i] = .5*( d1_p[nodesLocal[6],i] + 
                  a_x[1]      + a_x[2]*x6    + a_x[3]*y6 + 
                  a_x[4]*x6^2 + a_x[5]*x6*y6 + a_x[6]*y6^2 );
              d2_p[nodesLocal[6],i] = .5*( d2_p[nodesLocal[6],i] + 
                  a_y[1]      + a_y[2]*x6    + a_y[3]*y6 + 
                  a_y[4]*x6^2 + a_y[5]*x6*y6 + a_y[6]*y6^2 );
            end
          else
            fprintf("Problem identifying node at a vertex")
          end

        end
      end
    end
  end
#  Compare d1_p with the true gradient
#d_true=-2*x(:,1).*(1-x(:,2).^2); % elliptic_2da example
#disp([d1_p d_true])

#  Compare d2_p with the true gradient
#d_true=-2*x(:,2).*(1-x(:,1).^2);
#disp([d2_p d_true])
  computeElementError = false
  if ( computeElementError )
    #---------------------------------------------------------------------------
    #  Calculate the H1-seminorm of the error on each element
    #---------------------------------------------------------------------------
    eError = zeros(Float64,nElements,nFields);

    nQuadrature = 7;
    r,s,w   = twodQuadratureRule(nQuadrature);

    for nEl=1:nElements
      nodesLocal       = eConn[nEl,:];
      xLocal           = x[nodesLocal,:];
      uLocal           = u[nodesLocal,:];
      d1Local          = d1_p[nodesLocal,:];
      d2Local          = d2_p[nodesLocal,:];

      x_g,w_g,phi,p_x,p_y = twodShape(xLocal,r,s,w);

      # calculate the finite element derivatives and their projections at x_g
      dx_fe           = p_x*uLocal;
      dy_fe           = p_y*uLocal;

      dx_p            = phi*d1Local;
      dy_p            = phi*d2Local;

  # % exact derivatives as a test
  # dx_p = -2*x_g(:,1).*(1-x_g(:,2).^2);
  # dy_p = -2*x_g(:,2).*(1-x_g(:,1).^2);
      for i=1:nFields
        sqError     = (dx_fe[:,i]-dx_p[:,i]).^2 + (dy_fe[:,i]-dy_p[:,i]).^2;
        eError[nEl,i] = sqError'*w_g;

      end

    eError = sqrt.(eError);
    end
  end

  return d1_p, d2_p #, eError, node
end
