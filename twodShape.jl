function twodShape(x,r,s,w)
#----------------------------------------------------------------------------78-
#  twodShape - computes test functions and derivatives on an
#              element given element coordinates and Gauss points.
#
#              ! Note: optimized for straight-sided elements.  Use
#              ! `twod_shapeiso' for isoparametric elements.
#
#              local unknown numbering follows the conventions
#
#            3             3             3             3
#            |\            |\            |\            87
#            | \     or    6 5     or    6 5     or    | \   x is 10
#            |  \          |  \          |7 \          9x 6
#            1---2         1-4-2         1-4-2         14-52
#
#             P1            P2      Crouzeux-Raviart    P3
#
#
#  Copyright (c) 2015, Jeff Borggaard, Virginia Tech
#  Version: 1.1
#
#  Usage:    xg,wg,phi,p_x,p_y = twodShape(x,r,s,w)
#
#  Variables:     x
#                        Coordinates of the element nodes
#                 (r,s)
#                        Coordinates of Gauss points in unit triangle
#                 w
#                        Gauss weights associated with (r,s)
#
#                 xg
#                        Coordinates of Gauss points in the element
#                 wg
#                        Gauss weights scaled by the element Jacobian
#                 phi
#                        Value of element shape functions at xg
#                 p_x
#                 p_y
#                        First spatial derivatives of phi
#
#  Changes:
#      Crouzeux-Raviart element implemented by Alexander Hay, 2007
#-------------------------------------------------------------------------------
  n   = size(x,1);
  ng  = size(r,1);

  xg  = zeros(Float64,ng,2);
  phi = zeros(Float64,ng,n);
  p_x = zeros(Float64,ng,n);
  p_y = zeros(Float64,ng,n);

  wg  = zeros(Float64,ng);

  o   = ones(Float64,size(r));

  # Compute (r,s) -> (x,y) transformation for straight-sided elements
  c0 =  x[1,:];
  c1 = -x[1,:] + x[2,:];
  c2 = -x[1,:]          + x[3,:];

  xg[:,1] = c0[1]*o .+ c1[1]*r .+ c2[1]*s;
  xr      = c1[1]*o;
  xs      = c2[1]*o;

  xg[:,2] = c0[2]*o .+ c1[2]*r .+ c2[2]*s;
  yr      = c1[2]*o;
  ys      = c2[2]*o;

  # Compute the Jacobian of the (r,s) -> (x,y) transformation
  jac = xr.*ys .- yr.*xs;
  wg  = jac.*w;

  rx  = ys./jac;
  sx  =-yr./jac;
  ry  =-xs./jac;
  sy  = xr./jac;


  # Compute shape function and derivatives at Gauss points
  if n == 3
    phi[:,1] = 1.0*o - r  - s;
    phi[:,2] =         r     ;
    phi[:,3] =              s;
   
    p_x[:,1] =       -rx - sx;
    p_x[:,2] =        rx     ;
    p_x[:,3] =             sx;
   
    p_y[:,1] =      - ry - sy;
    p_y[:,2] =        ry     ;
    p_y[:,3] =             sy;

  elseif n == 6 
    phi[:,1] = 1.0*o - 3.0*r - 3.0*s + 2.0*r.*r + 4.0*r.*s + 2.0*s.*s;
    phi[:,2] =       - 1.0*r         + 2.0*r.*r                      ;
    phi[:,3] =               - 1.0*s                       + 2.0*s.*s;
    phi[:,4] =         4.0*r         - 4.0*r.*r - 4.0*r.*s           ;
    phi[:,5] =                                    4.0*r.*s           ;
    phi[:,6] =                 4.0*s            - 4.0*r.*s - 4.0*s.*s;
  
    p_x[:,1] = ( -3.0*o + 4.0*r + 4.0*s ).*rx + ( -3.0*o + 4.0*r + 4.0*s ).*sx;
    p_x[:,2] = ( -1.0*o + 4.0*r         ).*rx                                 ;
    p_x[:,3] =                                  ( -1.0*o         + 4.0*s ).*sx;
    p_x[:,4] = (  4.0*o - 8.0*r - 4.0*s ).*rx + (        - 4.0*r         ).*sx;
    p_x[:,5] = (                  4.0*s ).*rx + (          4.0*r         ).*sx;
    p_x[:,6] = (                - 4.0*s ).*rx + (  4.0*o - 4.0*r - 8.0*s ).*sx;
   
    p_y[:,1] = ( -3.0*o + 4.0*r + 4.0*s ).*ry + ( -3.0*o + 4.0*r + 4.0*s ).*sy;
    p_y[:,2] = ( -1.0*o + 4.0*r         ).*ry                                 ;
    p_y[:,3] =                                  ( -1.0*o         + 4.0*s ).*sy;
    p_y[:,4] = (  4.0*o - 8.0*r - 4.0*s ).*ry + (        - 4.0*r         ).*sy;
    p_y[:,5] = (                  4.0*s ).*ry + (          4.0*r         ).*sy;
    p_y[:,6] = (                - 4.0*s ).*ry + (  4.0*o - 4.0*r - 8.0*s ).*sy;
    
  elseif n == 7
    phi[:,1] = (1.0*o-r-s).*(2.0*(1.0-r-s)-1.0) +  3.0*(1.0-r-s).*r.*s; 
    phi[:,2] = r.*(2.0*r-1.0)                 +  3.0*(1.0-r-s).*r.*s; 
    phi[:,3] = s.*(2.0*s-1.0)                 +  3.0*(1.0-r-s).*r.*s; 
    phi[:,4] = 4.0*(1.0-r-s).*r               - 12.0*(1.0-r-s).*r.*s; 
    phi[:,5] = 4.0*r.*s                       - 12.0*(1.0-r-s).*r.*s; 
    phi[:,6] = 4.0*s.*(1.0-r-s)               - 12.0*(1.0-r-s).*r.*s; 
    phi[:,7] = 27.0*(1.0-r-s).*r.*s;                   

    p_r = zeros(Float64,ng,n);
    p_r[:,1] = -3.0 + 4.0*r + 7.0*s - 6.0*r.*s - 3.0*(s.^2);
    p_r[:,2] = -1.0 + 4.0*r + 3.0*s - 6.0*r.*s - 3.0*(s.^2);
    p_r[:,3] =                3.0*s - 6.0*r.*s - 3.0*(s.^2);
    p_r[:,4] =  4.0 - 8.0*r -16.0*s +24.0*r.*s +12.0*(s.^2);
    p_r[:,5] =              - 8.0*s +24.0*r.*s +12.0*(s.^2);
    p_r[:,6] =              -16.0*s +24.0*r.*s +12.0*(s.^2);
    p_r[:,7] =               27.0*s -54.0*r.*s -27.0*(s.^2);

    p_s = zeros(Float64,ng,n);
    p_s[:,1] = -3.0 + 7.0*r + 4.0*s - 6.0*r.*s - 3.0*(r.^2);
    p_s[:,2] =        3.0*r         - 6.0*r.*s - 3.0*(r.^2);
    p_s[:,3] = -1.0 + 3.0*r + 4.0*s - 6.0*r.*s - 3.0*(r.^2);
    p_s[:,4] =      -16.0*r         +24.0*r.*s +12.0*(r.^2);
    p_s[:,5] =      - 8.0*r         +24.0*r.*s +12.0*(r.^2);
    p_s[:,6] =  4.0 -16.0*r - 8.0*s +24.0*r.*s +12.0*(r.^2);
    p_s[:,7] =       27.0*r         -54.0*r.*s -27.0*(r.^2);

    p_x = p_r.*rx + p_s.*sx;
    p_y = p_r.*ry + p_s.*sy;

  elseif ( n==10 )
    # Compute shape function and derivatives at Gauss points
    t = ones(size(r)) - r - s;
    
    phi = zeros(Float64,ng,n);
    phi[:, 1] =  4.5*(t-1/3).*(t-2/3).*t;
    phi[:, 2] =  4.5*(r-1/3).*(r-2/3).*r;
    phi[:, 3] =  4.5*(s-1/3).*(s-2/3).*s;
    phi[:, 4] =  9.0*r-22.5*r.*(r+s)+13.5*r.*(r+s).^2;
    phi[:, 5] = -4.5*r+18.0*r.^2+4.5*r.*s-13.5*r.^3-13.5*(r.^2).*s;
    phi[:, 6] = -4.5*r.*s+13.5*(r.^2).*s;
    phi[:, 7] = -4.5*r.*s+13.5*r.*s.^2;
    phi[:, 8] = -4.5*s+4.5*r.*s+18.0*s.^2-13.5*r.*s.^2-13.5*s.^3;
    phi[:, 9] =  9.0*s-22.5*r.*s-22.5*s.^2+13.5*s.*(r+s).^2;
    phi[:,10] = 27.0*r.*s-27.0*(r.^2).*s-27.0*r.*s.^2;

    p_r = zeros(Float64,ng,n);
    p_r[:, 1] = -5.5 +18.0*(r+s)-13.5*(r+s).^2;
    p_r[:, 2] =  1.0 -9.0*r+13.5*r.^2;
#   p_r[:, 3] =  zeros(Float64,ng,1);
    p_r[:, 4] =  9.0 -45.0*r-22.5*s+40.5*r.^2+54.0*r.*s+13.5*s.^2;
    p_r[:, 5] = -4.5 +36.0*r+4.5*s-40.5*r.^2-27.0*r.*s;
    p_r[:, 6] = -4.5*s+27.0*r.*s;
    p_r[:, 7] = -4.5*s+13.5*s.*s;
    p_r[:, 8] =  4.5*s-13.5*s.*s;
    p_r[:, 9] =-22.5*s+27.0*(r.*s+s.*s);
    p_r[:,10] = 27.0*(s - 2.0*r.*s-s.*s);

    p_s = zeros(Float64,ng,n);
    p_s[:, 1] = -5.5+18.0*(r+s)-13.5*(r+s).^2;
#   p_s[:, 2] = zeros(Float64,ng,1);
    p_s[:, 3] = 1.0-9.0*s+13.5*s.^2;
    p_s[:, 4] = -22.5*r+27.0*(r.^2+r.*s);
    p_s[:, 5] = 4.5*r-13.5*r.^2;
    p_s[:, 6] =-4.5*r+13.5*r.^2;
    p_s[:, 7] =-4.5*r+27.0*r.*s;
    p_s[:, 8] =-4.5+4.5*r+36.0*s-27.0*r.*s-40.5*s.^2;
    p_s[:, 9] = 9.0-22.5*r-45.0*s+13.5*r.^2+54.0*r.*s+40.5*s.^2;
    p_s[:,10] = 27.0*(r - r.^2 - 2.0*r.*s);

    p_x = p_r.*rx + p_s.*sx;
    p_y = p_r.*ry + p_s.*sy;
    
#   else
#     warning('element not supported');
  end

  return xg, wg, phi, p_x, p_y

end
