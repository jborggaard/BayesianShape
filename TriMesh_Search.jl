function TriMesh_Search( x, eConn, eAdjacency, point ; verbose=0, startEl=rand(1:size(eConn,1)))
#=
   TriMesh_Search - Locates an element that contains a given point.  
   The function provides only one such element and makes an
   explicit assumption that the triangulation is Delauney.
 
   Usage:
     e, r, s = TriMesh_Search(x,eConn,eAdjacency,point)
 
   where:
     x    is an array of the nodal coordinates.
          Dimension: [ nNodes, 2 ]

     eConn  is the element connectivity.  The standard orientation is assumed
            so the first three nodes describe the vertices.
            Dimension: [ nElements, nLocalNoes ]

     eAdjacency  an [nElements, 3] array containing adjacent elements or -1
                 where -1 indicates a boundary (no adjacent element).

     point  coordinates of a point of interest in the mesh.


     e      the index of the element where "point" is located.

     r,s    the finite element coordinate of "point" in the reference element.

     startEl (optional) the element with which to start the search
   
     verbose (optional) set to >0 to increase verbosity
   
   Licensing:

     This code is distributed under the MIT license.

   Modified:

     18 January 2021

   Author:

     Jeff Borggaard
     
=#

  nElements = size(eConn,1)

  maxIterations=sqrt(nElements)
  iter = 0
  totalIter = 0

  element = startEl #define here so we can return it

  count = 0
  found = false
  iso1 = 0.0
  iso2 = 0.0
  iso3 = 0.0
  while ~found
    (totalIter > 0) && (element = rand(1:nElements));   # randomize initial guess if this is not our first iteration
    (verbose>0) && println("TriMesh_Search: Starting with element $(element)...");
    elementOld = element
    iter=0

    while ~found && iter<maxIterations
      iter += 1
      (verbose>1) && println("TriMesh_Search: Iteration $(iter) with element $(element)")
      xLocal = x[ eConn[element,1:3], : ]    # view?
    
#    plot( [ xLocal[1,1] xLocal[2,1] ],[ xLocal[1,2] xLocal[2,2]],...
#          [ xLocal[2,1] xLocal[3,1] ],[ xLocal[2,2] xLocal[3,2]],...
#          [ xLocal[3,1] xLocal[1,1] ],[ xLocal[3,2] xLocal[1,2]] )
      
      delta = ( xLocal[3,1]-xLocal[2,1] )*( xLocal[1,2]-xLocal[2,2] ) -
              ( xLocal[3,2]-xLocal[2,2] )*( xLocal[1,1]-xLocal[2,1] )
        
      iso1  = ( ( xLocal[2,1]-point[1] )*( xLocal[3,2]-point[2] ) -
                ( xLocal[3,1]-point[1] )*( xLocal[2,2]-point[2] ) ) / delta;
    
      iso2  = ( ( xLocal[3,1]-point[1] )*( xLocal[1,2]-point[2] ) -
                ( xLocal[1,1]-point[1] )*( xLocal[3,2]-point[2] ) ) / delta;
    
      iso3  = ( ( xLocal[1,1]-point[1] )*( xLocal[2,2]-point[2] ) -
                ( xLocal[2,1]-point[1] )*( xLocal[1,2]-point[2] ) ) / delta;
      
      if iso1>=0 && iso2>=0 && iso3>=0
        (verbose>0) && println("TriMesh_Search: FOUND at iteration $(iter)")
        found = true
      
      else
        count = count + 1;
        index = sortperm([iso1, iso2, iso3])
        
        if eAdjacency[element,index[1]]>0
          element_old = element;
          element = eAdjacency[element,index[1]]
          choice = 1
        elseif eAdjacency[element,index[2]]>0
          element_old = element
          element = eAdjacency[element,index[2]]
          choice = 2
        elseif eAdjacency[element,index[3]]>0
          element_old = element
          element = eAdjacency[element,index[3]]
          choice = 3;
        else
           # iso1
           # iso2
           # iso3
           # eAdjacency[element,:]
           # index
           error("no next element found, Delaunay mesh?")
        end
#       if ( index(choice)==1 )
#         plot( (xLocal(2,1)+xLocal(3,1))/2, (xLocal(2,2)+xLocal(3,2))/2, '*')
#       elseif ( index(choice)==2 )
#         plot( (xLocal(1,1)+xLocal(3,1))/2, (xLocal(1,2)+xLocal(3,2))/2, '*')
#       else
#         plot( (xLocal(1,1)+xLocal(2,1))/2, (xLocal(1,2)+xLocal(2,2))/2, '*')
#       end
        if count>nElements   # if we are somehow in a loop, get out...
          element = nElements
          count   = 0
        end
      end
    end
    totalIter += iter;
  end
  #(iter>=maxIterations) && @warn("Reached maximum number of iterations ($(maxIterations)).");
  (~found) && @warn("Element may not have been found. (iter=$(iter)).");
  (verbose>0) && println("TriMesh_Search: Total number of iterations: $(totalIter)");
  
  return element, iso2, iso3

end
