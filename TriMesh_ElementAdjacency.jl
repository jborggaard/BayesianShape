function TriMesh_ElementAdjacency(eConn)
#= -------------------------------------------------------------------------78--
 
   Discussion:
 
     TriMesh_ElementAdjacency determines the adjacent elements for each 
     triangular element and for each vertex in a triangular mesh.  This is based 
     on the function TRIANGULARION_TRIANGLE_NEIGHBORS written by John Burkardt.
 
   Usage:
 
     eAdjacency = TriMesh_ElementAdjacency(eConn)
     # planned: eAdjacency,vConn = TriMesh_ElementAdjacency(eConn)
 
     where:
 
     eConn     is the element connectivity.  Standard orientation is assumed
               so the first three nodes in the connectivity are the vertices.
               Dimension: [ nElements, nLocalNodes ]
 
     eAdjacency  is a [ n_elements, 3 ] array which contains the adjacent 
                 element to each triangle edge.  A value of -1 is given if that
                 edge has no adjacent element (e.g., on the boundary)
 
                                 3
                                | \
                         edge 2 |  \ edge 1
                                | e \
                                1----2
                                edge 3
 
 
     #vConn    is a structure of n_node arrays that list elements that are
               connected to each node.  vConn[1].e_list contains a list of the
               elements that include node 1, etc.  (planned)
 
   Licensing:
 
     This code is distributed under the MIT license.
 
   Modified:
 
     18 January 2021
 
   Authors:
 
     Jeff Borggaard and John Burkardt
 
##--------------------------------------------------------------------------78=#

  nElements, nelDOF = size(eConn)
  
#
#  Detect and correct 0-based indexing.
#
# minN = minimum(eConn)
# if (minN==0)
#   @printf('TriMesh_ElementAdjacency: element connectivity may be 0-based, ')
#   @printf('indexing\n will attempt to adjust and proceed\n')
#   eConn = eConn .+ 1
# end
  
  nNodes = maximum(eConn)
  
#
#  Create the triangle neighbor array.
#
 
#  Step 1.
#  From the list of nodes for triangle T, of the form: (I,J,K)
#  construct the three neighbor relations:
#
#    (I,J,3,T) or (J,I,3,T),
#    (J,K,1,T) or (K,J,1,T),
#    (K,I,2,T) or (I,K,2,T)
#
#  where we choose (I,J,3,T) if I < J, or else (J,I,3,T)
  
  faces = Array{Int64,2}(undef,4,3*nElements)
  for nElem = 1:nElements

    i = eConn[nElem,1];
    j = eConn[nElem,2];
    k = eConn[nElem,3];

    if i<j
      faces[:,1+3*(nElem-1)] = [ i, j, 3, nElem ]
    else
      faces[:,1+3*(nElem-1)] = [ j, i, 3, nElem ]
    end

    if j<k
      faces[:,2+3*(nElem-1)] = [ j, k, 1, nElem ]
    else
      faces[:,2+3*(nElem-1)] = [ k, j, 1, nElem ]
    end

    if k<i
      faces[:,3+3*(nElem-1)] = [ k, i, 2, nElem ]
    else
      faces[:,3+3*(nElem-1)] = [ i, k, 2, nElem ]
    end

  end
 
#  Step 2. Perform an ascending dictionary sort on the adjacency relations.
#  We only intend to sort on rows 1:2; the routine we call here
#  sorts on rows 1 through 4 but that won't hurt us.
#
#  What we need is to find cases where two triangles share a face.
#  By sorting the columns of the FACES array, we will put shared faces
#  next to each other.
#
#  this scrambles the element information.  
#   faces = sort(faces,1);
   
#   faces = sortslices(faces,dims=1)?
   index2 = sortperm(faces[2,:])
   faces = faces[:,index2]
   
   index1 = sortperm(faces[1,:])
   faces = faces[:,index1]

#  Step 3. Neighboring triangles show up as consecutive columns with
#  identical first two entries.  Whenever you spot this happening,
#  make the appropriate entries in eAdjacency.

  eAdjacency = fill(-1,nElements,3)

  face = 1
  while true

    if 3nElements<=face
      break
    end

    if ( faces[1:2,face] == faces[1:2,face+1] )
      face1 = faces[3,face]
      elem1 = faces[4,face]
      face2 = faces[3,face+1]
      elem2 = faces[4,face+1]
      eAdjacency[elem1,face1] = elem2
      eAdjacency[elem2,face2] = elem1
      face = face + 2
    else # there is no neighbor, move on
      face = face + 1;
    end

  end
    
  
#  if nargout==2
#    for n=1:nNodes
#      vConn[n].eList = []
#    end
#    for n=1:nElements  % for each node, add the current element to e_list
#      for m=1:nelDOF
#        vConn[eConn[n,m]].eList = [ vConn[eConn[n,m]].eList, n ]
#      end
#    end
#  end

  return eAdjacency #,v_conn

end # function TriMesh_ElementAdjacency
