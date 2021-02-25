function saveFEMasVTK( filename, 
                       x, 
                       eConn, 
                       scalarLabels, scalars, 
                       vectorLabels, vectors )
#  Saves a finite element solution using the WriteVTK package
#

#  using WriteVTK


  #  Get mesh dimensions
  nElements, nNodesPerElement = size(eConn)
  nNodes, nDim = size(x)
  
  #  Place the element connectivity in the correct format
  if (nNodesPerElement==3)
    cellType = VTKCellTypes.VTK_TRIANGLE       # =5 in the VTK 4.2 format
  elseif ( nNodesPerElement==6 )
    cellType = VTKCellTypes.VTK_QUADRATIC_TRIANGLE    # =22 in the VTK 4.2 format
  end

  cells = MeshCell[]
  for i=1:nElements
    c = MeshCell(cellType,eConn[i,:])
    cells = push!(cells, c)
  end


  #  Open the VTK file and write out the geometry data
  if ( nDim==2 )
    vtk = vtk_grid(filename, x[:,1], x[:,2], cells)
  elseif ( nDim==3 )
    vtk = vtk_grid(filename, x[:,1], x[:,2], x[:,3], cells)
  else
    # warning about mesh
  end

  #  Write out scalar data
  nScalars = length(scalarLabels)
  for i=1:nScalars
    vtk_point_data(vtk, scalars[:,i], scalarLabels[i])
  end

  #  Write out vector data
  nVectors = length(vectorLabels)
  if ( nDim==2 )
    for i=1:nVectors
      i1,i2 = 2*(i-1)+1, 2*(i-1)+2
      vec = (vectors[:,i1],vectors[:,i2])
      vtk_point_data(vtk, vec, vectorLabels[i])
    end
  else
    for i=1:nVectors
      i1,i2,i3 = 3*(i-1)+1, 3*(i-1)+2, 3*(i-1)+3
      vec = (vectors[:,i1],vectors[:,i2],vectors[:,i3])
      vtk_point_data(vtk, vec, vectorLabels[i])
    end
  end


  #  Close the file
  outfiles = String[]
  append!(outfiles, vtk_save(vtk))


end

  
