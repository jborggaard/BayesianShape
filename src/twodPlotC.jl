function twodPlotC(x::Array{Float64,2},eConn::Array{Int64,2},u::Array{Float64,1})
#  Plots a finite element solution in 2D

  if ( size(eConn,2)==3 )  # linear elements
    lConn = eConn-ones(size(eConn));
  elseif ( size(eConn,2)==6 ) # quadratic elements
    lConn = eConn[:,1:3]-ones(size(eConn[:,1:3]));
  end

  # plot_trisurf(x[:,1],x[:,2],lConn,u)
  cm = PyPlot.ColorMap("jet")
  PyPlot.plot_trisurf(x[:,1],x[:,2],lConn,u,cmap=cm,linewidth=0)

  # set_filename("test.svg")
  PyPlot.savefig("plot.svg")

  # set_filename("test.png")
  # printfigure("png")

end
