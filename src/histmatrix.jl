using HDF5, Plots
import Plots: _cycle

#Define recipe for 2d histogram (without NANs for zero values)
#@recipe function f(::Type{Val{:bins2d}}, x, y, z)
#    edge_x, edge_y, weights = x, y, z.surf
#
#    float_weights = float(weights)
#    if is(float_weights, weights)
#        float_weights = deepcopy(float_weights)
#    end
#    #for (i, c) in enumerate(float_weights)
#    #    if c == 0
#    #        float_weights[i] = NaN
#    #    end
#    #end
#
#    x := Plots._bin_centers(edge_x)
#    y := Plots._bin_centers(edge_y)
#    z := Surface(float_weights)
#
#    match_dimensions := true
#    seriestype := :heatmap
#    ()
#end
#Plots.@deps bins2d heatmap
@recipe function f(::Type{Val{:bins2d}}, x, y, z)
    edge_x, edge_y, weights = x, y, z.surf

    float_weights = float(weights)
    if float_weights === weights
        float_weights = deepcopy(float_weights)
    end
    #for (i, c) in enumerate(float_weights)
    #    if c == 0
    #        float_weights[i] = NaN
    #    end
    #end

    x := Plots._bin_centers(edge_x)
    y := Plots._bin_centers(edge_y)
    z := Surface(float_weights)

    match_dimensions := true
    seriestype := :heatmap
    ()
end
Plots.@deps bins2d heatmap



#Define HistMatrix recipe
@userplot HistMatrix

recipetype(::Val{:histmatrix}, args...) = HistMatrix(args)

@recipe function f(h::HistMatrix)
  if !(typeof(h.args[1]) <: AbstractMatrix)
    error("Marginal Histograms should be given a matrix  Got: $(typeof(h.args[1]))")
  end

  m = h.args[1]
  #m, rng, idx = h.args

  #min_x = rng[1]; max_x = rng[2];
  #min_y = rng[3]; max_y = rng[4];

  idx = 1:size(m,2);
  nIdx = length(idx);  

  layout := (nIdx^2);
  
  legend := false
#  ticks  := nothing

  labs = pop!(plotattributes, :label, [""])

  for i = idx
    for j = idx
      subplot := (nIdx-j)*nIdx + i
      plotattributes[:xguide] = (j==1 ? _cycle(labs,i) : "");
      plotattributes[:yguide] = (i==1 ? _cycle(labs,j) : "");

      #1d histograms along diagonal
      if i == j
        @series begin
          link := :x
          seriestype := :stephist
          #normed := true
          normalize := true
          xformatter --> ((j == 1) ? :auto : (x -> ""))
          yformatter --> (y -> "")
          m[:,i]
        end
      #2d histograms off-diagonal
      else
        @series begin
          link := :both
          aspect_ratio := :equal
          seriescolor := :viridis
          seriestype := :histogram2d
          xformatter --> ((j == 1) ? :auto : (x -> ""))
          yformatter --> ((i == 1) ? :auto : (y -> ""))
          m[:,i], m[:,j]
        end
      end #end if
    end
  end #end for
end #end recipe


