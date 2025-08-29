#using HDF5, Plots
import Plots: _cycle

# NOTE: This section used to not be commented, but it was causing precompilation errors like so:
#   WARNING: Method definition apply_recipe(Base.AbstractDict{Symbol, Any}, Type{Base.Val{:bins2d}}, Any, Any, Any) in module Plots at /home/jkrometi/.julia/packages/RecipesBase/BRe07/src/RecipesBase.jl:296 overwritten in module BayesianShape on the same line (check for duplicate calls to `include`).
#   ERROR: Method overwriting is not permitted during Module precompilation. Use `__precompile__(false)` to opt-out of precompilation.
# I believe this is because I was attempting to overwrite the meaning of bins2d (it's been a while since I wrote this).
# If we want to reenable it, I think we need to define a custom type and call that. See https://docs.juliaplots.org/dev/recipes/
#
# @recipe function f(::Type{Val{:bins2d}}, x, y, z)
#     edge_x, edge_y, weights = x, y, z.surf
# 
#     float_weights = float(weights)
#     if float_weights === weights
#         float_weights = deepcopy(float_weights)
#     end
#     #for (i, c) in enumerate(float_weights)
#     #    if c == 0
#     #        float_weights[i] = NaN
#     #    end
#     #end
# 
#     x := Plots._bin_centers(edge_x)
#     y := Plots._bin_centers(edge_y)
#     z := Surface(float_weights)
# 
#     match_dimensions := true
#     seriestype := :heatmap
#     ()
# end
# Plots.@deps bins2d heatmap



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


