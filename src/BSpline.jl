function BSpline(knots,index,order,t)
#  We assume knots, index, and t are zero-based arrays

  f = similar(Array{Float64},axes(t))

  if ( order==0 )
    for i=eachindex(t)
      if t[i]>=knots[index] && t[i]<knots[index+1]
        f[i] = 1.0;
      end
    end

    return f
  end

  f = (t.-knots[index]).*BSpline(knots,index,order-1,t)/(knots[index+order]-knots[index]) .+ 
      (knots[index+order+1].-t).*BSpline(knots,index+1,order-1,t)/(knots[index+order+1]-knots[index+1])

  return f
end
