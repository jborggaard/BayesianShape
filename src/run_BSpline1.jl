using OffsetArrays
#using CircularArrays
#include("CircularArrays.jl")
include("BSpline.jl")

#t = range(0.02*pi,2*pi,length=100);
#Ot = CircularArray(t);
t = range(0,2*pi,length=101);
Ot = OffsetArray(t,0:100);

#knots = range(0.2*pi,2*pi,length=10);
#Oknots = CircularArray(knots);
knots = range(0,2*pi,length=11);
Oknots = OffsetArray(knots,0:10);

B01 = BSpline(Oknots,0,1,Ot);
B11 = BSpline(Oknots,1,1,Ot);
B21 = BSpline(Oknots,2,1,Ot);
B31 = BSpline(Oknots,3,1,Ot);

using Plots
plot(Ot,B01)
plot!(Ot,B11)
plot!(Ot,B21)
plot!(Ot,B31)
