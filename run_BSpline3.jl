using Plots
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

B03 = BSpline(Oknots,0,3,Ot);
B13 = BSpline(Oknots,1,3,Ot);
B23 = BSpline(Oknots,2,3,Ot);
B33 = BSpline(Oknots,3,3,Ot);
B43 = BSpline(Oknots,4,3,Ot);
B53 = BSpline(Oknots,5,3,Ot);
B63 = BSpline(Oknots,6,3,Ot);

plot(Ot,B03)
plot!(Ot,B13)
plot!(Ot,B23)
plot!(Ot,B33)
plot!(Ot,B43)
plot!(Ot,B53)
plot!(Ot,B63)
