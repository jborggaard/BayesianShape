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

B00 = BSpline(Oknots,0,0,Ot);
B10 = BSpline(Oknots,1,0,Ot);
B20 = BSpline(Oknots,2,0,Ot);
B30 = BSpline(Oknots,3,0,Ot);
B40 = BSpline(Oknots,4,0,Ot);
B50 = BSpline(Oknots,5,0,Ot);

using Plots
plot(Ot,B00)
plot!(Ot,B10)
plot!(Ot,B20)
plot!(Ot,B30)
plot!(Ot,B40)
plot!(Ot,B50)

