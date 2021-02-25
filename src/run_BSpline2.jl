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

B02 = BSpline(Oknots,0,2,Ot);
B12 = BSpline(Oknots,1,2,Ot);
B22 = BSpline(Oknots,2,2,Ot);
B32 = BSpline(Oknots,3,2,Ot);
B42 = BSpline(Oknots,4,2,Ot);
B52 = BSpline(Oknots,5,2,Ot);
B62 = BSpline(Oknots,6,2,Ot);

plot(Ot,B02)
plot!(Ot,B12)
plot!(Ot,B22)
plot!(Ot,B32)
plot!(Ot,B42)
plot!(Ot,B52)
plot!(Ot,B62)
