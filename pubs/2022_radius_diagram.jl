using HDF5, Printf, Plots, Plots.Measures

include("../src/fourierBasis.jl");
include("../src/computeRadii.jl");
include("../src/plotSave.jl");

#labels
rMinLab = "\$r_{min}\$";      #minimum inner boundary
rMaxLab = "\$r_{max}\$";      #maximum inner boundary
rOutLab = "\$\\Gamma^{\\rm outer}\$ (radius \$R\$)";      #outer boundary
boundaryLab = "\$\\Gamma_b^{\\rm inner}\$ (radius \$c(b_0+b)\$)";  #inner boundary

rOutLabRad = "\$R\$";             #outer boundary
boundaryLabRad = "\$c(b_0+b)\$";  #inner boundary
unsquashedLab = "\$b_0+b\$";       #unsquashed boundary

#run setup to get radius squash (uses default :/ )
#also gets a0
include("../scenarios/svsector/setup.jl");

#file = "/projects/SIAllocation/stokes/svsector/svsector_029.h5";
file = "pubs/2022_radius_diagram.h5";

if isfile(file)
  println("Reading from: $file");
  f = h5open(file,"r");
  r    = read(f,"r");
  b    = read(f,"b");
  th   = read(f,"th");
  rMin = read(f,"rMin");
  rMax = read(f,"rMax");
  rOut = read(f,"rOut");
  close(f);
else
  # #radii
  # rMin = 0.5; 
  # rMax = 1.5; 
  rOut = 2.0; 
  
  # #read from file
  # f = h5open(file);
  # samples = read(f,"samples");
  # lpdfs   = read(f,"lpdfs");
  # close(f);
  # 
  # #get index
  # idx = argmax(lpdfs[:,3]);
  # samples = samples[idx,:];
  
  #angles
  dth = pi/180;
  th=0:dth:2*pi;
  #fourier basis
  # n  = length(samples);
  # nf = isodd(n) ? (n-1)÷2 : n÷2;
  # fb = fourierBasis(nf,th);
  fb = fourierBasis(unkDim,th);

  r = zeros(2*unkDim);
  b = zeros(2*unkDim);
  samples = zeros(2*unkDim);

  stopIter = false
  cnt = 1
  while !stopIter
    global stopIter, cnt, r, b, samples #julia is annoying
    println("Attempt $cnt");
    #draw randomly from prior
    samples = rand(mcmcP.prior);
    
    #compute radii
    r = computeRadii(samples,th);
    
    #compute unsquashed radii
    b  = a0 .+ fb'*samples;

    pct_at_boundary = ( sum(r .== rMax) + sum(r .== rMin) ) /length(r);

    stopIter = (pct_at_boundary < 0.25) && (cnt <= 50); 
    cnt += 1;
  end
end

# #get index
# idx = argmax(lpdfs[:,3]);
# r = r[idx,:];

#plot (radial)
sz = 500;
p = plot(proj=:polar,size=(sz,sz), margin=5Plots.mm);
plot!(p,th,[rMin],ls=:dash,lab=rMinLab);
plot!(p,th,r,lab=boundaryLab);
plot!(p,th,[rMax],ls=:dash,lab=rMaxLab);
plot!(p,th,[rOut],c=:black,lab=rOutLab);
plotSave(p,"pubs/2022_radius_diagram");

#plot (rectangular)
#set lc to manually match above
thDeg = th*180/pi;
sz = 500;
p = plot(size=(sz,sz), margin=5Plots.mm, xlab="Angle, \$\\phi\$ (Degrees)", ylab="Radius", xticks=0:90:360);
#plot!(p,th,[rMin],lc=1,ls=:dash,lab=rMinLab);
hline!(p,[rMin],lc=1,ls=:dash,lab=rMinLab);
plot!(p,thDeg,b,lc=4,lab=unsquashedLab);
plot!(p,thDeg,r,lc=2,lab=boundaryLabRad);
#plot!(p,th,[rMax],lc=3,ls=:dash,lab=rMaxLab);
#plot!(p,th,[rOut],c=:black,lab=rOutLabRad);
hline!(p,[rMax],lc=3,ls=:dash,lab=rMaxLab);
hline!(p,[rOut],c=:black,lab=rOutLabRad);
plot!(legend=:bottomright);
plotSave(p,"pubs/2022_radius_diagram_rect");

if !isfile(file)
  f = h5open(file,"w");
  write(f,"r",r);
  write(f,"b",b);
  write(f,"samples",samples);
  write(f,"th",collect(th));
  write(f,"rMin",rMin);
  write(f,"rMax",rMax);
  write(f,"rOut",rOut);
  close(f);
  println("Wrote data to: $file");
end
