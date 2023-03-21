using HDF5, Printf, Plots, Plots.Measures

include("../src/fourierBasis.jl");
include("../src/computeRadii.jl");
include("../src/plotSave.jl");

#labels
rMinLab = "\$r_{min}\$";      #minimum inner boundary
rMaxLab = "\$r_{max}\$";      #maximum inner boundary
rOutLab = "\$\\Gamma^o\$ (radius \$R\$)";      #outer boundary
boundaryLab = "\$\\Gamma_b^i\$ (radius \$c(b_0+b)\$)";  #inner boundary

#run setup to get radius squash (uses default :/ )
include("../scenarios/svsector/setup.jl");

#file = "/projects/SIAllocation/stokes/svsector/svsector_029.h5";
file = "pubs/2022_radius_diagram.h5";

if isfile(file)
  println("Reading from: $file");
  f = h5open(file,"r");
  r    = read(f,"r");
  th   = read(f,"th");
  rMin = read(f,"rMin");
  rMax = read(f,"rMax");
  rOut = read(f,"rOut");
  close(f);
else
  #radii
  rMin = 0.5; 
  rMax = 1.5; 
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
  
  #draw randomly from prior
  samples = rand(mcmcP.prior);
  
  #compute radii
  dth = pi/180;
  th=0:dth:2*pi;
  r = computeRadii(samples,th);
end

# #get index
# idx = argmax(lpdfs[:,3]);
# r = r[idx,:];

#plot
sz = 500;
p = plot(proj=:polar,size=(sz,sz), margin=5Plots.mm);
plot!(p,th,[rMin],ls=:dash,lab=rMinLab);
plot!(p,th,r,lab=boundaryLab);
plot!(p,th,[rMax],ls=:dash,lab=rMaxLab);
plot!(p,th,[rOut],c=:black,lab=rOutLab);

plotSave(p,"pubs/2022_radius_diagram");

if !isfile(file)
  f = h5open(file,"w");
  write(f,"r",r);
  write(f,"th",collect(th));
  write(f,"rMin",rMin);
  write(f,"rMax",rMax);
  write(f,"rOut",rOut);
  close(f);
  println("Wrote data to: $file");
end
