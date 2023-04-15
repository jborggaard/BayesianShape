using HDF5, Printf, Plots, Plots.Measures

include("../src/fourierBasis.jl");
include("../src/computeRadii.jl");
include("../src/plotSave.jl");

scen = "vortsensor";

#run setup to get radius squash (uses default :/ )
include("../scenarios/$(scen)/setup.jl");

dir   = "/projects/SIAllocation/stokes/$(scen)";
files = [ "$(scen)_0$(i).h5" for i=17:22 ];

#angles and fourier basis with which we'll take the inner product
dth = pi/180;
th=0:dth:2*pi;
fftb = exp.(im*th[1:end-1]);

#lists
nsamp = h5read(dir*"/"*files[1],"sampComplete"); #assume the same for all files
areas = zeros(nsamp,0);
phis  = zeros(nsamp,0);

#pphi  = plot(layout=length(files),xlab="Samples",ylab="Angle");
#parea = plot(layout=length(files),xlab="Samples",ylab="Area");
pphi  = plot(xlab="Samples",ylab="Mean Angle (Radians)");
parea = plot(xlab="Samples",ylab="Mean Area");
for i=1:length(files)
  file = files[i];
  #file   = "../svglobal/svglobal_040.h5";
  f = h5open(dir*"/"*file);
  samples = read(f,"samples");
  close(f);

  #compute radii
  r = computeRadii(samples,th);

  #area
  area = dth*sum(0.5*r[:,1:end-1].^2,dims=2)[:];
  
  #fft
  fftr = r[:,1:end-1]*fftb;
  #phi = atan.(imag.(fftr)./real.(fftr));
  phi = atan.(imag.(fftr),real.(fftr));

  # #test plot
  # sz = 500;
  # i = size(samples,1);
  # p = plot(proj=:polar,size=(sz,sz));#;kwargs...);
  # plot!(p, th, r[i,:]);
  # println("$(i): real=$(real(fftr[i])), imag=$(imag(fftr[i])), phi=$(phi[i]), phi(deg)=$(180*phi[i]/pi), area=$(area[i])");
  # plot!(p, [0,phi[i]],[0,maximum(r[i,:])]);
  # savefig("pubs/angle_test_$(file)_$(i).png");

  #running average
  mphi  = cumsum(phi)  ./ (1:length(phi));
  marea = cumsum(area) ./ (1:length(area));

  #running stdev
  sphi  = sqrt.( cumsum((phi.-mean(phi)).^2)   ./ (1:length(phi)) );
  sarea = sqrt.( cumsum((area.-mean(area)).^2) ./ (1:length(area)) );

  #plot attributes
  label = "Chain $i";

  #plot
  plot!(pphi ,1:length(mphi) ,mphi ,lab=label);
  plot!(parea,1:length(marea),marea,lab=label);
  #plot!(pphi,1:length(mphi),mphi,c=i,lab=label*" (mean)");
  #plot!(pphi,1:length(mphi),sphi,c=i,ls=:dash,lab=label*" (std)");
  #plot!(parea,1:length(marea),marea,c=i,lab=label*" (mean)");
  #plot!(parea,1:length(marea),marea+sarea,c=i,ls=:dash,lab=label*" (std)");
  #plot!(parea,1:length(marea),marea-sarea,c=i,ls=:dash,lab=false);#,lab=label*" (std)");

  global areas = hcat(areas,area);
  global phis  = hcat(phis, phi);
end

#plot!(pphi ,leg=:outerright,size=(1000,400),margin=5.0mm,ylims=(-1.00, 1.00));
#plot!(parea,leg=:outerright,size=(1000,400),margin=5.0mm,ylims=( 3.65, 3.85));
#plot!(pphi ,leg=:bottomright,margin=5.0mm,ylims=(-1.00, 1.00));
#plot!(parea,leg=:bottomright,margin=5.0mm,ylims=( 3.65, 3.85));
plot!(pphi ,leg=:bottomright);#,ylims=(-1.00, 1.00));
plot!(parea,leg=:bottomright,ylims=(0,5));#,ylims=( 3.65, 3.85));

#savefig(pphi,"pubs/2022_multiproposal_angle.png");
#println("Wrote: pubs/2022_multiproposal_angle.png");
#savefig(parea,"pubs/2022_multiproposal_area.png");
#println("Wrote: pubs/2022_multiproposal_area.png");
plotSave(pphi,"pubs/2022_stokes_compare_angle");
plotSave(parea,"pubs/2022_stokes_compare_area");



#see Gelman BDA3, pg284
function statGelmanRubin(psi)
  n,m = size(psi);
  psibar  = mean(psi);
  psibarj = mean(psi,dims=1)[:];
  
  B = n/(m-1)*norm(psibar.-psibarj)^2;

  s=zeros(m); 
  for j=1:m;
    for i=1:n; 
      s[j] += 1/(n-1)*(psi[i,j].-psibarj[j])^2; 
    end
  end
  W=mean(s);

  varp = (W*(n-1)+B)/n;
  Rhat = sqrt(varp/W);

  return Rhat;
end

#print some stats
grArea = statGelmanRubin(areas);
grPhi  = statGelmanRubin(phis);
println("Gelman-Rubin for area is : $grArea");
println("Gelman-Rubin for angle is: $grPhi");

for k=1:8
  comp = zeros(nsamp,length(files));
  for j=1:length(files)
    comp[:,j] = h5read(dir*"/"*files[j],"samples")[:,k];
  end
  grComp = statGelmanRubin(comp);
  println("Gelman-Rubin for component $k is : $grComp");
  if grComp > 1.1
    println("means for component $k are:");
    display(mean(comp,dims=1))
  end
end
