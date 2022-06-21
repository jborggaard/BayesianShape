using HDF5, Printf, Plots

include("../src/fourierBasis.jl");
include("../src/computeRadii.jl");
include("../src/plotSave.jl");

scen = "svsector";

#run setup to get radius squash (uses default :/ )
#include("../scenarios/$(scen)/setup.jl");

dir   = "/projects/SIAllocation/stokes/$(scen)";
files = [ "$(scen)_0$(i).h5" for i=18:29 ];

#first, build a list of beta values, which we'll use to color the figure
#betas = [];
#for file in files
#  append!(betas,[h5read(dir*"/"*file,"final_mcmc/beta")]);
#end
betas = [ h5read(dir*"/"*file,"final_mcmc/beta") for file in files ];
unique!(betas);


#angles and fourier basis with which we'll take the inner product
dth = pi/180;
th=0:dth:2*pi;
fftb = exp.(im*th[1:end-1]);

#pphi  = plot(layout=length(files),xlab="Samples",ylab="Angle");
#parea = plot(layout=length(files),xlab="Samples",ylab="Area");
pphi  = plot(xlab="Samples",ylab="Mean Angle (Radians)");
parea = plot(xlab="Samples",ylab="Mean Area");
for i=1:length(files)
  file = files[i];
  #file   = "../svglobal/svglobal_040.h5";
  f = h5open(dir*"/"*file);
  mcmc = read(f,"mcmc");
  beta = read(f,"final_mcmc/beta");
  samples = read(f,"samples");
  close(f);

  rho = sqrt(1-beta^2);

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
  #set label and line style based on mcmc type
  mcmcSp = split(mcmc,"|");
  if mcmcSp[1] == "pcn"
    label = "Vanilla pCN";
    ls = :dash
  else
    label = "$(mcmcSp[3]) Proposals";
    ls = :solid
  end
  label *= @sprintf(", \$\\rho\$=%5.3f",rho);
  #hack to only include unique labels - set to "" if it's already in the plot somewhere
  labelsSoFar = [ p.plotattributes[:label] for p in pphi.series_list ];
  label = (label in labelsSoFar) ? "" : label;
  #set color based on beta
  col = findall(==(beta), betas)[1];

  #plot
  plot!(pphi ,1:length(mphi) ,mphi ,c=col,ls=ls,lab=label);
  plot!(parea,1:length(marea),marea,c=col,ls=ls,lab=label);
  #plot!(pphi,1:length(mphi),mphi,c=i,lab=label*" (mean)");
  #plot!(pphi,1:length(mphi),sphi,c=i,ls=:dash,lab=label*" (std)");
  #plot!(parea,1:length(marea),marea,c=i,lab=label*" (mean)");
  #plot!(parea,1:length(marea),sarea,c=i,ls=:dash,lab=label*" (std)");
end

plot!(pphi ,leg=:outerright,size=(800,400),ylims=(-1.00, 1.00));
plot!(parea,leg=:outerright,size=(800,400),ylims=( 3.65, 3.85));

#savefig(pphi,"pubs/2022_multiproposal_angle.png");
#println("Wrote: pubs/2022_multiproposal_angle.png");
#savefig(parea,"pubs/2022_multiproposal_area.png");
#println("Wrote: pubs/2022_multiproposal_area.png");
plotSave(pphi,"pubs/2022_multiproposal_angle");
plotSave(parea,"pubs/2022_multiproposal_area");
