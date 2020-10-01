function sampleInnerGeometry(x=1)

  # for sane reproducibility in debugging, or parallel samples,
  # specify the random seed for this sample

  Random.seed!(x);

  a0 = 1.0;
  a  = [-1/2, -1/8, -1/8, 1/18, 1/18];
  b = [1/4, -1/8, 0.0, 1/4, -1/8];

  # ultimately we need to filter this to ensure the sample
  # does not lead to an illegal geometry (not too close to
  # zero--no radius or to two--the outer radius)

  return a0, a, b
end
