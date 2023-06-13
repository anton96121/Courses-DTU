function [W,Tw,h] = brownian_motion(h,N,dim, seed)

if isnumeric(seed)
  rng(seed);
end

W = sqrt(h)*randn(dim,N);
W = cumsum(W,2);
W = [zeros(dim,1) W];
Tw = 0:h:N;
end

