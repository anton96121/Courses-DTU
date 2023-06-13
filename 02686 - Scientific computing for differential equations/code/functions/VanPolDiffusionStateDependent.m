function g = VanPolDiffusionStateDependent(t,x,p)
    sigma = p(2);
    g = [0.0; sigma*(1+x(1)^2)];
    %g  = repmat(g,length(x)/2,1);
end