function g = VanPolDiffusion(t,x,p)
    sigma = p(2);
    g = [0.0; sigma];
    %g  = repmat(g,length(x)/2,1);
end