function [x,t,function_calls,hs,time] = ODEsolver(f,param,h,t0,T,x0, type,options)

if sum(strcmp(fieldnames(options), 'step_control')) == 1
    if islogical(options.step_control)
        
    else
        error("step_control must be a boolean")
    end
else
    if type == "Explicit Euler" || type == "RK4" || type == "DOPRI54" || type == "Implicit Euler"
        error("For the solvers Explicit Euler, Implicit Euler, RK4 and DOPRI54, step_control must be defined")   
    end
end

if sum(strcmp(fieldnames(options), 'control_type')) == 1
    if options.control_type == "I" || options.control_type == "PI" || options.control_type == "PID"
        
    else
        error("control_type must either I, PI or PID")
    end
else 
    options.control_type = "PI";
end

if sum(strcmp(fieldnames(options), 'initialStepSize')) == 1
    if islogical(options.initialStepSize)
        
    else
        error("initialStepSize must be a boolean")
    end
else
    if type == "Explicit Euler" || type == "RK4" || type == "DOPRI54"
        error("For the solvers Explicit Euler, RK4 and DOPRI54, initialStepSize must be defined")   
    end
end

if sum(strcmp(fieldnames(options), 'paths')) == 1
     if isnumeric(options.paths) && options.paths>0 && isaninteger(options.paths)
         
     else
         error("paths must be a positive integer")
     end
 else
     if type == "Explicit-Explicit" || type == "Implicit-Explicit"
         error("For the solvers Explicit-Explicit Implicit-Explicit, paths must be defined")   
     end
 end
 if sum(strcmp(fieldnames(options), 'g')) == 1
        % Cannot test if g is a function
 else
     if type == "Explicit-Explicit" || type == "Implicit-Explicit"
         error("For the solvers Explicit-Explicit or Implicit-Explicit, g must be defined")   
     end
 end
  if sum(strcmp(fieldnames(options), 'Jac')) == 1
        % Cannot test if Jac is a function
 else
     if type == "ESDIRK" || type == "Implicit-Explicit" || type == "Implicit Euler"
         error("For the solvers ESDIRK, Implicit-Explicit or Implicit Euler, Jac must be defined")   
     end
  end
  if sum(strcmp(fieldnames(options), 'Atol')) == 1
      if isnumeric(options.Atol) && 0<options.Atol
          Atol = options.Atol;
      else
          error("Atol must be a positive number")
      end
  else
     Atol = 0.00001;
  end
   if sum(strcmp(fieldnames(options), 'Rtol')) == 1
      if isnumeric(options.Rtol) && 0<options.Rtol
          Rtol = options.Rtol;
      else
          error("Rtol must be a positive number")
      end
  else
     Rtol = 0.00001;
   end
   if sum(strcmp(fieldnames(options), 'eps_tol')) == 1
      if isnumeric(options.eps_tol) && 0<options.eps_tol
          eps_tol = options.eps_tol;
      else
          error("eps_tol must be a positive number")
      end
  else
     eps_tol = 0.7;
   end
  if sum(strcmp(fieldnames(options), 'hmin')) == 1
      if isnumeric(options.hmin) && 0<options.hmin
          hmin = options.hmin;
      else
          error("hmin must be a positive number")
      end
  else
     hmin = 0.00001;
  end
  if sum(strcmp(fieldnames(options), 'hmax')) == 1
      if isnumeric(options.hmax) && 0<options.hmax
          hmax = options.hmax;
      else
          error("hmax must be a positive number")
      end
  else
     hmax = 5;
  end
  
    if sum(strcmp(fieldnames(options), 'newtonTolerance')) == 1
      if isnumeric(options.newtonTolerance) && 0<options.newtonTolerance
          newtonTolerance = options.newtonTolerance;
      else
          error("newtonTolerance must be a positive number")
      end
  else
     newtonTolerance = 1.0e-8;
   end
  
  if sum(strcmp(fieldnames(options), 'newtonMaxiterations')) == 1
      if isnumeric(options.newtonMaxiterations) && 0<options.newtonMaxiterations
          newtonMaxiterations = options.newtonMaxiterations;
      else
          error("newtonMaxiterations must be a positive number")
      end
  else
     newtonMaxiterations = 100;
  end
  if sum(strcmp(fieldnames(options), 'seed')) == 1
      if isaninteger(options.seed) && 0<options.seed
          seed = options.seed;
      else
          error("seed must be a positive integer")
      end
  else
        seed = randi([0 1000000],1);
  end
  
time = nan;

if type == "Explicit Euler"
    A = [0];
    b = [1]';
    c = [0]';
    p = 1;
    
    if options.step_control == 1 
          if options.initialStepSize == 1
              if options.control_type == "I"
                  tStart = cputime;
                  [x,t,function_calls,hs,rs] = explicitRungeKuttaDoubling(f,param,h,t0,T,x0,A,b,c,Atol,Rtol,hmin,hmax,eps_tol,p,true); % euler with everything
time = cputime - tStart;
              elseif options.control_type == "PI"
                  tStart = cputime;
                  [x,t,function_calls,hs,rs] = explicitRungeKuttaDoublingPI(f,param,h,t0,T,x0,A,b,c,Atol,Rtol,hmin,hmax,eps_tol,p,true); % euler with everything
time = cputime - tStart;
              else
                  tStart = cputime;
                  [x,t,function_calls,hs,rs] = explicitRungeKuttaDoublingPID(f,param,h,t0,T,x0,A,b,c,Atol,Rtol,hmin,hmax,eps_tol,p,true); % euler with everything
time = cputime - tStart;
              end
          else
              if options.control_type == "I"
                  tStart = cputime;
                  [x,t,function_calls,hs,rs] = explicitRungeKuttaDoubling(f,param,h,t0,T,x0,A,b,c,Atol,Rtol,hmin,hmax,eps_tol,p,false); % euler with everything
time = cputime - tStart;
              elseif options.control_type == "PI"
                  tStart = cputime;
                  [x,t,function_calls,hs,rs] = explicitRungeKuttaDoublingPI(f,param,h,t0,T,x0,A,b,c,Atol,Rtol,hmin,hmax,eps_tol,p,false); % euler with everything
time = cputime - tStart;
              else
                  tStart = cputime;
                  [x,t,function_calls,hs,rs] = explicitRungeKuttaDoublingPID(f,param,h,t0,T,x0,A,b,c,Atol,Rtol,hmin,hmax,eps_tol,p,false); % euler with everything
time = cputime - tStart;
              end
          end
    else
         tStart = cputime;
        [x,t,function_calls,hs] = explicitRungeKutta(f,param,h,t0,T,x0,A,b,c); % Standard euler
time = cputime - tStart;
    end

    
elseif type == "Implicit Euler"
    if options.step_control == 1 
          if options.initialStepSize == 1
                tStart = cputime;
                [x,t,function_calls,hs,rs] = implicitEulerDoubling(f,options.Jac,param,h,t0,T,x0,Atol,Rtol,hmin,hmax,eps_tol,options.initialStepSize,newtonTolerance, newtonMaxiterations); % euler with everything
time = cputime - tStart;
          else
                tStart = cputime;
                [x,t,function_calls,hs,rs] = implicitEulerDoubling(f,options.Jac,param,h,t0,T,x0,Atol,Rtol,hmin,hmax,eps_tol,options.initialStepSize,newtonTolerance, newtonMaxiterations); % euler with step control only
time = cputime - tStart;
          end
    else
        tStart = cputime;
        [x,t,function_calls,hs] = implicitEulerFixed(f,options.Jac,param,h,t0,T,x0,newtonTolerance,newtonMaxiterations); % Standard euler
time = cputime - tStart;
    end
    
    
elseif type == "RK4"
    A = [0 0 0 0; 0.5 0 0 0; 0 0.5 0 0; 0 0 1 0];
    b = [1/6 2/6 2/6 1/6]';
    c = [0 1/2 1/2 1]';
    p = 4;


    if options.step_control == 1 
          if options.initialStepSize == 1
              if options.control_type == "I"
                   tStart = cputime;
                  [x,t,function_calls,hs,rs] = explicitRungeKuttaDoubling(f,param,h,t0,T,x0,A,b,c,Atol,Rtol,hmin,hmax,eps_tol,p,true); % euler with everything
time = cputime - tStart;
              elseif options.control_type == "PI"
                   tStart = cputime;
                  [x,t,function_calls,hs,rs] = explicitRungeKuttaDoublingPI(f,param,h,t0,T,x0,A,b,c,Atol,Rtol,hmin,hmax,eps_tol,p,true); % euler with everything
time = cputime - tStart;
              else
                   tStart = cputime;
                  [x,t,function_calls,hs,rs] = explicitRungeKuttaDoublingPID(f,param,h,t0,T,x0,A,b,c,Atol,Rtol,hmin,hmax,eps_tol,p,true); % euler with everything
time = cputime - tStart;
              end
          else
              if options.control_type == "I"
                   tStart = cputime;
                  [x,t,function_calls,hs,rs] = explicitRungeKuttaDoubling(f,param,h,t0,T,x0,A,b,c,Atol,Rtol,hmin,hmax,eps_tol,p,false); % euler with everything
time = cputime - tStart;
              elseif options.control_type == "PI"
                   tStart = cputime;
                  [x,t,function_calls,hs,rs] = explicitRungeKuttaDoublingPI(f,param,h,t0,T,x0,A,b,c,Atol,Rtol,hmin,hmax,eps_tol,p,false); % euler with everything
time = cputime - tStart;
              else
                   tStart = cputime;
                  [x,t,function_calls,hs,rs] = explicitRungeKuttaDoublingPID(f,param,h,t0,T,x0,A,b,c,Atol,Rtol,hmin,hmax,eps_tol,p,false); % euler with everything
time = cputime - tStart;
              end
          end
    else
         tStart = cputime;
        [x,t,function_calls,hs] = explicitRungeKutta(f,param,h,t0,T,x0,A,b,c); % Standard euler
time = cputime - tStart;
    end
    
elseif type == "DOPRI54"
    c =[0 1/5 3/10 4/5 8/9 1 1]';
    A = [0 0 0 0 0 0 0 ;
     1/5 0 0 0 0 0 0
     3/40 9/40 0 0 0 0 0 ;
     44/45 -56/15 32/9 0 0 0 0 ;
     19372/6561 -25360/2187 64448/6561 -212/729 0 0 0;
     9017/3168 -355/33 46732/5247 49/176 -5103/18656 0 0;
     35/384 0 500/1113 125/192 -2187/6784 11/84 0];
    b = [35/384 0 500/1113 125/192 -2187/6784 11/84 0]';
    bhat = [5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40]';
    p = 5;
    

    if options.step_control == 1 
          if options.initialStepSize == 1
              if options.control_type == "I"
                   tStart = cputime;
                  [x,t,function_calls,hs,rs] = explicitRungeKuttaEmbedded(f,param,h,t0,T,x0,A,b,bhat,c,Atol,Rtol,hmin,hmax,eps_tol,p,true); % DOPRI54 with everything
time = cputime - tStart;
              elseif options.control_type == "PI"
                   tStart = cputime;
                  [x,t,function_calls,hs,rs] = explicitRungeKuttaEmbeddedPI(f,param,h,t0,T,x0,A,b,bhat,c,Atol,Rtol,hmin,hmax,eps_tol,p,true); % DOPRI54 with everything
time = cputime - tStart;
              else
                   tStart = cputime;
                  [x,t,function_calls,hs,rs] = explicitRungeKuttaEmbeddedPID(f,param,h,t0,T,x0,A,b,bhat,c,Atol,Rtol,hmin,hmax,eps_tol,p,true); % DOPRI54 with everything
time = cputime - tStart;
              end
                
          else
              if options.control_type == "I"
                   tStart = cputime;
                  [x,t,function_calls,hs,rs] = explicitRungeKuttaEmbedded(f,param,h,t0,T,x0,A,b,bhat,c,Atol,Rtol,hmin,hmax,eps_tol,p,false); % DOPRI54 with everything
time = cputime - tStart;
              elseif options.control_type == "PI"
                   tStart = cputime;
                  [x,t,function_calls,hs,rs] = explicitRungeKuttaEmbeddedPI(f,param,h,t0,T,x0,A,b,bhat,c,Atol,Rtol,hmin,hmax,eps_tol,p,false); % DOPRI54 with everything
time = cputime - tStart;
              else
                   tStart = cputime;
                  [x,t,function_calls,hs,rs] = explicitRungeKuttaEmbeddedPID(f,param,h,t0,T,x0,A,b,bhat,c,Atol,Rtol,hmin,hmax,eps_tol,p,false); % DOPRI54 with everything
time = cputime - tStart;
              end
          end
    else
       tStart = cputime;
        [x,t,function_calls,hs] = explicitRungeKutta(f,param,h,t0,T,x0,A,b,c); % Standard euler
time = cputime - tStart;
    end
    
elseif type == "ESDIRK23"
    % MANGLER MANGLER MANGLER MANGLER MANGLER MANGLER MANGLER MANGLER 
                   tStart = cputime;
    [t,x, ~ ,function_calls,hs] = ESDIRK(f,options.Jac,t0,T,x0,h,Atol,Rtol,param);
time = cputime - tStart;
elseif type == "Explicit-Explicit"
                   tStart = cputime;
    [x, t,function_calls,hs] = SDEExplicitExplicit(x0, f, options.g, h, t0, T, param, options.paths,seed);
time = cputime - tStart;
    
elseif type == "Implicit-Explicit"
    if options.paths>1
        error("Multiple paths cannot be evaluated on the same time with the Implicit-Explicit SDE solver. Use Explicit-Explicit or setup an extern loop.")
    end
                   tStart = cputime;
    [x, t,function_calls,hs] = SDEImplicitExplicit(x0, f,options.Jac, options.g, h, t0, T, param, options.paths, newtonTolerance, newtonMaxiterations,seed);
time = cputime - tStart;
else
    error("The ODE solver " + type + "does not exist; possible solvers are Explicit Euler, Implicit Euler, RK4, DOPRI54, ESDIRK23, Explicit-Explicit and Implicit-Explicit")
end

end

function [bol_val] = isaninteger(x)
    bol_val = isfinite(x) & x==floor(x);
end