%{
This is the driver for exercise 2. 
This file contains:
    - A test for correctness using the problem given in exercise 1.4.
      To test for this problem set the variable largeProb to false. The
      variable is found in the top of the script.
    - A test for efficiency under growing problem using a random EQP
      generator with bounds. To test for this problem set the 
      variable largeProb to true. The size of every increment is 
      controlled by n_large. The script will run 20 iterations
%}
%% Testing correctness
largeProb = false;
n_large = 50;


rng(12)
times = zeros(20, 3);
bs = [];
iterations = zeros(20, 3);
i = 1;
solution = zeros(20,1);
ldltimes = zeros(20,1);

%Plot settings
plotIterations = true;
plotTimes = true;
plotSolution = true;
plotldlvstime = false;
plotcvx = false;

%We test the QP solver for various values of b for the given problem.
for b1 = linspace(8.5, 18.68, 20)

    %Create a quadratic program
    if ~largeProb
    n = 5;
    H = [5.0 1.86 1.24 1.48 -0.46; ...
        1.86 3.0 0.44 1.12 0.52; ...
        1.24 0.44 3.8 1.56 -0.54; ...
        1.48 1.12 1.56 7.2 -1.12; ...
        -0.46 0.52 -0.54 -1.12 7.8];
    g = [-16.1; -8.5; -15.7; -10.02; -18.68];
    A = [16.1 8.5 15.7 10.02 18.68; 1.0 1.0 1.0 1.0 1.0]';
    b = [b1;1];


    l = zeros(5,1);
    u = ones(5,1);
    end

    if largeProb
        n = n_large*i;
        [H, g, A, b] = generateRandomEQP(n,n/2);
        l = zeros(n,1);
        u = ones(n,1);
    end

   % Solve by CVX
   if plotcvx
   start = cputime;

    cvx_begin quiet
        %cvx_precision low
        variable x(n)
        minimize( 1/2 * x' * H *x + g'*x)
        subject to
            A' * x == b
            l <= x
            x <= u
    cvx_end

    times(i,1) = cputime-start;
    iterations(i,1) = cvx_slvitr;
   end

    %solve by quadprog
   
    options = optimset('Display', 'off');
    options = optimset(options, 'Algorithm', 'interior-point-convex');
    start = cputime;
    [x2, optval, exitflag,output] = quadprog(H, g, [],[], A', b,l,u, 0, options);
    
    times(i,2) = cputime-start;
    iterations(i,2) = output.iterations;
    
    
    
    %solve by own solver
    x0 = zeros(n,1);
    s0 = ones(2*n,1);
    y0 = ones(length(b),1);
    z0 = ones(2*n,1);

    start = cputime;
    [x,y,z,s, iter, ldltime] = primalDualInteriorMethod_box(H,g,A,b,l,u,x0,y0,z0,s0);
    times(i,3) = cputime-start;
    iterations(i, 3) = iter;
    ldltimes(i,1) = ldltime;
    
    solution(i,1) = mean(sqrt((x-x2).^2));

    
    bs = [bs;b1];
    disp("Iteration " +i +"/" + 20);
    disp("Mean error: " + mean(sqrt((x-x2).^2)))
    disp("quadprog time: " + times(i,2))
    disp("our time: " + times(i,3))
    i = i+1;
    
    if exitflag ~= 1
        disp("Iteration with b(1)="+b1+" is infeasible!!")
    end
end


% PLOT ITERATIONS
if largeProb
    bs = n_large:n_large:(20*n_large);
end
if plotIterations
    figure;
    disp("Plotting # of iterations!")
    if plotcvx
        plot(bs, iterations(:,1))
    end
    hold on
    plot(bs, iterations(:,2))
    plot(bs, iterations(:,3))
    hold off
    if plotcvx
        legend(["cvx iterations", "quadprog interior-point iterations", "Our solver"])
    else
        legend(["quadprog interior-point iterations", "Our solver"])
    end
    
    ylabel("iterations")
    ylim([0 max(max(iterations))+2])
    if largeProb
         xlabel("n")
         title("iterations vs n")
    else
        xlabel("b(1)")
        title("iterations vs b(1)")
    end
    
end

if plotTimes
    disp("Plotting times!")
    figure;
    %Note: Time plotted is total wall time between starting and ending
    %solver
    % Not CPU time (since CPU time counts the number of cores used)
    if largeProb
        if plotcvx
            plot(bs, (times(:,1)))
        end
        hold on
        plot(bs, (times(:,2)))
        plot(bs, (times(:,3)))
        hold off
        if plotcvx
            legend(["cvx time", "quadprog interior-point time", "Our solver"])
        else
            legend(["quadprog interior-point time", "Our solver"])
        end

        ylabel("t [s]")
         xlabel("n")
         title("time vs n")
    else
        if plotcvx
            semilogy(bs, times(:,1))
        end
        hold on
        semilogy(bs, times(:,2))
        semilogy(bs, times(:,3))
        hold off
        if plotcvx
            legend(["cvx time", "quadprog interior-point time", "Our solver"])
        else
            legend(["quadprog interior-point time", "Our solver"])
        end

        ylabel("t [s]")
        xlabel("b(1)")
        title("time vs b(1)")
    end
end

if plotSolution
    figure;
    disp("Plotting solution!")
    plot(bs, solution)
    title("Correctness")
    ylabel("Mean Squared Error")
    legend("Error between quadprog and our solver")
    if largeProb
         xlabel("n")
    else
        xlabel("b(1)")
    end
end

if plotldlvstime
    disp('plotting ldl vs time')
    figure;
    hold off
    plot(bs, (times(:,2)))
    hold on
    plot(bs, (times(:,3)))
    plot(bs, (times(:,3))-(ldltimes(:,1)))
    hold off
    legend(["quadprog interior-point time", "Our solver", "Our solver without factorization"])
    ylabel("t [s]")
     xlabel("n")
     title("time vs n")
end

%% Large problem
%We look at the runtime of different algorithhms for a varying size problem
figure
largeProb = true;
n_large = 50;


rng(12)
times = zeros(20, 3);
bs = [];
iterations = zeros(20, 3);
i = 1;
solution = zeros(20,1);
ldltimes = zeros(20,1);

%Plot settings
plotIterations = true;
plotTimes = true;
plotSolution = true;
plotldlvstime = false;
plotcvx = false;

for b1 = linspace(8.5, 18.68, 20)

    %Create a quadratic program
    if ~largeProb
    n = 5;
    H = [5.0 1.86 1.24 1.48 -0.46; ...
        1.86 3.0 0.44 1.12 0.52; ...
        1.24 0.44 3.8 1.56 -0.54; ...
        1.48 1.12 1.56 7.2 -1.12; ...
        -0.46 0.52 -0.54 -1.12 7.8];
    g = [-16.1; -8.5; -15.7; -10.02; -18.68];
    A = [16.1 8.5 15.7 10.02 18.68; 1.0 1.0 1.0 1.0 1.0]';
    b = [b1;1];


    l = zeros(5,1);
    u = ones(5,1);
    end

    if largeProb
        n = n_large*i;
        [H, g, A, b] = generateRandomEQP(n,n/2);
        l = zeros(n,1);
        u = ones(n,1);
    end

   % Solve by CVX
   if plotcvx
   start = cputime;

    cvx_begin quiet
        %cvx_precision low
        variable x(n)
        minimize( 1/2 * x' * H *x + g'*x)
        subject to
            A' * x == b
            l <= x
            x <= u
    cvx_end

    times(i,1) = cputime-start;
    iterations(i,1) = cvx_slvitr;
   end

    %solve by quadprog
   
    options = optimset('Display', 'off');
    options = optimset(options, 'Algorithm', 'interior-point-convex');
    start = cputime;
    [x2, optval, exitflag,output] = quadprog(H, g, [],[], A', b,l,u, 0, options);
    
    times(i,2) = cputime-start;
    iterations(i,2) = output.iterations;
    
    
    
    %solve by own solver
    x0 = zeros(n,1);
    s0 = ones(2*n,1);
    y0 = ones(length(b),1);
    z0 = ones(2*n,1);

    start = cputime;
    [x,y,z,s, iter, ldltime] = primalDualInteriorMethod_box(H,g,A,b,l,u,x0,y0,z0,s0);
    times(i,3) = cputime-start;
    iterations(i, 3) = iter;
    ldltimes(i,1) = ldltime;
    
    solution(i,1) = mean(sqrt((x-x2).^2));

    
    bs = [bs;b1];
    disp("Iteration " +i +"/" + 20);
    disp("Mean error: " + mean(sqrt((x-x2).^2)))
    disp("quadprog time: " + times(i,2))
    disp("our time: " + times(i,3))
    i = i+1;
    
    if exitflag ~= 1
        disp("Iteration with b(1)="+b1+" is infeasible!!")
    end
end


% PLOT ITERATIONS
if largeProb
    bs = n_large:n_large:(20*n_large);
end
if plotIterations
    figure;
    disp("Plotting # of iterations!")
    if plotcvx
        plot(bs, iterations(:,1))
    end
    hold on
    plot(bs, iterations(:,2))
    plot(bs, iterations(:,3))
    hold off
    if plotcvx
        legend(["cvx iterations", "quadprog interior-point iterations", "Our solver"])
    else
        legend(["quadprog interior-point iterations", "Our solver"])
    end
    
    ylabel("iterations")
    ylim([0 max(max(iterations))+2])
    if largeProb
         xlabel("n")
         title("iterations vs n")
    else
        xlabel("b(1)")
        title("iterations vs b(1)")
    end
    
end

if plotTimes
    disp("Plotting times!")
    figure;
    %Note: Time plotted is total wall time between starting and ending
    %solver
    % Not CPU time (since CPU time counts the number of cores used)
    if largeProb
        if plotcvx
            plot(bs, (times(:,1)))
        end
        hold on
        plot(bs, (times(:,2)))
        plot(bs, (times(:,3)))
        hold off
        if plotcvx
            legend(["cvx time", "quadprog interior-point time", "Our solver"])
        else
            legend(["quadprog interior-point time", "Our solver"])
        end

        ylabel("t [s]")
         xlabel("n")
         title("time vs n")
    else
        if plotcvx
            semilogy(bs, times(:,1))
        end
        hold on
        semilogy(bs, times(:,2))
        semilogy(bs, times(:,3))
        hold off
        if plotcvx
            legend(["cvx time", "quadprog interior-point time", "Our solver"])
        else
            legend(["quadprog interior-point time", "Our solver"])
        end

        ylabel("t [s]")
        xlabel("b(1)")
        title("time vs b(1)")
    end
end

if plotSolution
    figure;
    disp("Plotting solution!")
    plot(bs, solution)
    title("Correctness")
    ylabel("Mean Squared Error")
    legend("Error between quadprog and our solver")
    if largeProb
         xlabel("n")
    else
        xlabel("b(1)")
    end
end

if plotldlvstime
    disp('plotting ldl vs time')
    figure;
    hold off
    plot(bs, (times(:,2)))
    hold on
    plot(bs, (times(:,3)))
    plot(bs, (times(:,3))-(ldltimes(:,1)))
    hold off
    legend(["quadprog interior-point time", "Our solver", "Our solver without factorization"])
    ylabel("t [s]")
     xlabel("n")
     title("time vs n")
end

