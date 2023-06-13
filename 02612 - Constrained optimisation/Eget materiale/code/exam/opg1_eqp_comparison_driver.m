%{
This is the driver for exercise 1. 
This file contains:
    - A test for correctness using the problem given in exercise 1.4.
    - A test for efficiency under growing problem using the recycling 
      problem given in week 5.
    - A benchmark test for matrix factorizations.
    - We vary the number of constraints instead of the number of variables
%}
%% Given problem
names = ["LDLdense", "LDLsparse", "LUdense", "LUsparse", "rangespace", "nullspace", "quadprog"];
tests = 20;
times = zeros(7,tests);
answers = zeros(5,tests,7);
bs = [];
i=1;
for b1 = linspace(8.5, 18.68, tests)
    bs = [bs;b1];
    n = 5;
    H = [5.0 1.86 1.24 1.48 -0.46; ...
        1.86 3.0 0.44 1.12 0.52; ...
        1.24 0.44 3.8 1.56 -0.54; ...
        1.48 1.12 1.56 7.2 -1.12; ...
        -0.46 0.52 -0.54 -1.12 7.8];
    g = [-16.1; -8.5; -15.7; -10.02; -18.68];
    A = [16.1 8.5 15.7 10.02 18.68; 1.0 1.0 1.0 1.0 1.0]';
    b = [b1;1];

    %We test all solvers one by one and save the answers
    start = cputime;
    [x, lambda] = EqualityQPSolver(H,g,A,b, "LDLdense");
    answers(:, i,1) = x;
    times(1, i) = cputime-start;
    start = cputime;
    [x, lambda] = EqualityQPSolver(H,g,A,b, "LDLsparse");
    answers(:, i,2) = x;
    times(2, i) = cputime-start;
    start = cputime;
    [x, lambda] = EqualityQPSolver(H,g,A,b, "LUdense");
    answers(:, i,3) = x;
    times(3, i) = cputime-start;
    start = cputime;
    [x, lambda] = EqualityQPSolver(H,g,A,b, "LUsparse");
    answers(:, i,4) = x;
    times(4, i) = cputime-start;
    start = cputime;
    [x, lambda] = EqualityQPSolver(H,g,A,b, "rangespace");
    answers(:, i,5) = x;
    times(5, i) = cputime-start;
    start = cputime;
    [x, lambda] = EqualityQPSolver(H,g,A,b, "nullspace");
    answers(:, i,6) = x;
    times(6, i) = cputime-start;
    options = optimset('Display', 'off');
    start = cputime;
    [x, lambda] = quadprog(H,g,[],[],A',b,[],[],[],options);
    answers(:, i,7) = x;
    times(7, i) = cputime-start;


    i = i+1;
end


%We compare answers to quadprog
answer_diff = zeros(6,20);
for i=1:6
    answer_diff(i,:) = mean(sqrt((answers(:,:,i)-answers(:,:,7)).^2));
end



hold off
for i=1:size(times,1)
    plot(bs, times(i,:))
    hold on
end
legend(names)
xlabel("b(1)")
ylabel("time [log s]")

%and plot the answers
figure
for i=1:6
    plot(bs, (answer_diff(i,:)))
    hold on
end
legend(names(1:6))
xlabel("b(1)")
ylabel("Error relative to quadprog")
    



%% Recycling problem
names = ["LDLdense", "LDLsparse", "LUdense", "LUsparse", "rangespace", "nullspace", "quadprog"];
tests = 20;
times = ones(7,tests)*tests;
ns = zeros(tests,1);
%We test every solver using the Recycling problem of different sizes n.
for i = 1:tests
    disp(i)
    n = i*(2000/tests);
    ns(i) = n;

    [H, g, A, b] = ProblemEQPRecycling(n, 0.2, 1);
    for k =1:1
        j=1;

        start = cputime;
        [x, lambda] = EqualityQPSolver(H,g,A,b, "LDLdense");
        times(j, i) = cputime-start;
        j = j+1;

        start = cputime;
        [x, lambda] = EqualityQPSolver(H,g,A,b, "LDLsparse");
        times(j, i) = cputime-start;
        j = j+1;
        %
        start = cputime;
        [x, lambda] = EqualityQPSolver(H,g,A,b, "LUdense");
        times(j, i) = cputime-start;
        j = j+1;
        %
        start = cputime;
        [x, lambda] = EqualityQPSolver(H,g,A,b, "LUsparse");
        times(j, i) = min(cputime-start, times(j,i));
        j = j+1;

        %
        start = cputime;
        [x, lambda] = EqualityQPSolver(H,g,A,b, "rangespace");
        times(j, i) = cputime-start;
        j = j+1;

        start = cputime;
        [x, lambda] = EqualityQPSolver(H,g,A,b, "nullspace");
        times(j, i) = cputime-start;
        j = j+1;

        options = optimset('Display', 'off');
        start = cputime;
        [x, lambda] = quadprog(H,g,[],[],A',b,[],[],[],options);
        times(j, i) = cputime-start;
        j = j+1;
        %
    end

    disp("n="+n)
end

%We plot the time for each method.
figure
hold off
for i=1:size(times,1)
    plot(ns, times(i,:))
    hold on
end
legend(names, 'Location','northwest')
xlabel("n")
ylabel("time [s]")
title("Growing size problem")


%% Factorization benchmark
tests = 10;
times = ones(5,tests)*100;
ns = zeros(tests,1);
%We bench LU vs LDL versus Cholesky of varying sizes
for i = 1:tests
    n = i*(5000/tests);
    ns(i) = n;

    [H,g,A,b,x,lambda] = generateRandomEQP(n,n);
    KKT = [H -A;-A', zeros(size(A,2), size(A,2))];
    for k =1:1
        j=1;
        
        start = cputime;
        x = lu(KKT,'vector');
        times(j, i) = min(cputime-start, times(j,i));
        j = j+1;

        start = cputime;
        x = lu(H,'vector');
        times(j, i) = min(cputime-start, times(j,i));
        j = j+1;

        start = cputime;
        x = ldl(KKT,'vector');
        times(j, i) = min(cputime-start, times(j,i));
        j = j+1;

        start = cputime;
        x = ldl(H,'vector');
        times(j, i) = min(cputime-start, times(j,i));
        j = j+1;

        start = cputime;
        x = chol(H);
        times(j, i) = min(cputime-start, times(j,i));
        j = j+1;


    end

    disp("n="+n)
end


%We then plot the time for each method and target.
hold off
for i=1:5%size(times,1)
    semilogy(ns, times(i,:))
    hold on
end
legend(["LU, KKT", "LU, H", "LDL, KKT", "LDL, H", "Cholesky, H"])
xlabel("n")
ylabel("time [s]")
title("Benchmark of factorizations")

%% m dependent
names = ["rangespace", "nullspace"];
tests = 10;
times = ones(2,tests)*100;
ms = zeros(tests,1);
top = 3000;
%To compare range space and null space, we solve random EQPs with varying
%number of constraints.
for i = 1:tests
    m = i*(top/tests);
    ms(i) = m;

    [H, g, A, b] = generateRandomEQP(top, m);
    for k =1:1
        j=1;


        start = cputime;
        [x, lambda, facTime_R] = EqualityQPSolver(H,g,A,b, "rangespace");
        times(j, i) = min(cputime-start, times(j,i));
        j = j+1;

        start = cputime;
        [x, lambda, facTime_N] = EqualityQPSolver(H,g,A,b, "nullspace");
        times(j, i) = min(cputime-start, times(j,i));
        j = j+1;

    end

    disp("m="+m)
end

figure
hold off
for i=1:2
    plot(ms, times(i,:))
    hold on
end
plot(ms, times(1,:)-facTime_R)
plot(ms, times(2,:)-facTime_N)
xline(1950)
legend([names,'range space -Cholesky','null space -QR','Theoretical tipping point'], 'Location','northwest')
xlabel("m")
ylabel("time [s]")
title("Growing constraints problem, n=3000")

