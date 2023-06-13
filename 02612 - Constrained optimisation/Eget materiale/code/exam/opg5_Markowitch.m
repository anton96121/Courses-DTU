close all

%%
%{
This is the driver for exercise 5. 
This file contains:
    - Optimal portfolio and its risk for exercise 5.3
    - Computation of the efficient frontier and optimal portfolio as a
    function of the return
    - Solution of efficient frontier of the bi-criterion problem using both
    EQP and QP solvers
    - Computation of the efficient frontier and optimal portfolio as a
    function of return when we add a risk-free asset
    - Finding and plotting optimal point for R=14 with a risk-free asset
%}

%% Exercise 5.1-5.3
%Find the solution when going for R=12, and the accompanying risk
returns = [16.1, 8.5, 15.7, 10.02, 18.68]; 
R = 12; 
covariance = [2.5 .93 .62 .74 -.23;
              0.93 1.5 0.22 0.56 0.26;
              .62  .22 1.9  .78  -0.27;
              .74 .56 .78 3.6 -0.56;
              -0.23 0.26 -0.27 -0.56 3.9];  

H = covariance;
f = [];

A1 = returns;
b1 = R;

A2 = [1,1,1,1,1];
b2 = 1;

Aeq = [A1; A2];
beq = [b1; b2];

Aineq = -eye(5);
bineq = zeros(5,1); 

x = quadprog( H, f, Aineq, bineq, Aeq, beq)

port_risk = x'*covariance*x

%% Exercise 5.4, Computing the efficient frontier

Rs = min(returns):0.01:max(returns);
port_risks = zeros(length(Rs),1);
options = optimset('Display', 'off');
optimalPorts = zeros(length(Rs),5);

%Compute the efficient frontier for each possible return
for i = 1:length(Rs)
beq = [Rs(i); b2];
x = quadprog( H, f, Aineq, bineq, Aeq, beq,[],[],[],options); 
port_risks(i) = x'*covariance*x;
optimalPorts(i,:) = x;

end

%% Exercise 5.4, Plotting the efficient frontier
figure;
hold on

opt_return1 = Rs(port_risks==min(port_risks));
opt_risk1 = min(port_risks);

effRs = Rs(find(Rs == opt_return1):end);
effRisk = port_risks(find(Rs == opt_return1):end);

%Plot efficient frontier for each possible return
plot(Rs,port_risks,'b', returns, diag(covariance),'ro')
h2 = plot(opt_return1,opt_risk1 ,'ko','MarkerFaceColor', 'g');
h1 = plot(effRs,effRisk , 'r');
title('The Efficient Frontier')
xlabel('Return [%]')
ylabel('Risk [Var]')
legend([h1,h2],{'The efficient frontier','Minimum efficient return'}, 'Location', 'northwest')
hold off

%% Exercise 5.4, Plotting the optimal port folio as a function of return
figure
h1 = plot(Rs,optimalPorts(:,1));
hold on
h2 = plot(Rs,optimalPorts(:,2));
h3 = plot(Rs,optimalPorts(:,3));
h4 = plot(Rs,optimalPorts(:,4));
h5 = plot(Rs,optimalPorts(:,5));
legend([h1,h2,h3,h4,h5], {'Security 1','Security 2','Security 3','Security 4','Security 5'})
xlabel('Return [%]')
ylabel('Percentage of the portfolio [%]')
title('Portfolio composition')

%% Exercise 5.4, Minimum risk
% Find the portfolio with the smallest possible risk.
beq = [opt_return1; b2];
x_min = quadprog( H, f, Aineq, bineq, Aeq, beq,[],[],[],options)
port_risks_min = x'*covariance*x

%% Exercise 5.5-5.7, Bi-criterion Optimization
%Setup for bi-criterion
H = covariance;
f = -returns;

Aeq = [1,1,1,1,1];
beq = 1;

Aineq = -eye(5);
bineq = zeros(5,1); 


%% Exercise 5.5-5.7, Solving for different alpha with and without short-selling
trials = 10;
x_short = zeros(trials,5,3);
x_nonshort = zeros(trials,5,3);

alphas = 1:1:trials;
alphas = alphas./(trials);
port_risk_short = zeros(trials,1,3);
port_risk_nonshort = zeros(trials,1,3);
port_return_short = zeros(trials,1,3);
port_return_nonshort = zeros(trials,1,3);

x0 = zeros(5,1);
s0 = ones(2*5,1);
y0 = ones(length(beq),1);
z0 = ones(2*5,1);
n = 5;
if trials<11
    cvx_solve = true;
else
    cvx_solve = false;
end

%optimal port. weights found by quadprog
for i = 1:trials
    if mod(i,trials/10)==0
        disp(i)
    end
    % Without shorting, i.e. an IQP problem
    x_nonshort(i,:,1) = quadprog( alphas(i).*H, (1-alphas(i)).*f', Aineq, bineq, Aeq, beq,[],[],[],options); 
    x_nonshort(i,:,2) = primalDualInteriorMethod_box(alphas(i).*H, (1-alphas(i)).*f',Aeq',beq,zeros(5,1),ones(5,1),x0,y0,z0,s0);
    if cvx_solve
        cvx_begin quiet
            %cvx_precision low
            variable x(n)
            minimize( 1/2 * x' *  (alphas(i).*H) *x + ((1-alphas(i)).*f)*x)
            subject to
                Aeq * x == beq
                zeros(5,1) <= x
                x <= ones(5,1)
        cvx_end
        x_nonshort(i,:,3) = x; 
    end
    
    % With shorting, i.e. an EQP problem
    x_short(i,:,1) = quadprog( alphas(i).*H, (1-alphas(i)).*f', [], [], Aeq, beq,[],[],[],options); 
    x_short(i,:,2) = EqualityQPSolver(alphas(i).*H, (1-alphas(i)).*f',Aeq', beq, "rangespace");
    if cvx_solve
        cvx_begin quiet
            %cvx_precision low
            variable x(n)
            minimize( 1/2 * x' *  (alphas(i).*H) *x + ((1-alphas(i)).*f)*x)
            subject to
                Aeq * x == beq
        cvx_end
        x_short(i,:,3) = x; 
    end
    
    port_risk_short(i,2) = x_short(i,:,2)*covariance*x_short(i,:,2)';
    port_risk_nonshort(i,2) = x_nonshort(i,:,2)*covariance*x_nonshort(i,:,2)';
    
    port_return_short(i,2) = -f*x_short(i,:,2)';
    port_return_nonshort(i,2) = -f*x_nonshort(i,:,2)';
end


%Plot the MSE btwn our method, quadprog, and cvx, solving the EQP problem
figure
h1 = plot((1:trials)./trials,log10(sum(sqrt((x_short(:,:,1)-x_short(:,:,2)).^2),2)));
hold on 
h2 = plot((1:trials)./trials,log10(sum(sqrt((x_short(:,:,3)-x_short(:,:,2)).^2),2)));
h3 = plot((1:trials)./trials,log10(sum(sqrt((x_short(:,:,3)-x_short(:,:,1)).^2),2)));
title('Comparison of found portfolios with shorting')
xlabel('alpha')
ylabel('Mean Squared Error, log10')
legend([h1,h2,h3],{'MSE log10, quadprog and our','MSE log10, CVX and our','MSE log10, CVX and quadprog'}, 'Location', 'northeast')

%And when solving the full QP problem
figure
h1 = plot((1:trials)./trials,log10(sum(sqrt((x_nonshort(:,:,1)-x_nonshort(:,:,2)).^2),2)));
hold on 
h2 = plot((1:trials)./trials,log10(sum(sqrt((x_nonshort(:,:,3)-x_nonshort(:,:,2)).^2),2)));
h3 = plot((1:trials)./trials,log10(sum(sqrt((x_nonshort(:,:,3)-x_nonshort(:,:,1)).^2),2)));
title('Comparison of found portfolios with no shorting')
xlabel('alpha')
ylabel('Mean Squared Error, log10')
legend([h1,h2,h3],{'MSE log10, quadprog and our','MSE log10, CVX and our','MSE log10, CVX and quadprog'}, 'Location', 'northeast')


%% Exercise 5.5-5.7, Plot the efficient frontiers
%Plot efficient frontiers and portfolios when allowing or disallowing
%shorting, solving the Bicriterion
clear log

trials = 1000;
x_short = zeros(trials,5,3);
x_nonshort = zeros(trials,5,3);

alphas = 1:1:trials;
alphas = alphas./(trials);
port_risk_short = zeros(trials,1,3);
port_risk_nonshort = zeros(trials,1,3);
port_return_short = zeros(trials,1,3);
port_return_nonshort = zeros(trials,1,3);

x0 = zeros(5,1);
s0 = ones(2*5,1);
y0 = ones(length(beq),1);
z0 = ones(2*5,1);
n = 5;
if trials<11
    cvx_solve = true;
else
    cvx_solve = false;
end

%optimal port. weights found by quadprog
for i = 1:trials
    if mod(i,trials/10)==0
        disp(i)
    end
    % Without shorting, i.e. an IQP problem
    x_nonshort(i,:,1) = quadprog( alphas(i).*H, (1-alphas(i)).*f', Aineq, bineq, Aeq, beq,[],[],[],options); 
    x_nonshort(i,:,2) = primalDualInteriorMethod_box(alphas(i).*H, (1-alphas(i)).*f',Aeq',beq,zeros(5,1),ones(5,1),x0,y0,z0,s0);
    if cvx_solve
        cvx_begin quiet
            %cvx_precision low
            variable x(n)
            minimize( 1/2 * x' *  (alphas(i).*H) *x + ((1-alphas(i)).*f)*x)
            subject to
                Aeq * x == beq
                zeros(5,1) <= x
                x <= ones(5,1)
        cvx_end
        x_nonshort(i,:,3) = x; 
    end
    
    % With shorting, i.e. an EQP problem
    x_short(i,:,1) = quadprog( alphas(i).*H, (1-alphas(i)).*f', [], [], Aeq, beq,[],[],[],options); 
    x_short(i,:,2) = EqualityQPSolver(alphas(i).*H, (1-alphas(i)).*f',Aeq', beq, "rangespace");
    if cvx_solve
        cvx_begin quiet
            %cvx_precision low
            variable x(n)
            minimize( 1/2 * x' *  (alphas(i).*H) *x + ((1-alphas(i)).*f)*x)
            subject to
                Aeq * x == beq
        cvx_end
        x_short(i,:,3) = x; 
    end
    
    port_risk_short(i,2) = x_short(i,:,2)*covariance*x_short(i,:,2)';
    port_risk_nonshort(i,2) = x_nonshort(i,:,2)*covariance*x_nonshort(i,:,2)';
    
    port_return_short(i,2) = -f*x_short(i,:,2)';
    port_return_nonshort(i,2) = -f*x_nonshort(i,:,2)';
end

%Plot the efficient frontier of the bi-criterion problem when shorting is
%allowed
figure
hold on
plot(port_return_short(:,2),port_risk_short(:,2),'b', returns, diag(covariance),'ro')
title('The Bi-Criterion when shorting is allowed')
xlabel('Return [%]')
ylabel('Risk [Var]')
hold off

%And the efficient frontier of the Bi-Criterion problem when shorting is
%not allowed
figure
hold on
plot(port_return_nonshort(:,2),port_risk_nonshort(:,2),'b', returns, diag(covariance),'ro')
title('The Bi-Criterion when shorting is not allowed')
xlabel('Return [%]')
ylabel('Risk [Var]')
hold off

% Plot comparing the solutions to the Bi-Criterion with no shorting and
% optimizing for risk given a return.
figure
hold on
h1 = plot(Rs,port_risks,'r');
plot(returns, diag(covariance),'ro')
h2 = plot(port_return_nonshort(:,2),port_risk_nonshort(:,2),'b');
plot(opt_return1,opt_risk1 ,'k|','MarkerSize', 14)
title('Comparison of static return and bi-criterion')
xlabel('Return [%]')
ylabel('Risk [Var]')
legend([h1,h2],{'Static return','Bi-criterion'}, 'Location', 'northwest')
hold off


%% Exercise 5.8-5.11, Introducing a risk free security with return 0
%Setup for the problem
returns = [16.1, 8.5, 15.7, 10.02, 18.68]; 
covariance = [2.5 .93 .62 .74 -.23;
              0.93 1.5 0.22 0.56 0.26;
              .62  .22 1.9  .78  -0.27;
              .74 .56 .78 3.6 -0.56;
              -0.23 0.26 -0.27 -0.56 3.9];  
          
returns = [returns, 0];

covariance = [covariance, zeros(5,1)];
covariance = [covariance; zeros(1,6)];

H = covariance;
f = []; 

A1 = returns;
b1 = min(A1);

A2 = [1,1,1,1,1,1];
b2 = 1;

Aeq = [A1; A2];
beq = [b1; b2];

Aineq = -eye(6);
bineq = zeros(6,1); 
%% Exercise 5.8-5.11, Computing the efficient frontier with a risk-free asset 

Rs_free = min(returns):0.01:max(returns);
port_risks_free = zeros(length(Rs_free),1);
optimalPorts = zeros(length(Rs),6);
%Just compute it with quadprog
for i = 1:length(Rs_free)
beq = [Rs_free(i); b2];
x = quadprog( H, f, Aineq, bineq, Aeq, beq,[],[],[],options); 
port_risks_free(i) = x'*covariance*x;
optimalPorts(i,:) = x;
end

%% Exercise 5.8-5.11, Plotting the efficient frontier with a risk-free asset 
%Then plot the efficient frontier with a risk-free asset
figure
hold on
plot(Rs_free,port_risks_free,'b', returns, diag(covariance),'ro')
title('The Efficient Frontier with a risk free asset')
xlabel('Return [%]')
ylabel('Risk [Var]')
hold off

lower_risk_pareto = port_risks(min(find(abs(round(port_risks_free(851:end)-port_risks,2))==0)));
lower_return_pareto = Rs(min(find(abs(round(port_risks_free(851:end)-port_risks,2))==0)));

figure
hold on
h1 = plot(Rs_free,port_risks_free,'b');
scatter( returns, diag(covariance),'ro')
h3 = plot(Rs,port_risks,'color',[0 0.5 0]);
plot(opt_return1, opt_risk1,'ko','MarkerFaceColor', 'k', 'MarkerSize',3);
plot(lower_return_pareto,lower_risk_pareto,'ko','MarkerFaceColor', 'k', 'MarkerSize',3);
h2 = plot(effRs,effRisk , 'r');
hold off
title('The Efficient Frontiers with and without a risk free asset')
xlabel('Return [%]')
ylabel('Risk [Var]')
legend([h1,h2, h3],{'With risk free asset','Without risk free asset and pareto optimal','Without risk free asset and not pareto optimal'}, 'Location', 'northwest')

%% Exercise 5.8-5.11, Portfolio compositions
%Plot the optimal portfolio when we have a risk-free security
figure
h1 = plot(Rs_free,optimalPorts(:,1));
hold on
h2 = plot(Rs_free,optimalPorts(:,2));
h3 = plot(Rs_free,optimalPorts(:,3));
h4 = plot(Rs_free,optimalPorts(:,4));
h5 = plot(Rs_free,optimalPorts(:,5));
h6 = plot(Rs_free,optimalPorts(:,6));
legend([h1,h2,h3,h4,h5, h6], {'Security 1','Security 2','Security 3','Security 4','Security 5', 'Risk free asset'})
xlabel('Return [%]')
ylabel('Percentage of the portfolio [%]')
title('Portfolio composition')


%% Exercise 5.8-5.11, Finding the optimum for R=14

beq = [14; b2];
x = quadprog( H, f, Aineq, bineq, Aeq, beq); 
opt_risk2 = x'*covariance*x;


%% Exercise 5.8-5.11, Ploting
%Plot the location R=14 with and without a risk-free asset, so we can
%compare
figure
hold on
opt_return2 = 14;
plot(Rs_free,port_risks_free,'b', returns, diag(covariance),'mo')
plot(Rs,port_risks,'color',[0 0.5 0])
h2 = plot(effRs,effRisk , 'r');
p1 = plot(opt_return2, opt_risk2,'ko','MarkerFaceColor', 'g', 'MarkerSize',4);
plot(opt_return1, opt_risk1,'k|','MarkerFaceColor', 'k', 'MarkerSize',8);
p2 = plot(Rs(find(Rs == 14)), port_risks(find(Rs == 14)),'ko','MarkerFaceColor', 'k', 'MarkerSize',4);
title('Comparison of Efficient Frontiers')
xlabel('Return [%]')
ylabel('Risk [Var]')
legend([p1,p2],{'E[R]=14 with risk free asset','E[R]=14 without risk free asset'},'Location','northwest') 
hold off



