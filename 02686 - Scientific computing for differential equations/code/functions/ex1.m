%% Explicit Euler
A = [0];
b = [1]';
c = [0]';
T = 6;
x0 = 1;
h = 0.6;
t0 = 0;
lambda = -1;
t_a = t0:0.01:T;
[x,t,function_calls,hs] = explicitRungeKutta(@testEq,lambda,h,t0,T,x0,A,b,c);

x_true = exp(lambda*t);
x_true_a = exp(lambda*t_a);
local_error = exp(lambda*h)*x(1:end-1)-x(2:end) ;
global_error = x_true'-x;

% Top two plots
fig1 = figure;
pos_fig1 = [500 300 800 600];
set(fig1,'Position',pos_fig1)

tiledlayout(2,2)
nexttile([1 2])
plot(t,x,t_a,x_true_a)
legend("Explicit Euler", "Analytical Solution",'FontSize',12,'Interpreter','latex')
xlabel("t",'FontSize',12,'Interpreter','latex')
ylabel("x(t)",'FontSize',12,'Interpreter','latex')
title("Explicit Euler",'FontSize',16,'Interpreter','latex')

nexttile
plot(t(1:end-1),local_error, 'o',"color",	'#EDB120')
hold on 
plot(t_a,x0*(R(lambda,h,b,A).^(t_a/h)*exp(h*lambda)-R(lambda,h,b,A).^((t_a+h)/h)))
legend("Emperical error", "Analytical error",'FontSize',12,'Interpreter','latex')
xlabel("t",'FontSize',12,'Interpreter','latex')
ylabel("e(h)",'FontSize',12,'Interpreter','latex')
title("Local truncation error",'FontSize',16,'Interpreter','latex')

nexttile
plot(t,global_error, 'o',"color", '#77AC30')
hold on 
plot(t_a,x0*(exp(t_a*lambda)-R(lambda,h,b,A).^((t_a)/h)))
legend("Emperical error", "Analytical error",'FontSize',12,'Interpreter','latex')
xlabel("t",'FontSize',12,'Interpreter','latex')
ylabel("E(h)",'FontSize',12,'Interpreter','latex')
title("Global truncation error",'FontSize',16,'Interpreter','latex')

%% Implicit Euler
A = [1];
b = [1]';
c = [1]';
T = 6;
x0 = 1;
h = 0.6;
t0 = 0;
lambda = -1;
t_a = t0:0.01:T;
options = struct('step_control',false, 'initialStepSize', false, 'Jac' , @testEqJac);
[x,t,function_calls,hs] = ODEsolver(@testEq,lambda,h,t0,T,x0, "Implicit Euler",options);

x_true = exp(lambda*t);
x_true_a = exp(lambda*t_a);
local_error = exp(lambda*h)*x(1:end-1)-x(2:end) ;
global_error = x_true-x;

% Top two plots
fig1 = figure;
pos_fig1 = [500 300 800 600];
set(fig1,'Position',pos_fig1)

tiledlayout(2,2)
nexttile([1 2])
plot(t,x,t_a,x_true_a)
legend("Implicit Euler", "Analytical Solution",'FontSize',12,'Interpreter','latex')
xlabel("t",'FontSize',12,'Interpreter','latex')
ylabel("x(t)",'FontSize',12,'Interpreter','latex')
title("Implicit Euler",'FontSize',16,'Interpreter','latex')

nexttile
plot(t(1:end-1),local_error, 'o',"color",	'#EDB120')
hold on 
plot(t_a,x0*(R(lambda,h,b,A).^(t_a/h)*exp(h*lambda)-R(lambda,h,b,A).^((t_a+h)/h)))
legend("Emperical error", "Analytical error",'Location',"southeast",'FontSize',12,'Interpreter','latex')
xlabel("t",'FontSize',12,'Interpreter','latex')
ylabel("e(h)",'FontSize',12,'Interpreter','latex')
title("Local truncation error",'FontSize',16,'Interpreter','latex')

nexttile
plot(t,global_error, 'o',"color", '#77AC30')
hold on 
plot(t_a,x0*(exp(t_a*lambda)-R(lambda,h,b,A).^((t_a)/h)))
legend("Emperical error", "Analytical error",'Location',"southeast",'FontSize',12,'Interpreter','latex')
xlabel("t",'FontSize',12,'Interpreter','latex')
ylabel("E(h)",'FontSize',12,'Interpreter','latex')
title("Global truncation error",'FontSize',16,'Interpreter','latex')

%% RK4
A = [0 0 0 0; 0.5 0 0 0; 0 0.5 0 0; 0 0 1 0];
b = [1/6 2/6 2/6 1/6]';
c = [0 1/2 1/2 1]';
T = 6;
x0 = 1;
h = 0.6;
t0 = 0;
lambda = -1;
t_a = t0:0.01:T;
[x,t,function_calls,hs] = explicitRungeKutta(@testEq,lambda,h,t0,T,x0,A,b,c);

x_true = exp(lambda*t);
x_true_a = exp(lambda*t_a);
local_error = exp(lambda*h)*x(1:end-1)-x(2:end) ;
global_error = x_true'-x;

ana_loc = x0*(R(lambda,h,b,A).^(t/h)*exp(h*lambda)-R(lambda,h,b,A).^((t+h)/h));

% Top two plots
fig1 = figure;
pos_fig1 = [500 300 800 600];
set(fig1,'Position',pos_fig1)

tiledlayout(2,2)
nexttile([1 2])
plot(t,x,t_a,x_true_a)
legend("Classic Runge-Kutta", "Analytical Solution",'FontSize',12,'Interpreter','latex')
xlabel("t",'FontSize',12,'Interpreter','latex')
ylabel("x(t)",'FontSize',12,'Interpreter','latex')
title("Classic Runge-Kutta",'FontSize',16,'Interpreter','latex')

nexttile
plot(t(1:end-1),local_error, 'o',"color",	'#EDB120')
hold on 
plot(t_a,x0*(R(lambda,h,b,A).^(t_a/h)*exp(h*lambda)-R(lambda,h,b,A).^((t_a+h)/h)))
legend("Emperical error", "Analytical error",'Location',"southeast",'FontSize',12,'Interpreter','latex')
xlabel("t",'FontSize',12,'Interpreter','latex')
ylabel("e(h)",'FontSize',12,'Interpreter','latex')
title("Local truncation error",'FontSize',16,'Interpreter','latex')

nexttile
plot(t,global_error, 'o',"color", '#77AC30')
hold on 
plot(t_a,x0*(exp(t_a*lambda)-R(lambda,h,b,A).^((t_a)/h)))
legend("Emperical error", "Analytical error",'Location',"southeast",'FontSize',12,'Interpreter','latex')
xlabel("t",'FontSize',12,'Interpreter','latex')
ylabel("E(h)",'FontSize',12,'Interpreter','latex')
title("Gloal truncation error",'FontSize',16,'Interpreter','latex')

%% Plot of local error vs step size
x0 = 1;
delta = 0.01;
zoom = 0.2;
h = 0:delta:5;
t0 = 0;
lambda = -1;

% RK4
A = [0 0 0 0; 0.5 0 0 0; 0 0.5 0 0; 0 0 1 0];
b = [1/6 2/6 2/6 1/6]';
c = [0 1/2 1/2 1]';
tau_RK4 = x0*(exp(h*lambda)-arrayfun(@(h_i) R(lambda,h_i,b,A),h));

% Implicit
A = [1];
b = [1]';
c = [1]';
tau_Imp = x0*(exp(h*lambda)-arrayfun(@(h_i) R(lambda,h_i,b,A),h));

% Explicit
A = [0];
b = [1]';
c = [0]';
tau_Exp = x0*(exp(h*lambda)-arrayfun(@(h_i) R(lambda,h_i,b,A),h));



fig1 = figure;
pos_fig1 = [500 300 800 600];
set(fig1,'Position',pos_fig1)
tiledlayout(2,2)
nexttile
plot(h,tau_RK4,h,tau_Imp,h,tau_Exp)
ylim([min(min(min(tau_RK4,tau_Imp),tau_Exp))*1.3,max(max(max(tau_RK4,tau_Imp),tau_Exp))*1.3])
legend("Classic Runge-Kutta", "Implicit Euler", "Explicit Euler",'Location','southwest','FontSize',10,'Interpreter','latex')
xlabel("h",'FontSize',12,'Interpreter','latex')
ylabel("e(h)",'FontSize',12,'Interpreter','latex')
title("Behaviour for large step sizes",'FontSize',16,'Interpreter','latex')


nexttile
plot(h,abs(tau_RK4),h,abs(tau_Imp),h,abs(tau_Exp))
ylim([0,max(max(max(tau_RK4,tau_Imp),tau_Exp))*1.3])
legend("Classic Runge-Kutta", "Implicit Euler", "Explicit Euler",'Location','northwest','FontSize',10,'Interpreter','latex')
xlabel("h",'FontSize',12,'Interpreter','latex')
ylabel("$|e(h)|$",'FontSize',12,'Interpreter','latex')
title("Behaviour for large step sizes",'FontSize',16,'Interpreter','latex')

nexttile
plot(h(1:zoom/delta+1),tau_RK4(1:zoom/delta+1),h(1:zoom/delta+1),tau_Imp(1:zoom/delta+1),h(1:zoom/delta+1),tau_Exp(1:zoom/delta+1))
ylim([min(min(min(tau_RK4(1:zoom/delta+1),tau_Imp(1:zoom/delta+1)),tau_Exp(1:zoom/delta+1)))*1.3,max(max(max(tau_RK4(1:zoom/delta+1),tau_Imp(1:zoom/delta+1)),tau_Exp(1:zoom/delta+1)))*1.3])
legend("Classic Runge-Kutta", "Implicit Euler", "Explicit Euler",'Location','northwest','FontSize',10,'Interpreter','latex')
xlabel("h",'FontSize',12,'Interpreter','latex')
ylabel("e(h)",'FontSize',12,'Interpreter','latex')
title("Behaviour for small step sizes",'FontSize',16,'Interpreter','latex')

nexttile
plot(h(1:zoom/delta+1),abs(tau_RK4(1:zoom/delta+1)),h(1:zoom/delta+1),abs(tau_Imp(1:zoom/delta+1)),h(1:zoom/delta+1),abs(tau_Exp(1:zoom/delta+1)))
ylim([0,max(max(max(tau_RK4(1:zoom/delta+1),tau_Imp(1:zoom/delta+1)),tau_Exp(1:zoom/delta+1)))*1.3])
legend("Classic Runge-Kutta", "Implicit Euler", "Explicit Euler",'Location','northwest','FontSize',10,'Interpreter','latex')
xlabel("h",'FontSize',12,'Interpreter','latex')
ylabel("$|e(h)|$",'FontSize',12,'Interpreter','latex')
title("Behaviour for small step sizes",'FontSize',16,'Interpreter','latex')

%% Global error for T=1 vs step size(h)
x0 = 1;
delta = 0.02;
zoom = 0.199;
h = 1*(0.0000001:delta:0.999999999);
T = 1*1;
lambda = -1;

% RK4
A = [0 0 0 0; 0.5 0 0 0; 0 0.5 0 0; 0 0 1 0];
b = [1/6 2/6 2/6 1/6]';
c = [0 1/2 1/2 1]';

tau_RK4 = x0*(exp(T*lambda)-arrayfun(@(h_i) R(lambda,h_i,b,A)^((T)./h_i),h));

% Implicit
A = [1];
b = [1]';
c = [1]';
tau_Imp = x0*(exp(T*lambda)-arrayfun(@(h_i) R(lambda,h_i,b,A)^((T)./h_i),h));

% Explicit
A = [0];
b = [1]';
c = [0]';
tau_Exp = x0*(exp(T*lambda)-arrayfun(@(h_i) R(lambda,h_i,b,A)^((T)/h_i),h));
tau_Exp = abs(tau_Exp);


fig1 = figure;
pos_fig1 = [500 300 800 600];
set(fig1,'Position',pos_fig1)
tiledlayout(2,2)
nexttile
plot(h,tau_RK4,h,tau_Imp,h,tau_Exp)
ylim([min(min(min(tau_RK4,tau_Imp),tau_Exp))*1.3,max(max(max(tau_RK4,tau_Imp),tau_Exp))*1.3])
legend("Classic Runge-Kutta", "Implicit Euler", "Explicit Euler",'Location','northwest','FontSize',10,'Interpreter','latex')
xlabel("h",'FontSize',12,'Interpreter','latex')
ylabel("E(h)",'FontSize',12,'Interpreter','latex')
title("Global error at T=1.0 for large step sizes",'FontSize',14,'Interpreter','latex')


nexttile
plot(h,abs(tau_RK4),h,abs(tau_Imp),h,abs(tau_Exp))
ylim([0,max(max(max(tau_RK4,tau_Imp),tau_Exp))*1.3])
legend("Classic Runge-Kutta", "Implicit Euler", "Explicit Euler",'Location','northwest','FontSize',10,'Interpreter','latex')
xlabel("h",'FontSize',12,'Interpreter','latex')
ylabel("$|E(h)|$",'FontSize',12,'Interpreter','latex')
title("Global error at T=1.0 for large step sizes",'FontSize',14,'Interpreter','latex')

nexttile
plot(h(1:zoom/delta+1),tau_RK4(1:zoom/delta+1),h(1:zoom/delta+1),tau_Imp(1:zoom/delta+1),h(1:zoom/delta+1),tau_Exp(1:zoom/delta+1))
ylim([min(min(min(tau_RK4(1:zoom/delta+1),tau_Imp(1:zoom/delta+1)),tau_Exp(1:zoom/delta+1)))*1.3,max(max(max(tau_RK4(1:zoom/delta+1),tau_Imp(1:zoom/delta+1)),tau_Exp(1:zoom/delta+1)))*1.3])
legend("Classic Runge-Kutta", "Implicit Euler", "Explicit Euler",'Location','northwest','FontSize',10,'Interpreter','latex')
xlabel("h",'FontSize',12,'Interpreter','latex')
ylabel("E(h)",'FontSize',12,'Interpreter','latex')
title("Global error at T=1.0 for small step sizes",'FontSize',14,'Interpreter','latex')

nexttile
plot(h(1:zoom/delta+1),abs(tau_RK4(1:zoom/delta+1)),h(1:zoom/delta+1),abs(tau_Imp(1:zoom/delta+1)),h(1:zoom/delta+1),abs(tau_Exp(1:zoom/delta+1)))
ylim([0,max(max(max(tau_RK4(1:zoom/delta+1),tau_Imp(1:zoom/delta+1)),tau_Exp(1:zoom/delta+1)))*1.3])
legend("Classic Runge-Kutta", "Implicit Euler", "Explicit Euler",'Location','northwest','FontSize',10,'Interpreter','latex')
xlabel("h",'FontSize',12,'Interpreter','latex')
ylabel("$|E(h)|$",'FontSize',12,'Interpreter','latex')
title("Global error at T=1.0 for small step sizes",'FontSize',14,'Interpreter','latex')

%% A stability (half plane)
A = [0];
b = [1]';
c = [0]';
T = 80;
x0 = [1;0];
h = 4;
t0 = 0;
lambda = 0;
[x,t,function_calls,hs] = explicitRungeKutta(@Astable_half,lambda,h,t0,T,x0,A,b,c);

options = struct('step_control',false, 'initialStepSize', false, 'Jac' , @Astable_half_Jac);
[x_imp,t,function_calls,hs] = ODEsolver(@Astable_half,lambda,h,t0,T,x0, "Implicit Euler",options);

[t,x_itra] = itrapez(@Astable_half,[t0 T],x0,T/h,@Astable_half_Jac,1e-4,100);

h = 0.1;
[t,x_itra_high] = itrapez(@Astable_half,[t0 T],x0,T/h,@Astable_half_Jac,1e-4,100);

%figure
%plot(x(:,1),x(:,2))
%xlim([-2.5,2.5])
%ylim([-2.5,2.5])

figure 
plot(x_imp(:,1),x_imp(:,2),'o')
hold on 
plot(x_itra_high(1,:),x_itra_high(2,:))
plot(x0(1),x0(2), 'x', 'MarkerSize', 20, 'color', 'r')
xlim([-1.3,1.7])
ylim([-1.3,1.3])
xlabel("$x_1$",'FontSize',14,'Interpreter','latex')
ylabel("$x_2$",'FontSize',14,'Interpreter','latex')
legend("Simulation", "Analytical solution", "x(0)=x0",'FontSize',12,'Interpreter','latex')

figure
plot(x_itra(1,:),x_itra(2,:), 'o', 'color', 'b')
hold on
plot(x_itra_high(1,:),x_itra_high(2,:))
plot(x0(1),x0(2), 'x', 'MarkerSize', 20, 'color', 'r')
xlim([-1.3,1.7])
ylim([-1.3,1.3])
xlabel("$x_1$",'FontSize',14,'Interpreter','latex')
ylabel("$x_2$",'FontSize',14,'Interpreter','latex')
legend("Simulation", "Analytical solution", "x(0)=x0",'AutoUpdate','off','FontSize',12,'Interpreter','latex')
plot(x_itra(1,:),x_itra(2,:), 'color', 'b')

%% L stability
A = [0];
b = [1]';
c = [0]';
T = 2*1.5;
x0 = [0];
h = 1.5/40;
t0 = 0;
lambda = 0;
[t_itra,x_itra] = itrapez(@Lstable,[t0 T],x0,T/h,@LstableJac,1e-4,100);

options = struct('step_control',false, 'initialStepSize', false, 'Jac' , @LstableJac);
[x_imp,t,function_calls,hs] = ODEsolver(@Lstable,lambda,h,t0,T,x0, "Implicit Euler",options);

figure
plot(t,x_imp,t_itra,x_itra)
hold on
plot(t0,x0, 'x', 'MarkerSize', 10, 'color', 'r')
xlim([-0.1,T])
xlabel("t",'FontSize',14,'Interpreter','latex')
ylabel("x(t)",'FontSize',14,'Interpreter','latex')
legend("Implicit Euler", "Trapezoidal", "x(0)=x0",'FontSize',12,'Interpreter','latex')


%%
precision = 0.02;
xlimit = [-4,4];
ylimit = [-4,4];

% Explicit Euler
A = [0];
b = [1]';
c = [0]';
stability(A,b,c,xlimit,ylimit,precision,"Explicit Euler")

%%

precision = 0.02;
xlimit = [-4,4];
ylimit = [-4,4];

% Implicit Euler
A = [1];
b = [1]';
c = [1]';
stability(A,b,c,xlimit,ylimit,precision,"Implicit Euler")


%%

precision = 0.02;
xlimit = [-4,4];
ylimit = [-4,4];

% RK4
A = [0 0 0 0; 0.5 0 0 0; 0 0.5 0 0; 0 0 1 0];
b = [1/6 2/6 2/6 1/6]';
c = [0 1/2 1/2 1]';
stability(A,b,c,xlimit,ylimit,precision,"Classic Runge-Kutta")

%%

precision = 0.02;
xlimit = [-4,4];
ylimit = [-4,4];

% RK4
A = [0 0; 1/2 1/2];
b = [1/2 1/2]';
c = [0 1]';
stability(A,b,c,xlimit,ylimit,precision,"Trapezoidal")

%% 
precision = 0.02;
xlimit = [-6,6];
ylimit = [-6,6];

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
stability(A,b,c,xlimit,ylimit,precision,"Dormand-Prince, 5th order")
%stabilityEmb(A,b,bhat,xlimit,ylimit,precision,"Dormand-Prince 5(4)")


%% 
precision = 0.04;
xlimit = [-20,20];
ylimit = [-20,20];

gamma = (2-sqrt(2))/2;
c =[0 2*gamma 1]';
A = [0 0 0;
    gamma gamma 0;
    (1-gamma)/2 (1-gamma)/2 gamma];

b = [(1-gamma)/2 (1-gamma)/2 gamma]';
bhat = [(6*gamma-1)/(12*gamma) 1/(12*gamma*(1-2*gamma)) (1-3*gamma)/(3*(1-2*gamma))]';
%stability(A,b,c,xlimit,ylimit,precision,"Dormand-Prince, 5th order")
%stabilityEmb(A,b,bhat,xlimit,ylimit,precision,"ESDIRK23")

stability(A,b,c,xlimit,ylimit,precision,"ESDIRK23")
%stability(A,bhat,c,xlimit,ylimit,precision,"Embedded")
%%
gamma = (2+sqrt(2))/2;
c =[0 2*gamma 1]';
A = [0 0 0;
    gamma gamma 0;
    (1-gamma)/2 (1-gamma)/2 gamma];

b = [(1-gamma)/2 (1-gamma)/2 gamma]';

seq_num = 0:10000:100000000;
seq_num = [-seq_num -Inf];
Lvals = arrayfun(@(x) R(x,1,b,A),seq_num);

figure
semilogy(seq_num,abs(Lvals))
%%
function [] = stability(A,b,c,xlimit,ylimit, precision,titletext)
    s = length(b);
    re = [xlimit(1):precision:xlimit(2)];
    im = [ylimit(1):precision:ylimit(2)];
    abs_matrix = zeros(length(re),length(im));
    fprintf('\n');
    reverseStr = '';
    for i = 1:length(re)
        for j = 1:length(im)
            abs_matrix(i,j) = abs(1+complex(re(i),im(j))*b'*inv(eye(s)-complex(re(i),im(j))*A)*ones(s,1));

        end
       percentDone = 100 * i / length(re);
       msg = sprintf('Percent done: %3.1f', percentDone);
       fprintf([reverseStr, msg]);
       reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    fprintf('\n');
    plot_matrix = abs_matrix;
    %plot_matrix = (abs_matrix > 1);
    plot_matrix(plot_matrix>1) = 1.0;
    %plot_matrix(round(plot_matrix,1)==1.0) = 0;

    figure
    hold on 
    %map = jet(256);
    %map(end,:) = [1,1,1]';
    map = [repmat([1 0 0],1000,1);[1,1,1]];
    colormap(map)
    %colorbar;
    stabMap = imagesc(re',im',plot_matrix');
    alpha(stabMap,0.6);
    grid on
    plot(re,zeros(1,length(re)),'k', 'LineWidth', 0.5)
    plot(zeros(1,length(im)),im,'k', 'LineWidth', 0.5)
    xlim([min(re) max(re)])
    ylim([min(im) max(im)])
    title(titletext,'FontSize',16,'Interpreter','latex')
    xlabel("Re($\lambda h$)",'FontSize',14,'Interpreter','latex')
    ylabel("Im($\lambda h$)",'FontSize',14,'Interpreter','latex')
    hgh = plot(1000,10000,'o','color','r',"MarkerFaceColor", "r" );
    legend(hgh, "Stable",'FontSize',12,'Interpreter','latex')


    if all(all(abs_matrix(re<0,:)<1))
        disp('The method is A-stable')
        if all(all(abs_matrix(re>0,:)>1))
            disp('Only Re(z)<0 behaves stable')
        end
        if b'*inv(A)*ones(s,1)-1 < 0.000001
            disp('The method is L-stable')
        end
    else
        disp('The method is not A-stable')
    end

end

function [] = stabilityEmb(A,b,b2,xlimit,ylimit, precision,titletext)
    s = length(b);
    re = [xlimit(1):precision:xlimit(2)];
    im = [ylimit(1):precision:ylimit(2)];
    abs_matrix = zeros(length(re),length(im));
    fprintf('\n');
    reverseStr = '';
    for i = 1:length(re)
        for j = 1:length(im)
            abs_matrix(i,j) = abs(1+complex(re(i),im(j))*b'*inv(eye(s)-complex(re(i),im(j))*A)*ones(s,1));

        end
       percentDone = 100 * i / length(re);
       msg = sprintf('Percent done: %3.1f', percentDone);
       fprintf([reverseStr, msg]);
       reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    fprintf('\n');
    plot_matrix = abs_matrix;
    %plot_matrix = (abs_matrix > 1);
    plot_matrix(plot_matrix>1) = 1.0;
    
    
    
    s = length(b2);
    re = [xlimit(1):precision:xlimit(2)];
    im = [ylimit(1):precision:ylimit(2)];
    abs_matrix = zeros(length(re),length(im));
    fprintf('\n');
    reverseStr = '';
    for i = 1:length(re)
        for j = 1:length(im)
            abs_matrix(i,j) = abs(1+complex(re(i),im(j))*b2'*inv(eye(s)-complex(re(i),im(j))*A)*ones(s,1));

        end
       percentDone = 100 * i / length(re);
       msg = sprintf('Percent done: %3.1f', percentDone);
       fprintf([reverseStr, msg]);
       reverseStr = repmat(sprintf('\b2'), 1, length(msg));
    end
    fprintf('\n');
    plot_matrix2 = abs_matrix;
    %plot_matrix = (abs_matrix > 1);
    plot_matrix2(plot_matrix2>1) = 1.0;
    
    %plot_matrix(round(plot_matrix,1)==1.0) = 0;

    fig1 = figure;
    pos_fig1 = [500 300 800 600];
    set(fig1,'Position',pos_fig1)
    hold on 
    %map = jet(256);
    %map(end,:) = [1,1,1]';
    map = [repmat([0 1 0],1000,1);[1,1,1]];
    colormap(map)
    %colorbar;
    stabMap = imagesc(re',im',plot_matrix');
    stabMap2 = imagesc(re',im',plot_matrix2');
    alpha(stabMap,0.6);
    alpha(stabMap2,0.2);
    grid on
    plot(re,zeros(1,length(re)),'k', 'LineWidth', 0.5)
    plot(zeros(1,length(im)),im,'k', 'LineWidth', 0.5)
    xlim([min(re) max(re)])
    ylim([min(im) max(im)])
    title(titletext,'FontSize',16,'Interpreter','latex')
    xlabel("Re($\lambda h$)",'FontSize',14,'Interpreter','latex')
    ylabel("Im($\lambda h$)",'FontSize',14,'Interpreter','latex')
    hgh = plot(1000,10000,'o','color','g',"MarkerFaceColor", "g" );
    hgh2 = plot(1000,10000,'o','color',[0.7 1 0.7],"MarkerFaceColor", [0.7 1 0.7] );
    legend([hgh,hgh2], {'Method', 'Embeded'},'FontSize',12,'Interpreter','latex')


    if all(all(abs_matrix(re<0,:)<1))
        disp('The method is A-stable')
        if all(all(abs_matrix(re>0,:)>1))
            disp('Only Re(z)<0 behaves stable')
        end
        if b'*inv(A)*ones(s,1)-1 < 0.000001
            disp('The method is L-stable')
        end
    else
        disp('The method is not A-stable')
    end

end


function[tau] = R(lambda,h,b,A)
    tau = (1+lambda*h*(b)'*inv(eye(length(b))-lambda*h*A)*ones(length(b),1));
end




