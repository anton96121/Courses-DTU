clear all

%% Stability
precision = 0.1;
xlimit = [-20,20];
ylimit = [-20,20];

% ESDIRK23
gamma = 1-1/sqrt(2);
        a31 = (1-gamma)/2;
        AT = [0 gamma a31;0 gamma a31;0 0 gamma];
        c  = [0; 2*gamma; 1];
        b  = AT(:,3);
stability(AT',b,c,xlimit,ylimit,precision,"ESDIRK23");


bhat = [    (6*gamma-1)/(12*gamma); ...
    1/(12*gamma*(1-2*gamma)); ...
    (1-3*gamma)/(3*(1-2*gamma))    ];

stability(AT',bhat,c,xlimit,ylimit,precision,"ESDIRK23, Embedded");
%% L-stability

% ESDIRK solution
mu = 0;
t0 = 0;
tfinal = 2*1.5;
x0 = [0];
h0 = 0.01;
options = struct('step_control',false, 'initialStepSize', false, 'Jac' , @LstableJac);
[Xout,Tout,function_calls1,hs1] = ODEsolver(@Lstable,[mu],h0,t0,tfinal,x0, "ESDIRK23", options);

t = Tout;
x = Xout;


A = [0];
b = [1]';
c = [0]';
T = 1;
x0 = [0];
h = 1.5/(40*2.25);
t0 = 0;
lambda = 0;
[t_itra,x_itra] = itrapez(@Lstable,[t0 T],x0,T/h,@LstableJac,1e-4,100);
% plot results
figure
plot(t,x,t_itra,x_itra)
hold on
plot(t0,x0, 'x', 'MarkerSize', 10, 'color', 'r')
xlim([-0.1,T])
xlabel("t",'FontSize',14,'Interpreter','latex')
ylabel("x(t)",'FontSize',14,'Interpreter','latex')
legend("ESDIRK23", "Trapezoidal", "x(0)=x0",'FontSize',12,'Interpreter','latex')

%% A stability
mu = 0;
t0 = 0;
tfinal = 8000;
x0 = [1;0];
h0 = 0.01;
options = struct('Jac' , @Astable_half_Jac,'Rtol',0.1,'Atol',0.1);
[Xout,Tout,function_calls1,hs1] = ODEsolver(@Astable_half,[mu],h0,t0,tfinal,x0, "ESDIRK23", options);
t = Tout;
x = Xout;

mu = 0;
t0 = 0;
tfinal = 800;
x0 = [1;0];
h0 = 0.01;
options = struct('Jac' , @Astable_half_Jac,'Rtol',0.0001,'Atol',0.0001);
[Xout,Tout,function_calls1,hs1] = ODEsolver(@Astable_half,[mu],h0,t0,tfinal,x0, "ESDIRK23", options);
t = Tout;
x2 = Xout;

h = 0.1;
[t,x_itra_high] = itrapez(@Astable_half,[t0 tfinal],x0,tfinal/h,@Astable_half_Jac,1e-4,100);


fig1 = figure;
pos_fig1 = [600 400 900 500];
set(fig1,'Position',pos_fig1)
tiledlayout(1,2)

nexttile
hold on
plot(x_itra_high(1,:),x_itra_high(2,:))
plot(x0(1),x0(2), 'x', 'MarkerSize', 20, 'color', 'r')
xlim([-1.3,1.7])
ylim([-1.3,1.3])
xlabel("$x_1$",'FontSize',14,'Interpreter','latex')
ylabel("$x_2$",'FontSize',14,'Interpreter','latex')
title("AbsTol = RelTol = 0.1",'FontSize',16,'Interpreter','latex')
plot(x(:,1),x(:,2), 'o', 'color', 'b')
legend("Analytical solution", "x(0)=x0","ESDIRK23",'AutoUpdate','off','FontSize',12,'Interpreter','latex')

nexttile
hold on
plot(x_itra_high(1,:),x_itra_high(2,:))
plot(x0(1),x0(2), 'x', 'MarkerSize', 20, 'color', 'r')
xlim([-1.3,1.7])
ylim([-1.3,1.3])
xlabel("$x_1$",'FontSize',14,'Interpreter','latex')
ylabel("$x_2$",'FontSize',14,'Interpreter','latex')
title("AbsTol = RelTol = 0.0001",'FontSize',16,'Interpreter','latex')
plot(x2(:,1),x2(:,2), 'o', 'color', 'b')
legend("Analytical solution", "x(0)=x0","ESDIRK23",'AutoUpdate','off','FontSize',12,'Interpreter','latex')




%% Van Pol
%% MU = 3, adaptive step size
h = 0.1;
t0 = 0;
mu = 3;
tend = 4*mu;
x0 = [1.0; 1.0];
Atol = 10^(-2);
Rtol = Atol;

% using the ESDIRK23method
options = struct('Jac' , @VanPolJac,'Rtol',Rtol,'Atol',Atol);
[X1,T1,function_calls7,hs1,rs7] = ODEsolver(@VanPol,[mu],h,t0,tend,x0, "ESDIRK23", options);
hs7 = hs1;

Atol = 10^(-4);
Rtol = Atol;
% using the ESDIRK23method
options = struct('Jac' , @VanPolJac,'Rtol',Rtol,'Atol',Atol);
[X2,T2,function_calls8,hs2,rs8] = ODEsolver(@VanPol,[mu],h,t0,tend,x0, "ESDIRK23", options);
hs8 = hs2;

Atol = 10^(-6);
Rtol = Atol;
% using the ESDIRK23method
options = struct('Jac' , @VanPolJac,'Rtol',Rtol,'Atol',Atol);
[X3,T3,function_calls9,hs3,rs9] = ODEsolver(@VanPol,[mu],h,t0,tend,x0, "ESDIRK23", options);
hs9 = hs3;

options = odeset('RelTol',0.000001,'AbsTol',0.000001);
[Tode45,Xode45] = ode15s(@VanPol,[t0 tend],x0,options,mu);

[Tode15s,Xode15s]=ode45(@VanPol,[t0 tend],x0,options,mu);

%
fig = figure;

subplot(3, 3, 1)
plot(T1,X1(:,1),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,1), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,1), '--','LineWidth', 1.5)
title('Atol = Rtol = $10^{-2}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_1$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 2)
plot(T2,X2(:,1),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,1), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,1), '--','LineWidth', 1.5)
title('Atol = Rtol = $10^{-4}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_1$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 3)
plot(T3,X3(:,1),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,1), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,1), '--','LineWidth', 1.5)
title('Atol = Rtol = $10^{-6}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_1$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 4)
plot(T1,X1(:,2),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,2), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 5)
plot(T2,X2(:,2),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,2), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 6)
plot(T3,X3(:,2),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,2), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 7)
plot(X1(:,1),X1(:,2),'LineWidth', 1.5)
hold on 
plot(Xode45(:,1),Xode45(:,2), '-.','LineWidth', 1.5)
plot(Xode15s(:,1),Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("$x_1$",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 8)
plot(X2(:,1),X2(:,2),'LineWidth', 1.5)
hold on 
plot(Xode45(:,1),Xode45(:,2), '-.','LineWidth', 1.5)
plot(Xode15s(:,1),Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("$x_1$",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 9)
plot(X3(:,1),X3(:,2),'LineWidth', 1.5)
hold on 
plot(Xode45(:,1),Xode45(:,2), '-.','LineWidth', 1.5)
plot(Xode15s(:,1),Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("$x_1$",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')


sgtitle('ESDIRK23 with adaptive step size, $\mu = 3$','Interpreter','latex','FontSize',20) 
 

% add a bit space to the figure
fig = gcf;
%fig.Position(3) = fig.Position(3) + 250;
% add legend
Lgnd = legend('show');
Lgnd.Position(1) = 0.449;
Lgnd.Position(2) = -0.02;
legend("ESDIRK23","ode45","ode15s",'NumColumns',3,'FontSize',16,'Interpreter','latex');


% add a bit space to the figure
fig = gcf;
%fig.Position(3) = fig.Position(3) + 250;
% add legend
Lgnd = legend('show');
Lgnd.Position(1) = 0.449;
Lgnd.Position(2) = -0.02;
legend("ESDIRK23","ode45","ode15s",'NumColumns',3,'FontSize',16,'Interpreter','latex');


fig = figure;
subplot(1, 3, 1)
[hs1, idx_reject] = reject_data(hs1);
semilogy(T1,hs1,'LineWidth', 1.5)
hold on
semilogy(T1(idx_reject),hs1(idx_reject),'x', 'color', 'r','MarkerSize',12)
title('Atol = Rtol = $10^{-2}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("h",'FontSize',16,'Interpreter','latex')
ylim([10^(-5),10])

subplot(1, 3, 2)
[hs2, idx_reject] = reject_data(hs2);
semilogy(T2,hs2,'LineWidth', 1.5)
hold on
semilogy(T2(idx_reject),hs2(idx_reject),'x', 'color', 'r','MarkerSize',12)
title('Atol = Rtol = $10^{-4}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("h",'FontSize',16,'Interpreter','latex')
ylim([10^(-5),10])

subplot(1, 3, 3)
[hs3, idx_reject] = reject_data(hs3);
semilogy(T3,hs3,'LineWidth', 1.5)
hold on
semilogy(T3(idx_reject),hs3(idx_reject),'x', 'color', 'r','MarkerSize',12)
title('Atol = Rtol = $10^{-6}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("h",'FontSize',16,'Interpreter','latex')
ylim([10^(-5),10])


sgtitle('ESDIRK23 with adaptive step size, $\mu = 3$','Interpreter','latex','FontSize',20) 
 

% add a bit space to the figure
fig = gcf;
%fig.Position(3) = fig.Position(3) + 250;
% add legend
Lgnd = legend('show');
Lgnd.Position(1) = 0.449;
Lgnd.Position(2) = 0;
legend("Step size","Rejections",'NumColumns',2,'FontSize',16,'Interpreter','latex');

%% MU = 20, adaptive step size
%% mu = 20, adaptive step size
h = 0.1;
t0 = 0;
mu = 20;
tend = 4*mu;
x0 = [1.0; 1.0];
Atol = 10^(-2);
Rtol = Atol;

% using the ESDIRK23method
options = struct('Jac' , @VanPolJac,'Rtol',Rtol,'Atol',Atol);
[X1,T1,function_calls10,hs1,rs10] = ODEsolver(@VanPol,[mu],h,t0,tend,x0, "ESDIRK23", options);
hs10 = hs1;

Atol = 10^(-4);
Rtol = Atol;
% using the ESDIRK23method
options = struct('Jac' , @VanPolJac,'Rtol',Rtol,'Atol',Atol);
[X2,T2,function_calls11,hs2,rs11] = ODEsolver(@VanPol,[mu],h,t0,tend,x0, "ESDIRK23", options);
hs11 = hs2;

Atol = 10^(-6);
Rtol = Atol;
% using the ESDIRK23method
options = struct('Jac' , @VanPolJac,'Rtol',Rtol,'Atol',Atol);
[X3,T3,function_calls12,hs3,rs12] = ODEsolver(@VanPol,[mu],h,t0,tend,x0, "ESDIRK23", options);
hs12 = hs3;

options = odeset('RelTol',0.000001,'AbsTol',0.000001);
[Tode45,Xode45] = ode15s(@VanPol,[t0 tend],x0,options,mu);

[Tode15s,Xode15s]=ode45(@VanPol,[t0 tend],x0,options,mu);

%
fig = figure;

subplot(3, 3, 1)
plot(T1,X1(:,1),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,1), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,1), '--','LineWidth', 1.5)
title('Atol = Rtol = $10^{-2}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_1$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 2)
plot(T2,X2(:,1),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,1), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,1), '--','LineWidth', 1.5)
title('Atol = Rtol = $10^{-4}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_1$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 3)
plot(T3,X3(:,1),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,1), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,1), '--','LineWidth', 1.5)
title('Atol = Rtol = $10^{-6}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_1$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 4)
plot(T1,X1(:,2),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,2), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 5)
plot(T2,X2(:,2),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,2), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 6)
plot(T3,X3(:,2),'LineWidth', 1.5)
hold on 
plot(Tode45,Xode45(:,2), '-.','LineWidth', 1.5)
plot(Tode15s,Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 7)
plot(X1(:,1),X1(:,2),'LineWidth', 1.5)
hold on 
plot(Xode45(:,1),Xode45(:,2), '-.','LineWidth', 1.5)
plot(Xode15s(:,1),Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("$x_1$",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 8)
plot(X2(:,1),X2(:,2),'LineWidth', 1.5)
hold on 
plot(Xode45(:,1),Xode45(:,2), '-.','LineWidth', 1.5)
plot(Xode15s(:,1),Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("$x_1$",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')

subplot(3, 3, 9)
plot(X3(:,1),X3(:,2),'LineWidth', 1.5)
hold on 
plot(Xode45(:,1),Xode45(:,2), '-.','LineWidth', 1.5)
plot(Xode15s(:,1),Xode15s(:,2), '--','LineWidth', 1.5)
xlabel("$x_1$",'FontSize',16,'Interpreter','latex')
ylabel("$x_2$",'FontSize',16,'Interpreter','latex')


sgtitle('ESDIRK23 with adaptive step size, $\mu = 20$','Interpreter','latex','FontSize',20) 
 

% add a bit space to the figure
fig = gcf;
%fig.Position(3) = fig.Position(3) + 250;
% add legend
Lgnd = legend('show');
Lgnd.Position(1) = 0.449;
Lgnd.Position(2) = -0.02;
legend("ESDIRK23","ode45","ode15s",'NumColumns',3,'FontSize',16,'Interpreter','latex');


% add a bit space to the figure
fig = gcf;
%fig.Position(3) = fig.Position(3) + 250;
% add legend
Lgnd = legend('show');
Lgnd.Position(1) = 0.449;
Lgnd.Position(2) = -0.02;
legend("ESDIRK23","ode45","ode15s",'NumColumns',3,'FontSize',16,'Interpreter','latex');


fig = figure;
subplot(1, 3, 1)
[hs1, idx_reject] = reject_data(hs1);
semilogy(T1,hs1,'LineWidth', 1.5)
hold on
semilogy(T1(idx_reject),hs1(idx_reject),'x', 'color', 'r','MarkerSize',12)
title('Atol = Rtol = $10^{-2}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("h",'FontSize',16,'Interpreter','latex')
ylim([10^(-5),10])

subplot(1, 3, 2)
[hs2, idx_reject] = reject_data(hs2);
semilogy(T2,hs2,'LineWidth', 1.5)
hold on
semilogy(T2(idx_reject),hs2(idx_reject),'x', 'color', 'r','MarkerSize',12)
title('Atol = Rtol = $10^{-4}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("h",'FontSize',16,'Interpreter','latex')
ylim([10^(-5),10])

subplot(1, 3, 3)
[hs3, idx_reject] = reject_data(hs3);
semilogy(T3,hs3,'LineWidth', 1.5)
hold on
semilogy(T3(idx_reject),hs3(idx_reject),'x', 'color', 'r','MarkerSize',12)
title('Atol = Rtol = $10^{-6}$','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("h",'FontSize',16,'Interpreter','latex')
ylim([10^(-5),10])


sgtitle('ESDIRK23 with adaptive step size, $\mu = 20$','Interpreter','latex','FontSize',20) 
 

% add a bit space to the figure
fig = gcf;
%fig.Position(3) = fig.Position(3) + 250;
% add legend
Lgnd = legend('show');
Lgnd.Position(1) = 0.449;
Lgnd.Position(2) = 0;
legend("Step size","Rejections",'NumColumns',2,'FontSize',16,'Interpreter','latex');


%% Timing

Stiff = zeros(5,3);
Stiff(1,1) = rs10;
Stiff(2,1) = function_calls10.nFun;
Stiff(3,1) = function_calls10.nJac;
Stiff(4,1) = function_calls10.nLU;
Stiff(5,1) = function_calls10.nStep;

Stiff(1,2) = rs11;
Stiff(2,2) = function_calls11.nFun;
Stiff(3,2) = function_calls11.nJac;
Stiff(4,2) = function_calls11.nLU;
Stiff(5,2) = function_calls11.nStep;

Stiff(1,3) = rs12;
Stiff(2,3) = function_calls12.nFun;
Stiff(3,3) = function_calls12.nJac;
Stiff(4,3) = function_calls12.nLU;
Stiff(5,3) = function_calls12.nStep;

NonStiff = zeros(5,3);

NonStiff(1,1) = rs7;
NonStiff(2,1) = function_calls7.nFun;
NonStiff(3,1) = function_calls7.nJac;
NonStiff(4,1) = function_calls7.nLU;
NonStiff(5,1) = function_calls7.nStep;

NonStiff(1,2) = rs8;
NonStiff(2,2) = function_calls8.nFun;
NonStiff(3,2) = function_calls8.nJac;
NonStiff(4,2) = function_calls8.nLU;
NonStiff(5,2) = function_calls8.nStep;

NonStiff(1,3) = rs9;
NonStiff(2,3) = function_calls9.nFun;
NonStiff(3,3) = function_calls9.nJac;
NonStiff(4,3) = function_calls9.nLU;
NonStiff(5,3) = function_calls9.nStep;




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