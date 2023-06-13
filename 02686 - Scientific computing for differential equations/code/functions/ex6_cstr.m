clear all

%% CSTR 1D
method = "DOPRI54";
options = struct('step_control', true, 'initialStepSize', true,'control_type', "PID");
optionsOde = odeset('RelTol',1e-6,'AbsTol',1e-6);

min = 60; %[s]
F = [0.7/min,0.6/min,0.5/min,0.4/min,0.3/min,0.2/min,0.3/min,0.4/min,0.5/min,0.6/min,0.7/min,0.7/min,0.2/min,0.2/min,0.7/min,0.7/min]; %[L/s]
Tin = 273.65; % [K]
V = .105; % [L]
beta = -(-560)/(1*4.186);%[mL/min]
Ca_in = 1.6/2; %[mol/L]
Cb_in = 2.4/2; %[mol/L]

t0 = 0; 
tend = 120; 
h = 1;

X = [];
Xode45 = [];
Xode15s = [];
T = [];
Fplot = [];
x0 = Tin;
x0ode45 = Tin;
x0ode15s = Tin;

for i=1:length(F)
    param = [F(i) Tin V beta Ca_in Cb_in];
    [X1,T1,function_calls1,hs1] = ODEsolver(@CSTR1D,param,h,t0,tend,x0, method,options);
   
    [~,X2]=ode45(@CSTR1D,T1,x0ode45,optionsOde, param);
    [~,X3]=ode15s(@CSTR1D,T1,x0ode15s,optionsOde,param);

    X = [X; X1];
    Xode45 = [Xode45; X2];
    Xode15s = [Xode15s; X3];

    T = [T T1];
    x0 = X1(end);
    x0ode45 = X2(end);
    x0ode15s = X3(end);

    t0 = T(end); 
    tend = T(end)+120; 
    Fplot = [Fplot, repelem(F(i), length(T1))];
end



fig = figure('Position', [10 10 1000 800]);
tend = 120;

subplot(2, 1, 1)
plot(T/60,X-273.15, 'LineWidth', 1.5)
hold on 
plot(T/60,Xode45(:,1)-273.15, '-.','LineWidth', 1.5)
plot(T/60,Xode15s(:,1)-273.15, '--','LineWidth', 1.5)
ylabel('T [$^\circ$C]', 'FontSize', 16, 'Interpreter','latex')
title('1D CSTR',  'FontSize', 20, 'Interpreter','latex')
grid on
legend(method,'ODE45','ODE15s','fontsize',12, 'Interpreter','latex')
hold off

subplot(2, 1, 2)
%plot([0 repelem((1:(length(F)-1))*tend,2) length(F)*tend]/60, repelem(F,2)*60*1000, 'LineWidth', 1.5)
plot(T/60, Fplot*60*1000, 'LineWidth', 1.5)
xlabel('t [min]', 'FontSize',16, 'Interpreter','latex')
ylabel('F [mL/min]', 'FontSize',16, 'Interpreter','latex')
ylim([0,1000])
grid on 

%% CSTR 3D
method = "DOPRI54";
Atol = 1e-10;
Rtol = 1e-10;
options = struct('step_control', true, 'initialStepSize', false, 'Rtol', Rtol, 'Atol', Atol,'control_type', "I", 'hmax', 6);
optionsOde = odeset('RelTol',Rtol,'AbsTol',Atol);

min = 60; %[s]
F = [0.7/min,0.6/min,0.5/min,0.4/min,0.3/min,0.2/min,0.3/min,0.4/min,0.5/min,0.6/min,0.7/min,0.7/min,0.2/min,0.2/min,0.7/min,0.7/min]; %[L/s]
Tin = 273.65; % [K]
V = .105; % [L]
beta = -(-560)/(1*4.186);%[mL/min]
Ca_in = 1.6/2; %[mol/L]
Cb_in = 2.4/2; %[mol/L]

t0 = 0; 
tend = 120; 
N = 30;
h = 1;

X = [];
Ca = [];
Cb = [];
Xode45 = [];
Caode45 = [];
Cbode45 = [];
Xode15s = [];
Caode15s = [];
Cbode15s = [];

Fplot = [];
T = [];
x0 = [Ca_in;Cb_in;Tin];
x0ode45 = x0;
x0ode15s = x0;

h_array = [];
for i=1:length(F)
    param = [F(i) Tin V beta Ca_in Cb_in];
    [X1,T1,function_calls1,hs1] = ODEsolver(@CSTR3D,param,h,t0,tend,x0, method,options);

    [~,X2]=ode45(@CSTR3D,T1,x0ode45,optionsOde, param);
    [~,X3]=ode15s(@CSTR3D,T1,x0ode15s,optionsOde,param);

    Ca = [Ca; X1(:,1)];
    Cb = [Cb; X1(:,2)];
    X = [X; X1(:,3)];
    
    Caode45 = [Caode45; X2(:,1)];
    Cbode45 = [Cbode45; X2(:,2)];
    Xode45 = [Xode45; X2(:,3)];

    Caode15s = [Caode15s; X3(:,1)];
    Cbode15s = [Cbode15s; X3(:,2)];
    Xode15s = [Xode15s; X3(:,3)];

    T = [T T1];
    h_array = [h_array hs1(:,1) hs1];
    x0 = X1(end,:);
    x0ode45 = X2(end,:);
    x0ode15s = X3(end,:);
    
    t0 = T(end); 
    tend = T(end)+120; 
    Fplot = [Fplot, repelem(F(i), length(T1))];
    h = hs1(1,end);
end


fig = figure('Position', [10 10 1000 800]);
subplot(2, 1, 1)
plot(T/60,X-273.15, 'LineWidth', 1.5)
hold on 
plot(T/60,Xode45(:,1)-273.15, '-.','LineWidth', 1.5)
plot(T/60,Xode15s(:,1)-273.15, '--','LineWidth', 1.5)
ylabel('T [$^\circ$C]', 'FontSize', 16, 'Interpreter','latex')
title('3D CSTR',  'FontSize', 20, 'Interpreter','latex')
grid on
legend(method,'ODE45','ODE15s','fontsize',12, 'Interpreter','latex')
hold off

subplot(2, 1, 2)
plot(T/60, Fplot*60*1000, 'LineWidth', 1.5)
xlabel('t [min]', 'FontSize',16, 'Interpreter','latex')
ylabel('F [mL/min]', 'FontSize',16, 'Interpreter','latex')
ylim([0,1000])
grid on 


fig = figure('Position', [10 10 1000 800]);
subplot(2, 1, 1)
plot(T/60,Ca, 'LineWidth', 1.5)
hold on
plot(T/60,Caode45, '-.','LineWidth', 1.5)
plot(T/60,Caode15s, '-.','LineWidth', 1.5)
ylabel('$C_A$ [$mol/L$]', 'FontSize', 16, 'Interpreter','latex')
title('3D CSTR',  'FontSize',16, 'Interpreter','latex')
legend(method,'ODE45','ODE15s','fontsize',12, 'Interpreter','latex')
grid on

subplot(2, 1, 2)
plot(T/60,Cb, 'LineWidth', 1.5)
hold on
plot(T/60,Cbode45, '-.','LineWidth', 1.5)
plot(T/60,Cbode15s, '-.','LineWidth', 1.5)
ylabel('$C_B$ [$mol/L$]', 'FontSize', 16, 'Interpreter','latex')
xlabel('t [min]', 'FontSize',16, 'Interpreter','latex')
legend(method,'ODE45','ODE15s','fontsize',12, 'Interpreter','latex')
grid on 


fig = figure;
T1 = T/60;
hs1 = h_array;
semilogy(T1,hs1(1,:),'LineWidth', 1.5)
hold on
semilogy(T1(hs1(2,:)>1),hs1(1,hs1(2,:)>1),'x', 'color', 'r','MarkerSize',12)
title('Used step sizes','FontSize',16,'Interpreter','latex')
xlabel("t",'FontSize',16,'Interpreter','latex')
ylabel("h",'FontSize',16,'Interpreter','latex')


 

% add a bit space to the figure
fig = gcf;
%fig.Position(3) = fig.Position(3) + 250;
% add legend
Lgnd = legend('show');
Lgnd.Position(1) = 0.449;
Lgnd.Position(2) = 0;
legend("Step size","Rejections",'NumColumns',2,'FontSize',16,'Interpreter','latex');