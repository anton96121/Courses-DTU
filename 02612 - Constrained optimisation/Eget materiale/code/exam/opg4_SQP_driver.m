close all
%%
%{
This is the driver for exercise 4. 
This file contains:
    - A contour plot of the Himmelblau problem (Exercise 4.4)
    - A comparison of fmincon, CasADi and our own SQP-BFGS (Exercise 4.5)
    - A test of the efficiency of subsolvers
    - A test of the correctness of subsolvers
    - A test of SQP BFGS 
    - A test of SQP BFGS Line Search
    - A test of non-monotone strategy
    - A test of SQP BFGS Line Search with infeasiblity handling
    - A test of SQP BFGS Trust Region
    - A test of SQP BFGS Trust Region from infeasible start

%}

%% Exercise 4.4, Himmelblau's test problem plotted

fig = figure('Position', [100, 0, 1000,1000]);
hold on
edges = 7;
%Plot the contours themselves
himmelblauContours(edges)

% Plot all the points of interest (max, min, saddle)
h1 = plot(-3.073025751,-0.08135304429,'d','MarkerFaceColor','r', ... 
      'MarkerEdgeColor','r','markersize',10);
h2 = plot(0.08667750456,2.884254701,'d','MarkerFaceColor','r', ... 
      'MarkerEdgeColor','r','markersize',10);
h3 = plot(3,2,'s','MarkerFaceColor','g', ... 
      'MarkerEdgeColor','r','markersize',15);
h4 = plot(-0.2983476136,2.895620844,'s','MarkerFaceColor','g', ... 
      'MarkerEdgeColor','r','markersize',15);
h5 = plot(-1.424243078,0.3314960331,'h','MarkerFaceColor','y', ... 
      'MarkerEdgeColor','r','markersize',15);
h6 = plot(-3.654605171,2.737718273,'s','MarkerFaceColor','g', ... 
      'MarkerEdgeColor','r','markersize',15);
h7 = plot(-3.548537388,-1.419414955,'s','MarkerFaceColor','g', ... 
      'MarkerEdgeColor','r','markersize',15);
h8 = plot(-0.4869360343,-0.1947744137,'h','MarkerFaceColor','y', ... 
      'MarkerEdgeColor','r','markersize',15);
h9 = plot(3.216440661,1.286576264,'h','MarkerFaceColor','y', ... 
      'MarkerEdgeColor','r','markersize',15);
 legend([h8,h7,h2],{'Maximum', 'Minimum', 'Saddle Point'})

% show constraints as dark areas over contour
x1lim = edges; x1 = -x1lim:0.01:x1lim;
x2lim = edges; x2 = -x1lim:0.01:x2lim;
x1f = -edges:0.01:edges;
h1 = fill(x1f, x1f.^2 + 4.*x1f + 4, [0 0 0], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h2 = patch('Faces',[1 2 3], 'Vertices', [5.348 x1lim; 4.324 -7; 100 -100], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h3 = patch('Faces',[1 2 3], 'Vertices', [-x1lim 0.4*-x1lim; x1lim 0.4*x1lim; x1lim -999], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h4 = patch('Faces',[1 2 3], 'Vertices', [-x1lim 0.4*-x1lim+7; x1lim 0.4*x1lim + 7; x1lim +999], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
set( get( get( h4, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h3, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );

h5 = patch('Faces',[1 2 3 4], 'Vertices', [-edges -5; edges -5; edges -edges; -edges -edges], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h6 = patch('Faces',[1 2 3 4], 'Vertices', [-edges edges; -5 edges; -5 -edges; -edges -edges], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h7 = patch('Faces',[1 2 3 4], 'Vertices', [-edges edges; edges edges; edges 5; -edges 5], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h8 = patch('Faces',[1 2 3 4], 'Vertices', [5 edges; edges edges; edges -edges; 5 -edges], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
set( get( get( h5, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h6, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h7, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h8, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
% set axis as needed
axis([-edges edges -edges edges])
title('Himmelblaus Test Problem')

%% Exercise 4.5, Comparison of fmincon, CasADi and our own SQP-BFGS
%Solve Himmelblau using fmincon, CasADi, and or own solver.
sympref('FloatingPointOutput',1);
l = [-5;-5];
u = [5;5];
cl = [0;0];
cu = [47;70];
addpath('C:\Users\anton\\OneDrive - Københavns Universitet\Uni\Uni\8. semester\Constrained optimization\Casadi')
import casadi.*

%Solve with fmincon
options = optimset('Display', 'off');
xfmin = fmincon(@objfminconHimmel, [0 0], [], [], [], [], l, u, @consfminconHimmel, options);

%Solve with our solver
options = struct('log',false, 'infesibility_handling', false, 'method', 'SQP', 'subsolver', 'own solver');
x0 = [0;0];
[xown,~] = SQPSolver(x0,@objHimmel,@consHimmel,l,u,cl,cu,options);


x1 = SX.sym('x1');
x2 = SX.sym('x2');
nlp = struct ('x',[x1;x2], 'f',(x1^2+x2-11)^2+(x1+x2^2-7)^2,'g',[(x1+2)^2-x2; -4*x1+10*x2]);

%Solve with CasADi
options = struct;
options.ipopt.print_level = 0;
options.print_time = 0;
S = nlpsol('S', 'ipopt', nlp,options);
r = S('x0', [0,0], 'lbg',0,'ubg',inf);
x_Cas =full(r.x);

%Print the solutions
fprintf('CasADi, solution: [%f,%f] \n', x_Cas(1), x_Cas(2)); 
fprintf('fmincon, solution: [%f,%f] \n', xfmin(1), xfmin(2)); 
fprintf('Own solver, solution: [%f,%f] \n', xown(1), xown(2)); 
 


%% Efficiency of subsolvers
%Compare quadprog and our internal solver stepwith growing size random EQP
% problems.
options = optimset('Display', 'off');
times = zeros(20, 2);
n_large = 50;

%Solve with both our and quadprog for several sizes
for i=1:20
        n = n_large*i;
        [H, g, A, b] = generateRandomEQP(n,n/2);
        C = A;
        dl = b-3;
        du = b+3;
        l = zeros(n,1);
        u = ones(n,1);
        x0 = zeros(n,1);
        
        start = cputime;
        [x,z,~,~] = intSQP(H,g,C,l,u,dl,du,x0);
        times(i,1) = cputime-start;
        
        start = cputime;
        [x,~,~,~,z] = quadprog(H,g,-[C'; -C'],-[dl;-du],[],[],l,u,x0, options);
        times(i,2) = cputime-start;
        
        disp("Iteration " +i +"/" + 20);
end


%Plot subsolver time spent on subproblem
figure
hold on
h1 = plot(n_large:n_large:20*n_large,times(:,1),'r');
h2 = plot(n_large:n_large:20*n_large,times(:,2),'b');
legend([h1,h2],{'Own solver', 'quadprog'})
title('Comparison of underlying subsolvers')
ylabel('time[s]')
xlabel('n')

l = [-5;-5];
u = [5;5];
cl = [0;0];
cu = [47;70];

%% Test of the correctness of subsolvers
%Compare our subsolver to quadprog subsolver, compare - they take the same
%path
close all
sympref('FloatingPointOutput',1);
l = [-5;-5];
u = [5;5];
cl = [0;0];
cu = [47;70];
% Comparison of own solver and the quadprog
x0 = [0.0;0.0];
options = struct('log',true, 'infesibility_handling', false, 'method', 'SQP', 'subsolver', 'own solver');
[x,~, Hist] = SQPSolver(x0,@objHimmel,@consHimmel,l,u,cl,cu,options);
xHist = Hist.xHist;
xHistLatex_1 = latex(sym(round(xHist,3)));
time_own = Hist.timePerformence;
%Trace using our solver
himmelPlot(x0,x,xHist, 'SQP', 'own solver','b',5,'')


options = struct('log',true, 'infesibility_handling', false, 'method', 'SQP', 'subsolver', 'quadprog');
[x,~, Hist] = SQPSolver(x0,@objHimmel,@consHimmel,l,u,cl,cu,options);
xHist = Hist.xHist;
xHistLatex_2 = latex(sym(round(xHist,3)));
time_quad = Hist.timePerformence;
%Trace using quadprog
himmelTrace(x0,x,xHist, 'r', 'quadprog','--')

legend show


%% Test of SQP BFGS 
%Test and plot four traces of SQP BFGS all with starting point near each
% their own minimum of the Himmelblau problem
close all
sympref('FloatingPointOutput',1);
l = [-5;-5];
u = [5;5];
cl = [0;0];
cu = [47;70];

% Several traces
options = struct('log',true, 'infesibility_handling', true, 'method', 'SQP', 'subsolver', 'own solver');
x0 = [0;0];
[x,~, Hist] = SQPSolver(x0,@objHimmel,@consHimmel,l,u,cl,cu,options);
Hist_1 = Hist;
xHist = Hist.xHist;
xHistLatex_1 = latex(sym(round(xHist,3)));
himmelPlot(x0,x,xHist, 'SQP', 'Trace 1')


x0 = [0;3.8416];
[x,~, Hist] = SQPSolver(x0,@objHimmel,@consHimmel,l,u,cl,cu,options);
Hist_2 = Hist;
xHist = Hist.xHist;
xHistLatex_2 = latex(sym(round(xHist,3)));
himmelTrace(x0,x,xHist, 'r', 'Trace 2')


x0 = [-4;4];
[x,~, Hist] = SQPSolver(x0,@objHimmel,@consHimmel,l,u,cl,cu,options);
Hist_3 = Hist;
xHist = Hist.xHist;
xHistLatex_3 = latex(sym(round(xHist,3)));
himmelTrace(x0,x,xHist, 'k', 'Trace 3')

x0 = [ -4.5;-1];
start = cputime;
[x,~, Hist] = SQPSolver(x0,@objHimmel,@consHimmel,l,u,cl,cu,options);
endtime_4 = cputime-start;
functionCalls_4 = Hist.functionCalls;
Hist_4 = Hist;
xHist = Hist.xHist;
xHistLatex_4 = latex(sym(round(xHist,3)));
himmelTrace(x0,x,xHist, 'g', 'Trace 4')

legend show

%% Test of SQP BFGS Line Search
%Test and plot four traces of line search converging to each their own
%minimum.
close all
sympref('FloatingPointOutput',1);
l = [-5;-5];
u = [5;5];
cl = [0;0];
cu = [47;70];
% Several traces
options = struct('log',true, 'infesibility_handling', false, 'method', 'SQP_ls', 'subsolver', 'own solver');
x0 = [0;0];
[x,~, Hist] = SQPSolver(x0,@objHimmel,@consHimmel,l,u,cl,cu,options);
Hist_1 = Hist;
xHist = Hist.xHist;
xHistLatex_1 = latex(sym(round(xHist,3)));
stepHist_1 = latex(sym(round(Hist.stepLength,3)));
himmelPlot(x0,x,xHist, 'SQP linesearch', 'Trace 1')


x0 = [0;3.8416];
[x,~, Hist] = SQPSolver(x0,@objHimmel,@consHimmel,l,u,cl,cu,options);
Hist_2 = Hist;
xHist = Hist.xHist;
xHistLatex_2 = latex(sym(round(xHist,3)));
stepHist_2 = latex(sym(round(Hist.stepLength,3)));
himmelTrace(x0,x,xHist, 'r', 'Trace 2')


x0 = [-4;4];
[x,~, Hist] = SQPSolver(x0,@objHimmel,@consHimmel,l,u,cl,cu,options);
Hist_3 = Hist;
xHist = Hist.xHist;
xHistLatex_3 = latex(sym(round(xHist,3)));
stepHist_3 = latex(sym(round(Hist.stepLength,3)));
himmelTrace(x0,x,xHist, 'k', 'Trace 3')

x0 = [ -4.5;-1];
[x,~, Hist] = SQPSolver(x0,@objHimmel,@consHimmel,l,u,cl,cu,options);
Hist_4 = Hist;
xHist = Hist.xHist;
xHistLatex_4 = latex(sym(round(xHist,3)));
stepHist_4 = latex(sym(round(Hist.stepLength,3)));
himmelTrace(x0,x,xHist, 'g', 'Trace 4')

legend show

%% Test of non-monotone strategy
%Compare non-monotone and monotone strategy for line search
close all
sympref('FloatingPointOutput',1);
l = [-5;-5];
u = [5;5];
cl = [0;0];
cu = [47;70];

%Non monotone line search
x0 = [ -4.5;-1];
options = struct('log',true, 'infesibility_handling', false, 'method', 'SQP_ls', 'subsolver', 'own solver', 'non_monotone', true);
[x,~, Hist] = SQPSolver(x0,@objHimmel,@consHimmel,l,u,cl,cu,options);
xHist = Hist.xHist;
Hist_1 = Hist;
steps1 = Hist.stepLength;
xHistLatex_1 = latex(sym(round(xHist,3)));
stepHist_1 = latex(sym(round(Hist.stepLength,3)));
himmelPlot(x0,x,xHist, 'SQP linesearch', 'With non monotone strategy')

%Monotone line search
options = struct('log',true, 'infesibility_handling', false, 'method', 'SQP_ls', 'subsolver', 'own solver', 'non_monotone', false);
[x,~, Hist] = SQPSolver(x0,@objHimmel,@consHimmel,l,u,cl,cu,options);
xHist = Hist.xHist;
Hist_2 = Hist;
steps2 = Hist.stepLength;
xHistLatex_2 = latex(sym(round(xHist,3)));
stepHist_2 = latex(sym(round(Hist.stepLength,3)));
himmelTrace(x0,x,xHist, 'r', 'Without non monotone strategy')


legend show


%% Test of SQP BFGS Line Search with infeasiblity handling
%Testing line search with infeasibility handling with just one trace.
close all
sympref('FloatingPointOutput',1);
l = [-5;-5];
u = [5;5];
cl = [0;0];
cu = [47;70];
% When infeasiblity handling is turned on, quadprog is the only solver for
%  the sub problem which is avalible.

x0 = [-9;6];
try
    options = struct('log',true, 'infesibility_handling', false, 'method', 'SQP_ls');
    [x,z, Hist] = SQPSolver(x0,@objHimmel,@consHimmel,l,u,cl,cu,options);
catch
    disp('We see that the subproblem is infeasible due to an infeasible linearization');
end
options = struct('log',true, 'infesibility_handling', true, 'method', 'SQP_ls');
[x,z, Hist] = SQPSolver(x0,@objHimmel,@consHimmel,l,u,cl,cu,options);
xHist = Hist.xHist;
xHistLatex_1 = latex(sym(round(xHist,3)));
stepHist_1 = latex(sym(round(Hist.stepLength,3)));
himmelPlot(x0,x,xHist, 'SQP linesearch with infeasiblity handling', 'Trace 1','r',9)
legend show

%% SQP BFGS Trust Region
%Four traces using Trust Region converging to each their own minimum
close all
sympref('FloatingPointOutput',1);
l = [-5;-5];
u = [5;5];
cl = [0;0];
cu = [47;70];
% When one uses the SQP based on a trust region, quadprog is the only solver for
%  the sub problem which is avalible.
options = struct('log',true, 'method', 'SQP_trust', 'trust_region',0.5);
x0 = [0;0];
[x,~, Hist] = SQPSolver(x0,@objHimmel,@consHimmel,l,u,cl,cu,options);
xHist = Hist.xHist;
Hist_1 = Hist;
xHistLatex_1 = latex(sym(round(xHist,3)));
himmelPlot(x0,x,xHist, 'SQP trust region', 'Trace 1')


x0 = [0;3.8416];
[x,~, Hist] = SQPSolver(x0,@objHimmel,@consHimmel,l,u,cl,cu,options);
xHist = Hist.xHist;
Hist_2 = Hist;
xHistLatex_2 = latex(sym(round(xHist,3)));
himmelTrace(x0,x,xHist, 'r', 'Trace 2')


x0 = [-4;4];
[x,z, Hist] = SQPSolver(x0,@objHimmel,@consHimmel,l,u,cl,cu,options);
xHist = Hist.xHist;
Hist_3 = Hist;
xHistLatex_3 = latex(sym(round(xHist,3)));
himmelTrace(x0,x,xHist, 'k', 'Trace 3')

x0 = [ -4.5;-1];
[x,~, Hist] = SQPSolver(x0,@objHimmel,@consHimmel,l,u,cl,cu,options);
xHist = Hist.xHist;
Hist_4 = Hist;
xHistLatex_4 = latex(sym(round(xHist,3)));
himmelTrace(x0,x,xHist, 'g', 'Trace 4')

legend show

%% Test of SQP BFGS Trust Region from infeasible start
%Testing just one trace from infeasible start using trust region and
%infeasibility handling
close all
sympref('FloatingPointOutput',1);
l = [-5;-5];
u = [5;5];
cl = [0;0];
cu = [47;70];
% When infeasiblity handling is turned on, quadprog is the only solver for
%  the sub problem which is avalible.

x0 = [-2.805118; 3.131312]; 
options = struct('log',true, 'method', 'SQP_trust', 'trust_region',5, 'penalty',1000);
[x,z, Hist] = SQPSolver(x0,@objHimmel,@consHimmel,l,u,cl,cu,options);
xHist = Hist.xHist;
xHistLatex_1 = latex(sym(round(xHist,3)));
himmelPlot(x0,x,xHist, 'SQP trust region with infeasible start', 'Trace 1','r',9)
legend show


%% objective, constraints, derivitives and plotting functions
function [f,df] = objHimmel(x)

x1 = x(1);
x2 = x(2);

temp1 = x1^2+x2-11;
temp2 = x1+x2^2-7;

f = temp1^2+temp2^2;

df = zeros(2,1);
df(1) = 4*x1*temp1+2*temp2;
df(2) = 2*temp1+4*x2*temp2;

end

function [c,dc] = consHimmel(x,full,d)

if nargin<2
    full = false;
    d = 0;
end

if nargin<3 && full
    error('One also needs to input d to det the full dc')
end  


x1 = x(1);
x2 = x(2);

c = zeros(2,1);
c(1) = (x1+2)^2-x2;
c(2) = -4*x1+10*x2;




dc = zeros(2,2);
dc(1,1) = 2*x1+4;
dc(1,2) = -1;
dc(2,1) = -4;
dc(2,2) = 10;

dc = dc';

if full
    c = [x; -x; c; -c];     
    c = c-d;

    dc = [eye(2) -eye(2) dc -dc];
end

end

function [f] = objfminconHimmel(x)

x1 = x(1);
x2 = x(2);

temp1 = x1^2+x2-11;
temp2 = x1+x2^2-7;

f = temp1^2+temp2^2;

df = zeros(2,1);
df(1) = 4*x1*temp1+2*temp2;
df(2) = 2*temp1+4*x2*temp2;

end

function [c, ceq] = consfminconHimmel(x)

x1 = x(1);
x2 = x(2);

c = zeros(2,1);
c(1) = (x1+2)^2-x2;
c(2) = -4*x1+10*x2;
c(1) = -c(1);
c(2) = - c(2);
ceq = [0];

%dc = zeros(2,2);
%dc(1,1) = 2*x1+4;
%dc(1,2) = -1;
%dc(2,1) = -4;
%dc(2,2) = 10;

%dc = dc';

end


function [] = himmelTrace(x0,x,xHist, tracecolor, legendname,linepattern)
if nargin<6
    linepattern = '-.';
end
hold on 
h0 = plot(x0(1,:),x0(2,:),'k.','markersize',40);
h1 = plot(x(1,:),x(2,:),'r.','markersize',40);
plot(xHist(1,:),xHist(2,:),[tracecolor,linepattern],'linewidth',2,'DisplayName',legendname);
set( get( get( h0, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
%legend(h2, legendname)
hold off
end


function [] = himmelPlot(x0,x,xHist, titlename,legendname, color, edges ,linepattern)
if nargin <6
    color = 'b';
end
if nargin <7
    edges = 5;
end
if nargin<8
    linepattern = '-.';
end
fig = figure('Position', [100, 0, 1000,1000]);
hold on
himmelblauContours(edges)
h1 = plot(x0(1,:),x0(2,:),'k.','markersize',40);
h2 = plot(x(1,:),x(2,:),'r.','markersize',40);
plot(xHist(1,:),xHist(2,:),strcat(color,linepattern),'linewidth',2,'DisplayName',legendname);
%legend(h1, legendname)
set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
% show constraints
x1lim = edges; x1 = -x1lim:0.01:x1lim;
x2lim = edges; x2 = -x1lim:0.01:x2lim;
x1f = -edges:0.01:edges;
h3 = fill(x1f, x1f.^2 + 4.*x1f + 4, [0 0 0], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
set( get( get( h3, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
h4 = patch('Faces',[1 2 3], 'Vertices', [-x1lim 0.4*-x1lim; x1lim 0.4*x1lim; x1lim -999], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
set( get( get( h4, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );

h5 = patch('Faces',[1 2 3 4], 'Vertices', [-edges -5; edges -5; edges -edges; -edges -edges], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h6 = patch('Faces',[1 2 3 4], 'Vertices', [-edges edges; -5 edges; -5 -edges; -edges -edges], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h7 = patch('Faces',[1 2 3 4], 'Vertices', [-edges edges; edges edges; edges 5; -edges 5], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h8 = patch('Faces',[1 2 3 4], 'Vertices', [5 edges; edges edges; edges -edges; 5 -edges], 'FaceColor', 'black', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
set( get( get( h5, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h6, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h7, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
set( get( get( h8, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
% set axis as needed
axis([-edges edges -edges edges])
title(titlename)
hold off
end

