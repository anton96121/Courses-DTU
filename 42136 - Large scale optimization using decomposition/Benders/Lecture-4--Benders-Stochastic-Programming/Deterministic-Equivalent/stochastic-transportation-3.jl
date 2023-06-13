using JuMP
using CPLEX
using GLPK

#m = Model(CPLEX.Optimizer)
#set_parameter(m, "CPXPARAM_Benders_Strategy",3)
m=Model(GLPK.Optimizer)

pi = [0.25 0.50 0.25]

demandScenario =[150 160 170;
                 100 120 135;
                 250 270 300;
                 300 325 350;
                 600 700 800]

capacity = [500 450 650]

cost = [2.49 5.21 3.76 4.85 2.07;
        1.46 2.54 1.83 1.86 4.76;
        3.26 3.08 2.60 3.76 4.45]

xvalues =[0.0 0.0 0.0 0.0 500.0;
        125.0 0.0 0.0 325.0 0.0;
        35.0 118.75 272.5 0.0 200.0]

numScenarios = 3
numCustomers = 5
numFactories = 3

@variable(m, x[i=1:numFactories,j=1:numCustomers]>=0)
@variable(m, z[i=1:numFactories]>=0)
@variable(m, p[j=1:numCustomers, w=1:numScenarios]>=0)
@variable(m, wa[j=1:numCustomers, w=1:numScenarios]>=0)

prodCost = 14
salePrice = 24
wasteCost = 4

@objective(m, Max, sum(salePrice*pi[w]*p[j,w] for w=1:numScenarios for j=1:numCustomers)
                   -sum(cost[i,j]*x[i,j] for i=1:numFactories for j=1:numCustomers)
                   -sum(wasteCost*pi[w]*wa[j,w] for w=1:numScenarios, j=1:numCustomers)
                   -sum(prodCost*z[i] for i=1:numFactories))

@constraint(m, [i=1:numFactories], z[i] - sum(x[i,j] for j=1:numCustomers) == 0)
@constraint(m, [i=1:numFactories], z[i] <= capacity[i])
@constraint(m, [j=1:numCustomers, w=1:numScenarios],  sum(x[i,j] for i=1:numFactories)-p[j,w]-wa[j,w]==0)
@constraint(m, [j=1:numCustomers, w=1:numScenarios], p[j,w]<=demandScenario[j,w])

optimize!(m)

println("")
println("Objective value: ",objective_value(m))
