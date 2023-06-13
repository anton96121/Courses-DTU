using JuMP
using CPLEX
using GLPK

#m = Model(CPLEX.Optimizer)
m=Model(GLPK.Optimizer)

pi = [0.25 0.50 0.25;
      0.25 0.50 0.25;
      0.25 0.50 0.25;
      0.30 0.40 0.30;
      0.30 0.40 0.30]

demandScenario =[150 160 170;
                 100 120 135;
                 250 270 300;
                 300 325 350;
                 600 700 800]

capacity = [500 450 650]

cost = [2.49 5.21 3.76 4.85 2.07;
        1.46 2.54 1.83 1.86 4.76;
        3.26 3.08 2.60 3.76 4.45]

numCustomers = 5
numFactories = 3
numDemands = 3

demand = zeros(numCustomers)

for i=1:numCustomers
    sum = 0.0
    for w=1:numDemands
        sum=sum+pi[i,w]*demandScenario[i,w]
    end
    demand[i]=sum
end

println("New Demand Values")
println(demand)

@variable(m, x[i=1:numFactories,j=1:numCustomers]>=0)
@variable(m, z[i=1:numFactories]>=0)
@variable(m, p[j=1:numCustomers]>=0)
@variable(m, w[j=1:numCustomers]>=0)

prodCost = 14
salePrice = 24
wasteCost = 4

@objective(m, Max, sum(salePrice*p[j] for j=1:numCustomers)
                   -sum(cost[i,j]*x[i,j] for i=1:numFactories for j=1:numCustomers)
                   -sum(wasteCost*w[j] for j=1:numCustomers)
                   -sum(prodCost*z[i] for i=1:numFactories))

@constraint(m, [i=1:numFactories], z[i] - sum(x[i,j] for j=1:numCustomers) == 0)
@constraint(m, [i=1:numFactories], z[i] <= capacity[i])
@constraint(m, [j=1:numCustomers],  sum(x[i,j] for i=1:numFactories)-p[j]-w[j]==0)
@constraint(m, [j=1:numCustomers], p[j]<=demand[j])

optimize!(m)

println("Optimal Objective Value: ", objective_value(m))
