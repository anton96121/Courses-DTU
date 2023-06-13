using JuMP
using CPLEX
using GLPK

#m = Model(CPLEX.Optimizer)
m = Model(GLPK.Optimizer)

numNodes = 13

#node probability - nodes are number left to right and top to bottom within a stage
pi = [1 0.25 0.50 0.25 0.15 0.075 0.025 0.10 0.30 0.10 0.075 0.075 0.10]

#predecessor node
predNode = [-1 1 1 1 2 2 2 3 3 3 4 4 4]

numCustomers = 5
numFactories = 3

#node demand at each scenario tree node for each customer
demand =[0 150 160 170 150 160 170 150 160 170 150 160 170;
         0 100 120 135 100 120 135 100 120 135 100 120 135;
         0 250 270 300 250 270 300 250 270 300 250 270 300;
         0 300 325 350 300 325 350 300 325 350 300 325 350;
         0 600 700 800 600 700 800 600 700 800 600 700 800]

#Transportation cost per unit
cost = [2.49 5.21 3.76 4.85 2.07;
        1.46 2.54 1.83 1.86 4.76;
        3.26 3.08 2.60 3.76 4.45]

#capacity
capacity = [500 450 650]

#variables
@variable(m, x[i=1:numFactories, j=1:numCustomers, l=1:numNodes]>=0)
@variable(m, b[j=1:numCustomers, l=1:numNodes]>=0)
@variable(m, i[j=1:numCustomers, l=1:numNodes]>=0)

#backorder and inventory cost by node
backCost = [0 5 5 5 24 24 24 24 24 24 24 24 24]
invCost  = [0 1 1 1 12 12 12 12 12 12 12 12 12]

#objective
@objective(m, Min, sum(pi[l]*cost[i,j]*x[i,j,l] for i=1:numFactories for j=1:numCustomers for l=1:numNodes)
                  +sum(pi[l]*backCost[l]*b[j,l] for j=1:numCustomers for l=1:numNodes)
                  +sum(pi[l]*invCost[l]*i[j,l]  for j=1:numCustomers for l=1:numNodes))

##capcity and demand constraints
@constraint(m, [i=1:numFactories, l=1:numNodes], sum(x[i,j,l] for j=1:numCustomers) <= capacity[i])
@constraint(m, [j=1:numCustomers, l=2:numNodes],  i[j, predNode[l]] + b[j,l] + sum(x[i,j,predNode[l]] for i=1:numFactories) == demand[j,l] + b[j,predNode[l]]+ i[j,l])

#initial conditions - no inventory no backorders
@constraint(m, [j=1:numCustomers], i[j,1]==0)
@constraint(m, [j=1:numCustomers], b[j,1]==0)

optimize!(m)
println("")
println("Optimal Objective Value: ",objective_value(m))
