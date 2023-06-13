using JuMP
using CPLEX
using GLPK

#m = Model(CPLEX.Optimizer)
m = Model(GLPK.Optimizer)

#probability of each scenario - scenarios are numbered top to bottom
pi = [0.15 0.075 0.025 0.10 0.30 0.10 0.075 0.075 0.10 ]

numScenarios = 9
numTimeSteps = 2
numCustomers = 5
numFactories = 3

demand =[150 160 170;
         100 120 135;
         250 270 300;
         300 325 350;
         600 700 800]

# generate demand scenarios
demandScenario = Array{Array{Float64,1},2}(undef,numScenarios, numTimeSteps)

for i=1:9
   for j=1:2
       demandScenario[i,j]=[]
    end
end

count = 1;
for i=1:numFactories
   for j=1:3
      demandScenario[count,1]=demand[:,i]
      demandScenario[count,2]=demand[:,j]
      global count+=1
   end
end



#capacity = [500 450 650]
capacity = [500 450 650]

#Transportation cost per unit
cost = [2.49 5.21 3.76 4.85 2.07;
        1.46 2.54 1.83 1.86 4.76;
        3.26 3.08 2.60 3.76 4.45]

#variables
@variable(m, x[i=1:numFactories, j=1:numCustomers, t=1:numTimeSteps, w=1:numScenarios]>=0)
@variable(m, b[j=1:numCustomers, t=0:numTimeSteps, w=1:numScenarios]>=0)
@variable(m, i[j=1:numCustomers, t=0:numTimeSteps, w=1:numScenarios]>=0)

#backorder and inventory cost by time step
backCost = [5 24]
invCost  = [1 12]

@objective(m, Min,  sum(pi[w]*cost[i,j]*x[i,j,t,w] for i=1:numFactories for j=1:numCustomers for w=1:numScenarios for t=1:numTimeSteps)+
                    sum(pi[w]*invCost[t]*i[j,t,w] for j=1:numCustomers for w=1:numScenarios for t=1:numTimeSteps)+
                    sum(pi[w]*backCost[t]*b[j,t,w] for j=1:numCustomers for w=1:numScenarios for t=1:numTimeSteps))

@constraint(m, [i=1:numFactories, t=1:numTimeSteps, w=1:numScenarios], sum(x[i,j,t, w] for j=1:numCustomers) <= capacity[i])
@constraint(m, [j=1:numCustomers, t=1:numTimeSteps, w=1:numScenarios],  i[j, t-1, w] + b[j, t, w] + sum(x[i,j,t, w] for i=1:numFactories) == demandScenario[w,t][j] + b[j,t-1,w]+ i[j,t,w])

@constraint(m, [j=1:numCustomers, w=1:numScenarios], i[j,0,w]==0)
@constraint(m, [j=1:numCustomers, w=1:numScenarios], b[j,0,w]==0)

# non-anticipativity constraints
@constraint(m, [i=1:numFactories, j=1:numCustomers, w=2:numScenarios], x[i,j,1,w] ==  x[i,j,1,w-1] )
@constraint(m, [i=1:numFactories, j=1:numCustomers, w=2:3], x[i,j,2,w] ==  x[i,j,2,w-1] )
@constraint(m, [j=1:numCustomers, w=2:3], b[j,1,w] == b[j,1,w-1] )
@constraint(m, [j=1:numCustomers, w=2:3], i[j,1,w] == i[j,1,w-1] )
@constraint(m, [i=1:numFactories, j=1:numCustomers, w=5:6], x[i,j,2,w] ==  x[i,j,2,w-1] )
@constraint(m, [j=1:numCustomers, w=5:6], b[j,1,w] == b[j,1,w-1])
@constraint(m, [j=1:numCustomers, w=5:6], i[j,1,w] == i[j,1,w-1] )
@constraint(m, [i=1:numFactories, j=1:numCustomers, w=8:numScenarios], x[i,j,2,w] ==  x[i,j,2,w-1] )
@constraint(m, [j=1:numCustomers, w=8:numScenarios], b[j,1,w] == b[j,1,w-1])
@constraint(m, [j=1:numCustomers, w=8:numScenarios], i[j,1,w] == i[j,1,w-1] )

optimize!(m)

println("")
println("Optimal objective: ", objective_value(m))
