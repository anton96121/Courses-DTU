using JuMP
using GLPK
using CPLEX

# Data
pi = [0.25 0.5 0.25]

demandScenario =[150 160 170;
                 100 120 135;
                 250 270 300;
                 300 325 350;
                 600 700 800]

capacity = [500 450 650]

cost = [2.49 5.21 3.76 4.85 2.07;
        1.46 2.54 1.83 1.86 4.76;
        3.26 3.08 2.60 3.76 4.45]

prodCost = 14
salePrice = 24
wasteCost = 4

numFactories = 3
numCustomers = 5
numScenarios = 3
#-------------------------------------------------------------------

masterproblem=Model(CPLEX.Optimizer)
#set_optimizer_attribute(masterproblem, "CPXPARAM_ScreenOutput", 0)
#masterproblem=Model(GLPK.Optimizer)

# Variables - single cut subproblem
@variable(masterproblem, q )

# x-variables prod variables
@variable(masterproblem, x[i=1:numFactories,j=1:numCustomers]>=0)
@variable(masterproblem, z[i=1:numFactories]>=0)

#master constraints
@constraint(masterproblem, [i=1:numFactories], z[i] - sum(x[i,j] for j=1:numCustomers) == 0)
@constraint(masterproblem, [i=1:numFactories], z[i] <= capacity[i])

@objective(masterproblem, Max, -sum(cost[i,j]*x[i,j] for i=1:numFactories for j=1:numCustomers) -sum(prodCost*z[i] for i=1:numFactories)+q)

# solve master problem
function solve_master(alphabar, betabar)
    # Add Constraints
    @constraint(masterproblem,
                sum(alphabar[j,w]*sum(x[i,j] for i=1:numFactories) for j=1:numCustomers for w=1:numScenarios)
                +sum(betabar[j,w]*demandScenario[j,w] for j=1:numCustomers for w=1:numScenarios) >= q)

    optimize!(masterproblem)
    return objective_value(masterproblem)
end
#-------------------------------------------------------------------

#-------------------------------------------------------------------
# Sub problem
function solve_subproblem(xbar)

    subproblem=Model(CPLEX.Optimizer)
    #set_optimizer_attribute(subproblem, "CPXPARAM_ScreenOutput", 0)
    #subproblem=Model(GLPK.Optimizer)

    @variable(subproblem, alpha[j=1:numCustomers, w=1:numScenarios])
    @variable(subproblem, beta[j=1:numCustomers, w=1:numScenarios]>=0)

    @objective(subproblem, Min, sum(alpha[j,w]*sum(xbar[i,j] for i=1:numFactories) for j=1:numCustomers for w=1:numScenarios)
                                +sum(beta[j,w]*demandScenario[j,w] for j=1:numCustomers for w=1:numScenarios))

    @constraint(subproblem, [j=1:numCustomers, w=1:numScenarios], +alpha[j,w]+beta[j,w] >= salePrice*pi[w])
    @constraint(subproblem, [j=1:numCustomers, w=1:numScenarios], alpha[j,w] >= -wasteCost*pi[w])

    optimize!(subproblem)

    return (objective_value(subproblem), value.(alpha), value.(beta))

end
#-------------------------------------------------------------------

#-------------------------------------------------------------------
# main code
let
    UB=Inf
    LB=-Inf
    Delta=1e-6
    #initial solution - set all shipment values zero
    xbar=zeros(3,5)
    #initial solution - set all production variable values to zero
    zbar=zeros(Float64,3)
    sub_obj=zeros(Float64,3)
    alpha=zeros(Float64,3)
    beta=zeros(Float64,3)
    it=1
    while (UB-LB>Delta)
        (sub_obj,alpha,beta) = solve_subproblem(xbar)
        mas_component = sum(-cost[i,j]*xbar[i,j] for i=1:numFactories for j=1:numCustomers) -sum(prodCost*zbar[i] for i=1:numFactories)
        LB=max(LB,sub_obj + mas_component)
        mas_obj=solve_master(alpha, beta)
        xbar = value.(x)
        zbar = value.(z)
        UB=mas_obj
        println("It: $(it) UB: $(UB) LB: $(LB)  Sub: $(sub_obj)")
        it+=1
    end
end
println("Optimal Solution Found")