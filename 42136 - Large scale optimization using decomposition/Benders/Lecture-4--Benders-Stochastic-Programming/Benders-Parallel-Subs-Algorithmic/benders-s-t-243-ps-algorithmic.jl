using JuMP
using CPLEX
using GLPK

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

# determine all possible scenarios
function scenarioGen(n::Int64, j::Int64, scen::Array{Int64,1}, result::Array{Array{Int64,1},1})
      for w=1:n
            push!(scen, w)
            if j< 5
                  scenarioGen(n, j+1, scen, result)
            else
                  push!(result,deepcopy(scen))
            end
            pop!(scen)
      end
end

scen=Int64[]
result = Array{Int64,1}[]

scenarioGen(size(demandScenario,2),1,scen, result)
fullDemandScenarios = Array{Int64,1}[[] for i=1:length(result)]

# generate scenario probabilities
prob = zeros(length(result),1)
for i=1:length(result)
      sum = 1.0
      for j=1:length(result[i])
            idx = result[i][j]
            sum = sum*pi[j,idx]
            push!(fullDemandScenarios[i],demandScenario[j,idx])
      end
      prob[i]=sum

end

capacity = [500 450 650]

cost = [2.49 5.21 3.76 4.85 2.07;
        1.46 2.54 1.83 1.86 4.76;
        3.26 3.08 2.60 3.76 4.45]

numScenarios = length(prob)

prodCost = 14
salePrice = 24
wasteCost = 4

#-------------------------------------------------------------------
# Master problem
#masterproblem=Model(CPLEX.Optimizer)
#set_optimizer_attribute(masterproblem, "CPXPARAM_ScreenOutput", 0)
masterproblem=Model(GLPK.Optimizer)

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
                sum(sum(sum(alphabar[w][j]*(x[i,j]) for i=1:numFactories) for j=1:numCustomers)
                +sum(betabar[w][j]*fullDemandScenarios[w][j] for j=1:numCustomers) for w=1:numScenarios) >= q)

    optimize!(masterproblem)
    return objective_value(masterproblem)
end
#-------------------------------------------------------------------

#-------------------------------------------------------------------
# Sub problem
function solve_subproblem(xbar, w)

    alpha = zeros(5)
    beta = zeros(5)

    for j=1:numCustomers
        sum = 0.0
        for i=1:numFactories
            sum += xbar[i,j]
        end
        if(sum > fullDemandScenarios[w][j])
            alpha[j] = -wasteCost*prob[w]
            beta[j] = salePrice*prob[w]-alpha[j]
        else
            alpha[j] = salePrice*prob[w]
            beta[j] = 0
        end
    end
    #println(alpha)
    #println(beta)


    subobj=(sum(alpha[j]*sum(xbar[i,j] for i=1:numFactories) for j=1:numCustomers)
                                +sum(beta[j]*fullDemandScenarios[w][j] for j=1:numCustomers))

    return (subobj, alpha, beta)

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
    sub_obj=zeros(Float64,numScenarios)
    alpha = Array{Float64,1}[[] for i=1:numScenarios]
    #alpha=zeros(Float64,3)
    beta = Array{Float64,1}[[] for i=1:numScenarios]
    it=1
    while (UB-LB>Delta)
        for w=1:numScenarios
            (sub_obj[w],alpha[w],beta[w]) = solve_subproblem(xbar, w)
        end

        mas_component = sum(-cost[i,j]*xbar[i,j] for i=1:numFactories for j=1:numCustomers) -sum(prodCost*zbar[i] for i=1:numFactories)
        LB=max(LB, sum(sub_obj[w] for w=1:numScenarios) + mas_component)
        mas_obj=solve_master(alpha, beta)
        xbar = value.(x)
        zbar = value.(z)
        UB=mas_obj
        println("It: $(it) UB: $(UB) LB: $(LB)")#" #Sub: $(sub_obj)")
        it+=1
    end
end
println("Optimal Solution Found")
