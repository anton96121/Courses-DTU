using JuMP
using CPLEX
####################################################################################

############################         EXERCISE 1         ############################

####################################################################################
# Data
demand=[ 12, 14, 16, 18, 20, 22, 24, 26, 28, 30] # Demand of newspapers in each scenario
S=length(demand)
prob=[0.05, 0.10, 0.10, 0.10, 0.15, 0.15, 0.10, 0.10, 0.10, 0.05 ] # probability of scenario
c=20 # purchase price
p=70 # selling price
h=10 # scrap value

# Deterministic-equivalent
direct = Model(CPLEX.Optimizer)
@variable(direct, y >= 0, Int) # bought newspapers
@variable(direct, x[1:S] >= 0) # sold newspapers
# maximize profit from sales and scrap minus cost
@objective(direct, Max, sum( prob[s]*((p-h)*x[s] + (h-c)*y) for s=1:S))
# cant sell more than the demand
@constraint(direct, [s=1:S], x[s] <= demand[s])
# cant sell more newspapers than purchased
@constraint(direct, [s=1:S], x[s] <= y)
optimize!(direct)
println("Obj. Val.: ", objective_value(direct))
println("x = ", value.(x), ", y = ", value(y))

####################################################################################

############################         EXERCISE 2         ############################

####################################################################################


# Data
demand=[ 12, 14, 16, 18, 20, 22, 24, 26, 28, 30] # Demand of newspapers in each scenario
S=length(demand)
prob=[0.05, 0.10, 0.10, 0.10, 0.15, 0.15, 0.10, 0.10, 0.10, 0.05 ] # probability of scenario
c=20 # purchase price
p=70 # selling price
h=10 # scrap value

# direct problem
direct = Model(GLPK.Optimizer)
ybar=21
@variable(direct, x[1:S] >= 0) # sold newspapers

# maximize profit from sales and scrap minus cost
@objective(direct, Max, sum( prob[s]*((p-h)*x[s] + (h-c)*ybar) for s=1:S))
# cant sell more than the demand
@constraint(direct, dem[s=1:S], x[s] <= demand[s])
# cant sell more newspapers than purchased
@constraint(direct, bought[s=1:S], x[s] <= ybar)
optimize!(direct)
println("Obj. Val.: ", objective_value(direct))
println("x = ", value.(x))
println("dem dual = ", dual.(dem))
println("bought dual = ",dual.(bought))


####################################################################################

############################         EXERCISE 4         ############################

####################################################################################

#-------------------------------------------------------------------
#************************************************************************
# Data
demand=[ 12, 14, 16, 18, 20, 22, 24, 26, 28, 30] # Demand of newspapers in each scenario
S=length(demand)
prob=[0.05, 0.10, 0.10, 0.10, 0.15, 0.15, 0.10, 0.10, 0.10, 0.05 ] # probability of scenario
c=20 # purchase price
p=70 # selling price
h=10 # scrap value
#************************************************************************
#-------------------------------------------------------------------
# Master problem
mas=Model(CPLEX.Optimizer)
# Variables
@variable(mas, q )
@variable(mas, 0 <= y, Int)
@objective(mas, Max, (h-c)*y + q)

function solve_master( alphabar, betabar )
    # Add Constraints
    @constraint(mas, sum( demand[s]*alphabar[s] for s=1:S) + sum( y*betabar[s] for s=1:S) >= q)
    optimize!(mas)
    return objective_value(mas)
end
#-------------------------------------------------------------------
#-------------------------------------------------------------------
# Sub problem
function solve_sub( ybar )
    sub = Model(CPLEX.Optimizer)
    @variable(sub, alpha[1:S] >= 0)
    @variable(sub, beta[1:S] >= 0)
    @objective(sub, Min, sum( demand[s]*alpha[s] for s=1:S) + sum( ybar*beta[s] for s=1:S) )
    @constraint(sub, [s=1:S], alpha[s] + beta[s] >= prob[s]*(p-h) )
    solution = optimize!(sub)
    return (objective_value(sub), value.(alpha), value.(beta) )
end
#-------------------------------------------------------------------
#-------------------------------------------------------------------
# main code
let
    UB=Inf
    LB=-Inf
    Delta=0
    ybar=0
    it=1
    while (UB-LB>Delta)
        (sub_obj, alpha, beta)=solve_sub(ybar)

        LB=max(LB,sub_obj + (h-c)*ybar )
        mas_obj=solve_master(alpha,beta)
        ybar=value(y)
        UB=mas_obj
        println("It: $(it) UB: $(UB) LB: $(LB) Sub: $(sub_obj)")
        println("Beta: $(beta) Alpha: $(alpha)")
        it+=1
    end
end
println("Correct Ending")

####################################################################################

############################         EXERCISE 5         ############################

####################################################################################

# Data
demand=[ 12, 14, 16, 18, 20, 22, 24, 26, 28, 30] # Demand of newspapers in each scenario
S=length(demand)
prob=[0.05, 0.10, 0.10, 0.10, 0.15, 0.15, 0.10, 0.10, 0.10, 0.05 ] # probability of scenario
c=20 # purchase price
p=70 # selling price
h=10 # scrap value
#************************************************************************
#-------------------------------------------------------------------
# Master problem
mas=Model(CPLEX.Optimizer)
# Variables
@variable(mas, q )
@variable(mas, 0 <= y <= 30, Int)
@objective(mas, Max, (h-c)*y + q)
function solve_master( alphabar, betabar )
    # Add Constraints
    @constraint(mas, sum( demand[s]*alphabar[s] for s=1:S) + sum( y*betabar[s] for s=1:S) >= q)
    optimize!(mas)
    return objective_value(mas)
end
#-------------------------------------------------------------------
#-------------------------------------------------------------------
# Sub problem
function solve_sub_scenario( ybar, s)
    sub = Model(CPLEX.Optimizer)
    @variable(sub, alpha >= 0)
    @variable(sub, beta >= 0)
    @objective(sub, Min, demand[s]*alpha + ybar*beta )
    @constraint(sub, alpha + beta >= prob[s]*(p-h) )
    solution = optimize!(sub)
    return (objective_value(sub), value(alpha), value(beta) )
end
#-------------------------------------------------------------------
#-------------------------------------------------------------------
# main code
let
    UB=Inf
    LB=-Inf
    Delta=0
    ybar=0
    sub_obj=zeros(Float64,S)
    alpha=zeros(Float64,S)
    beta=zeros(Float64,S)
    it=1
    while (UB-LB>Delta)
        for s=1:S
            (sub_obj[s], alpha[s], beta[s])=solve_sub_scenario(ybar,s)
        end
        LB=max(LB,sum(sub_obj[s] for s=1:S) + (h-c)*ybar )
        mas_obj=solve_master(alpha,beta)
        ybar=value(y)
        UB=mas_obj
        println("It: $(it) UB: $(UB) LB: $(LB) Sub: $(sub_obj)")
        it+=1
    end
end
println("Correct Ending")