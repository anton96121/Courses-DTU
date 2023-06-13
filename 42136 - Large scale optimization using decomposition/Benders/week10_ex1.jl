# Benders template (no rays)

#-------------------------------------------------------------------
# Intro definitions
using JuMP
using CPLEX
#-------------------------------------------------------------------
c = 5
f = -3
A = [1; 2; 1]
B = [2; -1; -3]
b = [4; 0; -13]

#-------------------------------------------------------------------
# Master problem
mas=Model(CPLEX.Optimizer)

# Variables
@variable(mas, q )
@variable(mas, 10>=y>=0 , Int)

@objective(mas, Min, f*y + q)

function solve_master(ustar)
    # Add Constraints
    @constraint(mas, (b-B*y)'*ustar <= q)

    optimize!(mas)
    
    return objective_value(mas)
    
end
#-------------------------------------------------------------------


#-------------------------------------------------------------------
# Sub problem
function solve_sub( ystar )
    sub=Model(CPLEX.Optimizer)

    # Variables
    @variable(sub, u[i = 1:3]>=0 , base_name = "dual_var")

    # Objective
    @objective(sub, Max,  (b-B*ystar)'*u)

    # Constraints
    @constraint(sub, A'*u<=c )

    optimize!(sub)

    return (objective_value(sub),  value.(sub[:u]))
end    
#-------------------------------------------------------------------




#-------------------------------------------------------------------
# main code
let
    UB=Inf
    LB=-Inf
    Delta=0.1
    ybar= 0
    it=1
    while (UB-LB>Delta)
        
        (sub_obj, ustar )=solve_sub(ybar)
        
        UB=min(UB,sub_obj + f*ybar)
        
        mas_obj=solve_master( ustar )
        
        ybar=value.(y)
        LB=mas_obj
        
        println("It: $(it) UB: $(UB) LB: $(LB)  Sub: $(sub_obj)")
        it+=1
    end
end
println("Correct Ending")
