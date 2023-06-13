using JuMP
using GLPK

#ex1
function one_d_minkowski_weyl()
    m = Model(with_optimizer(GLPK.Optimizer))
    @variable(m, l1 >= 0)
    @variable(m, l2 >= 0)
    @variable(m, l3 >= 0)
    @variable(m, l4 >= 0)
    @objective(m, Max, (2*l2+6*l3)+3*(3*l3+6*l4))
    @constraint(m, l1 + l2 + l3 + l4 == 1)

    optimize!(m)

    if termination_status(m) == MOI.OPTIMAL
        println("Objective value: ", JuMP.objective_value(m))
        println("l1 = ", JuMP.value(l1))
        println("l2 = ", JuMP.value(l2))
        println("l3 = ", JuMP.value(l3))
        println("l4 = ", JuMP.value(l4))
    else
        println("Failed optimization. Return code: ", termination_status(m) )
    end
end

function one_d_normal()
    m = Model(with_optimizer(GLPK.Optimizer))
    @variable(m, x1 >= 0)
    @variable(m, x2 >= 0)
    @objective(m, Max, x1+3*x2)
    @constraint(m, (x1-2)/4 - x2/3 <= 0)
    @constraint(m, (x1/2) - 6 + x2 <= 0)

    optimize!(m)

    if termination_status(m) == MOI.OPTIMAL
        println("Objective value: ", JuMP.objective_value(m))
        println("x1 = ", JuMP.value(x1))
        println("x2 = ", JuMP.value(x2))
    else
        println("Failed optimization. Return code: ", termination_status(m) )
    end
end


#ex2
function two()
    m = Model(with_optimizer(GLPK.Optimizer))
    @variable(m, l[i=1:3] >= 0)
    x1 = l[1]+l[2]+2*l[3]
    x2 = l[1]+2*l[2]+2*l[3]
    @objective(m, Max, x1 + 2*x2)
    @constraint(m, 3*x1 + x2 <= 13)
    @constraint(m, -x1 -3*x2 <= -7)
    @constraint(m, sum(l) == 1)

    optimize!(m)

    if termination_status(m) == MOI.OPTIMAL
        println("Objective value: ", JuMP.objective_value(m))
        println("l = ", JuMP.value.(l))
    else
        println("Failed optimization. Return code: ", termination_status(m) )
    end
end

#ex3
function threeBin()
    m = Model(with_optimizer(GLPK.Optimizer))
    @variable(m, x[i=1:4], Bin)
    @objective(m, Max, x[1]+x[2]+x[3]+x[4])
    @constraint(m, 3*x[1]-2*x[2]+2*x[3]-x[4] >= 3)
    @constraint(m, x[1]+x[2]+4*x[3]+x[4] <= 3)

    optimize!(m)

    if termination_status(m) == MOI.OPTIMAL
        println("Objective value: ", JuMP.objective_value(m))
        println("x = ", JuMP.value.(x))
    else
        println("Failed optimization. Return code: ", termination_status(m) )
    end
end

function threeLin()
    m = Model(with_optimizer(GLPK.Optimizer))
    @variable(m, 0<=x[i=1:4]<=1)
    @objective(m, Max, x[1]+x[2]+x[3]+x[4])
    @constraint(m, 3*x[1]-2*x[2]+2*x[3]-x[4] >= 3)
    @constraint(m, x[1]+x[2]+4*x[3]+x[4] <= 3)

    optimize!(m)

    if termination_status(m) == MOI.OPTIMAL
        println("Objective value: ", JuMP.objective_value(m))
        println("x = ", JuMP.value.(x))
    else
        println("Failed optimization. Return code: ", termination_status(m) )
    end
end

function threeDWc1()
    m = Model(with_optimizer(GLPK.Optimizer))
    @variable(m, l[i=1:4] >= 0)
    x = [
        l[1]+l[2]+l[3]+l[4],
        l[3],
        l[1]+l[2]+l[3],
        l[2]
        ]
    @objective(m, Max, x[1]+x[2]+x[3]+x[4])
    @constraint(m, x[1]+x[2]+4*x[3]+x[4] <= 3)
    @constraint(m, sum(l) == 1)

    optimize!(m)

    if termination_status(m) == MOI.OPTIMAL
        println("Objective value: ", JuMP.objective_value(m))
        println("l = ", JuMP.value.(l))
    else
        println("Failed optimization. Return code: ", termination_status(m) )
    end
end

function threeDWc2()
    m = Model(with_optimizer(GLPK.Optimizer))
    @variable(m, 0 <= l[i=1:8])
    x = [l[2]+l[4]+l[6]+l[8], l[3]+l[4]+l[7]+l[8], 0, l[5]+l[6]+l[7]+l[8]]
    @objective(m, Max, sum(x))
    @constraint(m, 3*x[1] - 2*x[2] + 2*x[3] -x[4] >= 3)
    @constraint(m, sum(l) == 1)

    optimize!(m)

    if termination_status(m) == MOI.OPTIMAL
        println("Objective value: ", JuMP.objective_value(m))
        println("l = ", JuMP.value.(l))
    else
        println("Failed optimization. Return code: ", termination_status(m) )
    end
end