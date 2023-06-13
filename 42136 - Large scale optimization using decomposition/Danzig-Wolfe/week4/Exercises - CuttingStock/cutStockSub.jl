module test

using JuMP, GLPK, LinearAlgebra

W=62 
w = [14 20 27]
pi = [0.25 0.3333333333333333 0.5]

n = length(w)
sub = Model(GLPK.Optimizer)

@variable(sub, x[1:n] >= 0, Int )

@objective(sub, Min, 1-sum(pi[i] * x[i] for i=1:n) )
@constraint(sub, sum(w[i] * x[i] for i=1:n) <= W)

optimize!(sub)

if termination_status(sub) == MOI.OPTIMAL
    println("Objective value: ", JuMP.objective_value(sub))
    println("x = ", JuMP.value.(x))
else
    println("Optimize was not succesful. Return code: ", termination_status(m))
end

end
