module test

using JuMP, GLPK, LinearAlgebra

b = [3;2;3] 
# our patterns
Xprime = [ 4 0 0; 0 3 0; 0 0 2; 3 1 0; 1 1 1 ]'

(m,n) = size(Xprime)

master = Model(GLPK.Optimizer)

@variable(master, lambda[1:n] >= 0, Int )
@objective(master, Min, sum(lambda[j] for j=1:n) )
@constraint(master, cons[i=1:m], sum(Xprime[i,j]*lambda[j] for j=1:n) >= b[i] )

optimize!(master)

if termination_status(master) == MOI.OPTIMAL
    println("")
    println("Objective value: ", JuMP.objective_value(master))
    println("lambda = ", value.(lambda))
    #println("pi = ", dual.(cons))
else
    println("Optimize was not succesful. Return code: ", termination_status(master))
end

end
