using JuMP
using GLPK
using Base
using LinearAlgebra

#### ex1, 1
N = [1 1 1 0 0 0 0 0 0; 
     0 0 0 1 1 1 0 0 0;
     0 0 0 0 0 0 1 1 1]
p = [2,12,6,1,13,7,1,10,4]
w = [1,6,3,1,8,4,1,7,3]
Q = 10


function multipleChoiceKnapsack(p,w,Q,N)

    m = Model(with_optimizer(GLPK.Optimizer))
    @variable(m, x[i=1:length(p)], Bin)

    @objective(m, Max, p'*x)

    @constraint(m,N*x .==ones(size(N)[1],1))
    
    @constraint(m, w'*x <= Q)

    optimize!(m)

    if termination_status(m) == MOI.OPTIMAL
        println("Objective value: ", JuMP.objective_value(m))
        println("x = ", JuMP.value.(x))
    else
        println("Failed optimization. Return code: ", termination_status(m) )
    end
end

multipleChoiceKnapsack(p,w,Q,N)

#### ex1, 2
function multipleChoiceKnapsackLP(p,w,Q,N)

    m = Model(with_optimizer(GLPK.Optimizer))
    # We just remove ',Bin' to get it in LP mode
    @variable(m, x[i=1:length(p)])

    @objective(m, Max, p'*x)

    @constraint(m,N*x .==ones(size(N)[1],1))

    @constraint(m, w'*x <= Q)

    # We add this constraint due to the lp relaxation
    @constraint(m, 0 .<= x .<= 1)

    optimize!(m)

    if termination_status(m) == MOI.OPTIMAL
        println("Objective value: ", JuMP.objective_value(m))
        println("x = ", JuMP.value.(x))
    else
        println("Failed optimization. Return code: ", termination_status(m) )
    end
end

multipleChoiceKnapsackLP(p,w,Q,N)


#### ex1, 3
# First we find all permutations
k = 9
permuts = [ [ bit == '1' ? 1 : 0 for bit in Base.bin(UInt16(n),k,false) ] for n in 0:2^k-1 ]

# Then we find all valid permutations
X̄¹  = permuts[[w'*permuts[i]<=Q for i in 1:2^k]];
X̄¹ = mapreduce(permutedims, vcat, X̄¹)';
start_point = findall(x-> (x == 1) ,[col == [1,0,0,1,0,0,1,0,0] for col in eachcol(X̄¹)]);

function RMP(p,w,Q,N, X̄¹, verbose)

    m = Model(with_optimizer(GLPK.Optimizer))
    @variable(m, λ[i=1:size(X̄¹)[2]])

    @objective(m, Max, p'*(X̄¹*λ))
    @constraint(m, A0,N*(X̄¹*λ) .==ones(size(N)[1],1))
    @constraint(m, conv, sum(λ) == 1)
    @constraint(m, λ.>= 0)

    optimize!(m)
    μ = -JuMP.dual.(A0)
    κ = -JuMP.dual.(conv)
    if verbose
        if termination_status(m) == MOI.OPTIMAL
            println("Objective value: ", JuMP.objective_value(m))
            println("λ = ", JuMP.value.(λ))
            println("μ = ", -JuMP.dual.(A0))
            println("κ = ", -JuMP.dual.(conv))
        else
            println("Failed optimization. Return code: ", termination_status(m) )
        end
    end

    return(JuMP.value.(λ),μ,κ)
end

function SubProb(p,w,Q,N, μ, κ, verbose)

    m = Model(with_optimizer(GLPK.Optimizer))
    @variable(m, x[i=1:length(p)], Bin)

    @objective(m, Max, p'*x-dot(μ'*N,x)-κ)
    @constraint(m,w'*x<= Q)

    optimize!(m)

    out = JuMP.value.(x)

    if verbose
        if termination_status(m) == MOI.OPTIMAL
            println("Objective value: ", JuMP.objective_value(m))
            println("x = ", JuMP.value.(x))
        else
            println("Failed optimization. Return code: ", termination_status(m) )
        end
    end     
    return(out,JuMP.objective_value(m))
end

flag = true
verbose = false
set_in = copy(start_point)
i = 1
while (flag)
    λ, μ, κ =  RMP(p,w,Q,N, X̄¹[:,set_in],verbose);
    answer = copy(X̄¹[:,set_in]*λ);
    println([round(x,digits = 3) for x in X̄¹[:,set_in]*λ])
    x_add, val  = SubProb(p,w,Q,N, μ, κ,verbose);
    flag = val>0.0

    add_point = findall(x-> (x == 1) ,[col == x_add for col in eachcol(X̄¹)]);
    if flag
        append!( set_in, add_point );
    end
end
[abs(round(x,digits = 3)) for x in answer]


# ex2
# First we find all permutations
k = 9
permuts = [ [ bit == '1' ? 1 : 0 for bit in Base.bin(UInt16(n),k,false) ] for n in 0:2^k-1 ]

# Then we find all valid permutations
X̄¹  = permuts[[all(N*permuts[i] .==ones(size(N)[1],1)) for i in 1:2^k]];
X̄¹ = mapreduce(permutedims, vcat, X̄¹)';
start_point = findall(x-> (x == 1) ,[col == [1,0,0,1,0,0,1,0,0] for col in eachcol(X̄¹)]);

function RMP2(p,w,Q,N, X̄¹, verbose)

    m = Model(with_optimizer(GLPK.Optimizer))
    @variable(m, λ[i=1:size(X̄¹)[2]])

    @objective(m, Max, p'*(X̄¹*λ))
    @constraint(m, A0,w'*(X̄¹*λ) <= Q)
    @constraint(m, conv, sum(λ) == 1)
    @constraint(m, λ.>= 0)

    optimize!(m)
    μ = -JuMP.dual.(A0)
    κ = -JuMP.dual.(conv)
    if verbose
        if termination_status(m) == MOI.OPTIMAL
            println("Objective value: ", JuMP.objective_value(m))
            println("λ = ", JuMP.value.(λ))
            println("μ = ", -JuMP.dual.(A0))
            println("κ = ", -JuMP.dual.(conv))
        else
            println("Failed optimization. Return code: ", termination_status(m) )
        end
    end

    return(JuMP.value.(λ),μ,κ)
end

function SubProb2(p,w,Q,N, μ, κ, verbose)

    m = Model(with_optimizer(GLPK.Optimizer))
    @variable(m, x[i=1:length(p)], Bin)

    @objective(m, Max, p'*x-dot(μ'*w',x)-κ)
    @constraint(m,N*x.== ones(size(N)[1],1))

    optimize!(m)

    out = JuMP.value.(x)

    if verbose
        if termination_status(m) == MOI.OPTIMAL
            println("Objective value: ", JuMP.objective_value(m))
            println("x = ", JuMP.value.(x))
        else
            println("Failed optimization. Return code: ", termination_status(m) )
        end
    end     
    return(out,JuMP.objective_value(m))
end

flag = true
verbose = false
set_in = copy(start_point)
i = 1
while (flag)
    λ, μ, κ =  RMP2(p,w,Q,N, X̄¹[:,set_in],verbose);
    answer = copy(X̄¹[:,set_in]*λ);
    println([round(x,digits = 3) for x in X̄¹[:,set_in]*λ])
    x_add, val  = SubProb2(p,w,Q,N, μ, κ,verbose);
    flag = val>0.0

    add_point = findall(x-> (x == 1) ,[col == x_add for col in eachcol(X̄¹)]);
    if flag
        append!( set_in, add_point );
    end
end
[abs(round(x,digits = 3)) for x in answer]