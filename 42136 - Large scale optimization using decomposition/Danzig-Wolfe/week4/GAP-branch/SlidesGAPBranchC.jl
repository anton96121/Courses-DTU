include("ColGenGeneric.jl")

using JuMP, GLPK

function GAPMIPbranch(w, p, c)
    (m,n) = size(w)
    myModel = Model(GLPK.Optimizer)
    @variable(myModel, 0 <= x[1:m,1:n] <= 1, Int)
    @objective(myModel, Max, sum(p[i,j]*x[i,j] for i=1:m for j=1:n))
    # each job must be served
    @constraint(myModel, [j=1:n],sum(x[i,j] for i=1:m) == 1)
    @constraint(myModel, capacity[i=1:m],sum(w[i,j]*x[i,j] for j=1:n) <= c[i])
    @constraint(myModel, branch1, x[1,5] <= 0)

    # We have to construct the blocks in a slightly more complicated way since they now are going to contain constraints of differerent types (both <= and >=)
    # We need to make sure that we store "ConstraintRef"s which is the general reference to a constraint.
    blocks = [ConstraintRef[] for i=1:m]
    for i=1:m
        push!(blocks[i], capacity[i])
    end
    # we add the branching constraint to the subproblem belonging to machine 1
    push!(blocks[1], branch1)
    
    return myModel, blocks
end

function run()
    # w(i,j) = capacity used when assigning job j to machine i
    # p(i,j) = profit of assigning job j to machine i
    w = [
    8	1	9	4	10	10	9	1	8	6	4	5	10	3	4
    2	1	8	10	6	1	10	8	8	7	3	10	2	1	9
    5	5	4	8	10	4	6	6	2	8	8	2	1	9	2
    3	7	3	8	8	2	4	4	7	1	1	6	10	6	6    ]
    p = [
    3	6	8	2	5	1	2	4	4	2	8	2	4	2	6
    7	2	4	6	3	10	4	5	7	10	2	6	3	4	1
    3	9	9	5	4	6	10	9	3	4	4	5	2	7	1
    9	4	9	6	6	9	7	8	6	8	10	10	5	5	4
    ]
    c = [23 23 23 23]

    gapModel, blocks = GAPMIPbranch(w, p, c)
    DW_ColGen.DWColGenEasy(gapModel, blocks)
end

run()
