cd("C:\\Users\\anton\\OneDrive - Københavns Universitet\\Uni\\Uni\\10. semester\\decomposition\\Danzig-Wolfe\\week4\\Exercises - GAP-branch")
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
    # We have to construct the blocks in a slightly more complicated way since they now are
    # going to contain constraints of differerent types (both <= and >= once we start branching)
    # We need to make sure that we store "ConstraintRef"s which is the general reference to a constraint.
    blocks = [ConstraintRef[] for i=1:m]
    for i=1:m
        push!(blocks[i], capacity[i])
    end

    return myModel, blocks
end

function run()
    # w(i,j) = capacity used when assigning job j to machine i
    # p(i,j) = profit of assigning job j to machine i
    w = [
    8 6 1 7 7 7 5 5 7 7 3 3 8 5 4
    5 5 3 8 5 6 5 9 9 6 1 9 5 6 3
    3 2 4 4 9 1 7 3 3 3 5 3 7 7 1
    6 7 9 9 3 5 2 1 5 5 4 4 6 2 1]
    p = [
    9 9 1 4 6 4 2 3 5 1 9 9 7 6 3
    9 1 8 6 7 8 2 6 6 5 3 4 7 5 3
    4 8 1 8 1 4 1 4 2 6 2 1 2 5 1
    7 9 7 8 3 6 2 3 4 7 9 9 3 9 1]
    c = [22 18 18 19]

    gapModel, blocks = GAPMIPbranch(w, p, c)
    DW_ColGen.DWColGenEasy(gapModel, blocks)
end

run()

# We find the solution is fractional with obj_val = 102.5. We choose x[3,10] to branch on

# This is the x[3,10]>=1 branch
function GAPMIPbranch(w, p, c)
    (m,n) = size(w)
    myModel = Model(GLPK.Optimizer)
    @variable(myModel, 0 <= x[1:m,1:n] <= 1, Int)
    @objective(myModel, Max, sum(p[i,j]*x[i,j] for i=1:m for j=1:n))
    # each job must be served
    @constraint(myModel, [j=1:n],sum(x[i,j] for i=1:m) == 1)
    @constraint(myModel, capacity[i=1:m],sum(w[i,j]*x[i,j] for j=1:n) <= c[i])
    @constraint(myModel, branch1, x[3,10] >= 1)
    # We have to construct the blocks in a slightly more complicated way since they now are going to contain constraints of differerent types (both <= and >=)
    # We need to make sure that we store "ConstraintRef"s which is the general reference to a constraint. 
    blocks = [ConstraintRef[] for i=1:m]
    for i=1:m
        push!(blocks[i], capacity[i])
    end
    # we add the branching constraint to the subproblem belonging to machine 1
    push!(blocks[3], branch1)

    return myModel, blocks
end

function run()
    # w(i,j) = capacity used when assigning job j to machine i
    # p(i,j) = profit of assigning job j to machine i
    w = [
    8 6 1 7 7 7 5 5 7 7 3 3 8 5 4
    5 5 3 8 5 6 5 9 9 6 1 9 5 6 3
    3 2 4 4 9 1 7 3 3 3 5 3 7 7 1
    6 7 9 9 3 5 2 1 5 5 4 4 6 2 1]
    p = [
    9 9 1 4 6 4 2 3 5 1 9 9 7 6 3
    9 1 8 6 7 8 2 6 6 5 3 4 7 5 3
    4 8 1 8 1 4 1 4 2 6 2 1 2 5 1
    7 9 7 8 3 6 2 3 4 7 9 9 3 9 1]
    c = [22 18 18 19]

    gapModel, blocks = GAPMIPbranch(w, p, c)
    DW_ColGen.DWColGenEasy(gapModel, blocks)
end

run()
# We find the solution is fractional with obj_val = 102.25.

# This is the x[3,10]<=0 branch
function GAPMIPbranch(w, p, c)
    (m,n) = size(w)
    myModel = Model(GLPK.Optimizer)
    @variable(myModel, 0 <= x[1:m,1:n] <= 1, Int)
    @objective(myModel, Max, sum(p[i,j]*x[i,j] for i=1:m for j=1:n))
    # each job must be served
    @constraint(myModel, [j=1:n],sum(x[i,j] for i=1:m) == 1)
    @constraint(myModel, capacity[i=1:m],sum(w[i,j]*x[i,j] for j=1:n) <= c[i])
    @constraint(myModel, branch1, x[3,10] <= 0)
    # We have to construct the blocks in a slightly more complicated way since they now are going to contain constraints of differerent types (both <= and >=)
    # We need to make sure that we store "ConstraintRef"s which is the general reference to a constraint. 
    blocks = [ConstraintRef[] for i=1:m]
    for i=1:m
        push!(blocks[i], capacity[i])
    end
    # we add the branching constraint to the subproblem belonging to machine 1
    push!(blocks[3], branch1)

    return myModel, blocks
end

function run()
    # w(i,j) = capacity used when assigning job j to machine i
    # p(i,j) = profit of assigning job j to machine i
    w = [
    8 6 1 7 7 7 5 5 7 7 3 3 8 5 4
    5 5 3 8 5 6 5 9 9 6 1 9 5 6 3
    3 2 4 4 9 1 7 3 3 3 5 3 7 7 1
    6 7 9 9 3 5 2 1 5 5 4 4 6 2 1]
    p = [
    9 9 1 4 6 4 2 3 5 1 9 9 7 6 3
    9 1 8 6 7 8 2 6 6 5 3 4 7 5 3
    4 8 1 8 1 4 1 4 2 6 2 1 2 5 1
    7 9 7 8 3 6 2 3 4 7 9 9 3 9 1]
    c = [22 18 18 19]

    gapModel, blocks = GAPMIPbranch(w, p, c)
    DW_ColGen.DWColGenEasy(gapModel, blocks)
end

run()

# We find the solution is fractional with obj_val = 102.5.

## We can continue branching like this until an optimal solution is found. Right know we have the branching tree shown below.
#                      | UB = 102.5 |
                      /              \
#          x[3,10]=1 /                \ x[3,10]=0 
                    /                  \
#            | UB = 102.25 |      | UB = 102.5 |
# 
#