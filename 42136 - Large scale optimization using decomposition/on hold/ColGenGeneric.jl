include("JumpModelToMatrix.jl")
include("DecompMatrix.jl")

module DW_ColGen

export GAPMIP, DWColGenEasy

using JuMP, GLPK
using LinearAlgebra, SparseArrays
using ..DecompMatrix, ..JumpModelToMatrix

EPSVAL = 0.00001

function setupSub(sub::JuMP.Model, Asub, bsub, sense, subvars, varLB, varUB, vecIsInt)
    (mAsub, n) = size(Asub)

    @variable(sub, xVars[1:n])
    # define upper and lower bounds as well as integrality constraints
    for j=1:n
         set_lower_bound(xVars[j], varLB[subvars[j]])
         set_upper_bound(xVars[j], varUB[subvars[j]])
         if vecIsInt[subvars[j]]
             set_integer(xVars[j])
         end
     end

    # objective is not here. We define once dual variables become known
    @objective(sub, Max, 0 )

    for i=1:mAsub
        if sense[i]==LEQ
            @constraint(sub, dot(Asub[i,:],xVars) <= bsub[i] )
        elseif sense[i]==GEQ
            @constraint(sub, dot(Asub[i,:],xVars) >= bsub[i] )
        else
            @constraint(sub, dot(Asub[i,:],xVars) == bsub[i] )
        end
    end
    return xVars
end

function setupMaster(master::JuMP.Model, A0, b0, nSub, sense)
    (mA0, nA0) = size(A0[1])
    K = 1

    # In this case we do not use a starting set of extreme points.
    # we just use a dummy starting column

    @variable(master, lambda[1:K] >= 0 )
    # Remember to consider if we need to maximize or minimize
    @objective(master, Max, sum(- 1000 * lambda[j] for j=1:K) )
    JuMPConstraintRef = JuMP.ConstraintRef
    consref = Vector{JuMPConstraintRef}(undef,mA0)
    for i=1:mA0
        if sense[i]==LEQ
            consref[i] = @constraint(master, b0[i]*lambda[1] <= b0[i] )
        elseif sense[i]==GEQ
            consref[i] = @constraint(master, b0[i]*lambda[1] >= b0[i] )
        else
            consref[i] = @constraint(master, b0[i]*lambda[1] == b0[i] )
        end
    end
    @constraint(master, convexityCons[k=1:nSub], lambda[1]  == 1)

    return consref, convexityCons, lambda
end

# x is the new extreme point we wish to add to the master problem
function addColumnToMaster(master::JuMP.Model, cVec, A0, x, consRef, convexityCons, subNumber)
    (mA0, n) = size(A0[subNumber])
    A0x = A0[subNumber]*x

    oldvars = JuMP.all_variables(master)
    new_var = @variable(master, base_name="lambda_$(length(oldvars))", lower_bound=0)
    JuMP.set_objective_coefficient(master, new_var, dot(cVec[subNumber],x))

    for i=1:mA0
        # only insert non-zero elements (this saves memory and may make the master problem easier to solve)
        if A0x[i] != 0
            set_normalized_coefficient(consRef[i], new_var, A0x[i])
        end
    end
    # add variable to convexity constraint.
    set_normalized_coefficient(convexityCons[subNumber], new_var, 1)

    return new_var
end

function solveSub(sub, myPi, myKappa, cVec, A0, xVars,subNumber)
    c_k = cVec[subNumber]
    A0Sub = A0[subNumber]
    xVarsSub = xVars[subNumber]

    # set objective.
    # NOTE: Remember to consider if maximization or minimization is needed
    @objective(sub[subNumber], Max, dot(c_k,xVarsSub) - dot(myPi*A0Sub,xVarsSub) - myKappa[subNumber] )

    optimize!(sub[subNumber])
    println(termination_status(sub[subNumber]))
    if termination_status(sub[subNumber]) != MOI.OPTIMAL
        throw("Error: Non-optimal sub-problem status")
    end

    return JuMP.objective_value(sub[subNumber]), value.(xVarsSub)
end

function DWColGen(A0,ASub,b0,bSub, senseA0, senseSub, pPerSub,subvars, mip::MIP)
    varNames = mip.varNames
    varLB = mip.varLB
    varUB = mip.varUB
    vecIsInt = mip.vecIsInt

    nSub = length(bSub)
    sub = Vector{JuMP.Model}(undef, nSub)

    for k=1:nSub
        sub[k] = Model(GLPK.Optimizer)
    end
    master = Model(GLPK.Optimizer)
    xVars = []
    for k=1:nSub
        push!(xVars, setupSub(sub[k], ASub[k], bSub[k], senseSub[k], subvars[k], varLB, varUB, vecIsInt))
    end
    (consRef, convexityCons, lambdas) = setupMaster(master, A0, b0, nSub, senseA0)
    # extremePoints records the extreme points. extremePoints[p] corresponds to variable lambda[p]
    # extremePointForSub[p] records which sub-problem the p'th extreme point "belongs to"
    # first extreme point is a dummy one.
    extremePoints = [[]]
    extremePointForSub = [-1]
    done = false
    iter = 1
    while !done
        optimize!(master)
        if termination_status(master) != MOI.OPTIMAL
            throw("Error: Non-optimal master-problem status")
        end
        myPi = -dual.(consRef)
        # ensure that myPi and myKappa are  row vectors
        myPi = reshape(myPi, 1, length(myPi))
        myKappa = -dual.(convexityCons)
        myKappa = reshape(myKappa, 1, length(myKappa))
        println("myPi = $myPi")
        println("myKappa = $myKappa")
        done = true
        println("iteration: $iter, objective value = $(JuMP.objective_value(master))")
        bestRedCost = -1
        for k=1:nSub
            redCost, xVal = solveSub(sub, myPi, myKappa, pPerSub, A0, xVars, k)
            print("sub $k red cost = $redCost, ")
            if redCost > bestRedCost
                bestRedCost = redCost
            end
            if redCost > EPSVAL
                newVar = addColumnToMaster(master, pPerSub, A0, xVal, consRef, convexityCons,k)
                push!(lambdas,newVar)
                push!(extremePoints, xVal)
                push!(extremePointForSub, k)
                done = false
            end
        end
        iter += 1
        println("best reduced cost = $bestRedCost")
    end
    println("Done after $(iter-1) iterations. Objective value = $(JuMP.objective_value(master))")
    # compute values of original variables
    origVarValSub = []
    for s = 1:length(subvars)
        push!(origVarValSub,zeros(length(subvars[s])))
    end
    lambdaVal = value.(lambdas)
    for p=1:length(lambdaVal)
        if lambdaVal[p] > 0.0001
            println("lambda_$p=", lambdaVal[p], ", sub=$(extremePointForSub[p]), extr.point=$(extremePoints[p])")
            origVarValSub[extremePointForSub[p]] += lambdaVal[p]*extremePoints[p]
        end
    end
    for s = 1:length(subvars)
        #println("var val for sub problem $s: $(origVarValSub[s])")
        for t=1:length(origVarValSub[s])
            if abs(origVarValSub[s][t]) > 0.0001
                println("$(varNames[subvars[s][t]])=$(origVarValSub[s][t])")
            end
        end
    end
    
end

function convertBlocks(blocks, constraintRefToRowIdDict)
    blockInt = Vector{Vector{Int64}}()
    for block in blocks
        push!(blockInt, Vector{Int64}())
        for constraintRef in block
            push!(blockInt[end], constraintRefToRowIdDict[constraintRef])
        end
        println("block $(length(blockInt)): ", blockInt[end])
    end
    return blockInt
end

function DWColGenEasy(jumpModel, blocks::Vector{Vector{ConstraintRef}})
    mip, constraintRefToRowIdDict = getConstraintMatrix(jumpModel)
    blocksAsRowIdx = convertBlocks(blocks, constraintRefToRowIdDict)
    A0, b0, A0Sense, Asub, bSub, subSense, subVars, cPerSub = constructSubMatrices(mip, blocksAsRowIdx)

    timeStart = time()
    DWColGen(A0,Asub,b0,bSub, A0Sense, subSense, cPerSub,subVars,mip)
    println("Elapsed time: $(time()-timeStart) seconds")
end

function GAPMIP(w, p, c)
    (m,n) = size(w)
    myModel = Model(GLPK.Optimizer)
    @variable(myModel, x[1:m,1:n], Bin)
    @objective(myModel, Max, sum(p[i,j]*x[i,j] for i=1:m for j=1:n))
    # each job must be served
    @constraint(myModel, [j=1:n],sum(x[i,j] for i=1:m) == 1)
    @constraint(myModel, capacity[i=1:m],sum(w[i,j]*x[i,j] for j=1:n) <= c[i])

    # define blocks (Each block becomes a sub-problem)
    blocks = [[capacity[i]] for i=1:m]

    return myModel, blocks
end

function test2()
    # w(i,j) = capacity used when assigning job j to machine i
    # p(i,j) = profit of assigning job j to machine i
    w = [
    8 6 1 7 7 7 5 5 7 7 3 3 8 5 4 ;
    5 5 3 8 5 6 5 9 9 6 1 9 5 6 3 ;
    3 2 4 4 9 1 7 3 3 3 5 3 7 7 1 ;
    6 7 9 9 3 5 2 1 5 5 4 4 6 2 1
    ]
    p =[
    9 9 1 4 6 4 2 3 5 1 9 9 7 6 3 ;
    9 1 8 6 7 8 2 6 6 5 3 4 7 5 3 ;
    4 8 1 8 1 4 1 4 2 6 2 1 2 5 1 ;
    7 9 7 8 3 6 2 3 4 7 9 9 3 9 1]
    cap = [22 18 18 19]

    gapModel, blocks = GAPMIP(w, p, cap)
    DWColGenEasy(gapModel, blocks)
end

function exampleSlides()
    # w(i,j) = capacity used when assigning job j to machine i
    # p(i,j) = profit of assigning job j to machine i
    w = [
    7 2 8;
    8 7 6;
    9 1 9;
    ]
    p =[
    6 4 6
    1 3 4
    1 2 8
    ]
    cap = [9 7 10]

    gapModel, blocks = GAPMIP(w, p, cap)
    DWColGenEasy(gapModel, blocks)
end


#test2()

end

# run with
#DW_ColGen.test2()
#DW_ColGen.exampleSlides()
