module JumpModelToMatrix

using JuMP, SparseArrays

export getConstraintMatrix, Sense, GenLinConstr, LEQ, EQ, GEQ, MIP

# constraint types:
@enum Sense LEQ=1 GEQ=2 EQ=3

struct GenLinConstr
    vars::Vector{Int64}
    coeffs::Vector{Float64}
    sense::Sense
    rhs::Float64
    #index::Int64
	ref::ConstraintRef
end

mutable struct MIP
    A::SparseMatrixCSC{Float64,Int64}
    b::Vector{Float64}
    c::Vector{Float64}
    varLB::Vector{Float64}
    varUB::Vector{Float64}
    vecIsInt::Vector{Bool}
    consType::Vector{Sense}
    varNames::Vector{String}
    minimization::Bool
    fileName::String
end


function parseAffExpr(affExpr::GenericAffExpr{Float64,VariableRef})
    vecVars = Vector{Int64}()
    vecCoeffs = Vector{Float64}()
    for varMultPair in affExpr.terms
        push!(vecVars, varMultPair[1].index.value)
        push!(vecCoeffs, varMultPair[2])
    end
    return vecVars, vecCoeffs, affExpr.constant
end

function getLinearConstraints(constraintRefs)
    vecConstraints = Vector{GenLinConstr}()
    for constref in constraintRefs
        constr = constraint_object(constref)
        vecVars, vecCoeffs, constantTerm = parseAffExpr(constr.func)
        #println("Parsing $constref")
        #println("Got: vecVars=$vecVars, vecCoeffs=$vecCoeffs, constantTerm=$constantTerm")
        if constantTerm != 0.0
            throw("Expeced a constraint function without a constant term. Stefan needs to fix this :-)")
        end
        if typeof(constr.set) == MOI.LessThan{Float64}
            sense = LEQ
            RHS = constr.set.upper
        end
        if typeof(constr.set) == MOI.EqualTo{Float64}
            sense = EQ
            RHS = constr.set.value
        end
        if typeof(constr.set) == MOI.GreaterThan{Float64}
            sense = GEQ
            RHS = constr.set.lower
        end

        linConstr = GenLinConstr(vecVars, vecCoeffs,sense, RHS, constref)
        #println(linConstr)
        push!(vecConstraints, linConstr)
    end
    return vecConstraints
end


function getMatrixForm(vecLinConstraints::Vector{GenLinConstr})
    b = Vector{Float64}(undef, length(vecLinConstraints))
    sense = Vector{Sense}(undef, length(vecLinConstraints))
    constraintRefToRowIdDict = Dict{ConstraintRef, Int64}()
    AColIdx = Vector{Int64}()
    ARowIdx = Vector{Int64}()
    ANZval = Vector{Float64}()

    for i=1:length(vecLinConstraints)
        #index = vecLinConstraints[i].index
        #rowNo = index-minIndex+1
        b[i] = vecLinConstraints[i].rhs
        sense[i] = vecLinConstraints[i].sense
		constraintRefToRowIdDict[vecLinConstraints[i].ref] = i
        vecVars = vecLinConstraints[i].vars
        vecCoeffs = vecLinConstraints[i].coeffs
        for j=1:length(vecVars)
            push!(ARowIdx, i)
            push!(AColIdx, vecVars[j])
            push!(ANZval, vecCoeffs[j])
        end
    end
    #println("ARowIdx: ", ARowIdx)
    #println("AColIdx: ", AColIdx)
    #println("ANZval: ", ANZval)
    A = sparse(ARowIdx, AColIdx, ANZval)
    #println("A:")
    #I,J,V = findnz(A)
    #println("I: ", I)
    #println("J: ", J)
    #println("V: ", V)
    return A, b, sense, constraintRefToRowIdDict
end

function getConstraintMatrix(model)
    #println("getConstraintMatrix input model")
    #println(model)
    #println("Constraint types: ", list_of_constraint_types(model))
    # Get all variables used in the model
    vars = JuMP.all_variables(model)
    varIndices = Vector{Int64}()
    varNames = Vector{String}()
    for varRef in vars
        push!(varIndices, varRef.index.value)
        push!(varNames, name(varRef))
    end
    #println("varIndices: $varIndices")
    #println("varNames: $varNames")

    minIndex = minimum(varIndices)
    maxIndex = maximum(varIndices)
    for i=1:maxIndex
        if varIndices[i] != i
            throw("getConstraintMatrix(...): varIndices array is out of order or is missing entries. This was not expected  :-(")
        end
    end
    nInt = num_constraints(model, VariableRef, MOI.Integer)
    println("number of integer variables: $nInt")
    nBin = num_constraints(model, VariableRef, MOI.ZeroOne)
    println("number of binary variables: $nBin")

    vecVarLB = -Inf*ones(Float64, maxIndex)
    vecVarUB = Inf*ones(Float64, maxIndex)
    println("Number of variable lower bounds: $(num_constraints(model, VariableRef, MOI.GreaterThan{Float64}))")
    println("Number of variable upper bounds: $(num_constraints(model, VariableRef, MOI.LessThan{Float64}))")
    println("Number of variable interval constraints: $(num_constraints(model, VariableRef, MOI.Interval{Float64}))")

    for varLBConsRef in all_constraints(model, VariableRef, MOI.GreaterThan{Float64})
         constr = constraint_object(varLBConsRef)
         varRef = constr.func
         LB = constr.set.lower
         vecVarLB[varRef.index.value] = LB
     end
     for varUBConsRef in all_constraints(model, VariableRef, MOI.LessThan{Float64})
          constr = constraint_object(varUBConsRef)
          varRef = constr.func
          UB = constr.set.upper
          vecVarUB[varRef.index.value] = UB
      end
      for varIntervalConsRef in all_constraints(model, VariableRef, MOI.Interval{Float64})
           constr = constraint_object(varIntervalConsRef)
           varRef = constr.func
           LB = constr.set.lower
           UB = constr.set.upper
           vecVarLB[varRef.index.value] = LB
           vecVarUB[varRef.index.value] = UB
       end
      vecIsInt = [false for i=1:maxIndex]
      for intConstraintRef in all_constraints(model, VariableRef, MOI.Integer)
          constr = constraint_object(intConstraintRef)
          vecIsInt[constr.func.index.value] = true
      end
      for binConstraintRef in all_constraints(model, VariableRef, MOI.ZeroOne)
          constr = constraint_object(binConstraintRef)
          vecIsInt[constr.func.index.value] = true
          # binary variables does not occur in the (VariableRef, MOI.GreaterThan{Float64}) and (VariableRef, MOI.LessThan{Float64}) constraints
          # so we have to set their bounds explicitely
          vecVarLB[constr.func.index.value] = 0
          vecVarUB[constr.func.index.value] = 1
      end

      # get all the <= constraints
      vec1 = getLinearConstraints(all_constraints(model, GenericAffExpr{Float64,VariableRef}, MOI.LessThan{Float64}))
      vec2 = getLinearConstraints(all_constraints(model, GenericAffExpr{Float64,VariableRef}, MOI.EqualTo{Float64}))
      vec3 = getLinearConstraints(all_constraints(model, GenericAffExpr{Float64,VariableRef}, MOI.GreaterThan{Float64}))
      vecLinConstraints = vcat(vec1,vec2,vec3)
      objVars, objCoeffs, objConst =  parseAffExpr(objective_function(model))
      if objConst != 0.0
          throw("getConstraintMatrix(...): Didn't expect objective function to have a constant term!")
      end
      vecObj = zeros(Float64, maxIndex)
      for i=1:length(objVars)
          vecObj[objVars[i]] = objCoeffs[i]
      end
      A, b, consType, constraintRefToRowIdDict = getMatrixForm(vecLinConstraints)
      #println("b=$b")
      #println("consType=$consType")
      mip = MIP(A,b,vecObj, vecVarLB, vecVarUB, vecIsInt, consType, varNames, objective_sense(model)==MOI.MIN_SENSE, "")
      return mip, constraintRefToRowIdDict
end



end
