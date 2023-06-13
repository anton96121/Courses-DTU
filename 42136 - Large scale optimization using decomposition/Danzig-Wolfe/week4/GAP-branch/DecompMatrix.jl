module DecompMatrix

using SparseArrays
using ..JumpModelToMatrix

export constructSubMatrices

function addVarsFromRow(Arow, varBitSet)
    #println("typeof(Arow): ", typeof(Arow))
    #println("Arow: ", Arow, ", Arow.n: ",Arow.n, "Arow.nzval: ",Arow.nzval, "Arow.nzind: ",Arow.nzind)
    for idx=1:length(Arow.nzval)
        if Arow.nzval[idx] != 0
            push!(varBitSet, Arow.nzind[idx])
        end
    end
end

function constructSubMatrices(mip::MIP,blocks)
    A = mip.A
    b = mip.b
    c = mip.c
    consType = mip.consType
    #println("blocks = $blocks")
    # find out which variables that belong to which block:
    #println("A:",A)
    #=for i=1:length(mip.varNames)
        println("var $i: ", mip.varNames[i])
    end=#
    nBlocks = length(blocks)
    if nBlocks < 1
        throw("error. We expect at least one block")
    end
    (m,n) = size(A)
    varsInSub=Vector{BitSet}(undef,nBlocks)
    for k=1:nBlocks
        #println("Processing block $k: $(blocks[k])")
        varsInSub[k] = BitSet()
        for row in blocks[k]
            #println("row $row:", A[row,:])
            addVarsFromRow(A[row,:], varsInSub[k])
        end
    end
    # check that the blocks are not sharing any variables
    for k1=1:nBlocks
        for k2=k1+1:nBlocks
            if !isempty(intersect(varsInSub[k1],varsInSub[k2]))
                println("Vars in sub $k1: ", varsInSub[k1])
                println("Vars in sub $k2: ", varsInSub[k2])
                throw("Blocks $k1 and $k2 share variables. This is not supported.")
            end
        end
    end
    #println("varsInSub: $varsInSub")
    # Check that all variables belong to a block
    setUnion = deepcopy(varsInSub[1])
    for k=2:nBlocks
         union!(setUnion, varsInSub[k])
     end
     for j=1:n
         if !(j in setUnion)
             throw("variable $j does not belong to any block. This is not supported")
         end
    end
    # Find the rows that are left in master problem:
    masterRows = BitSet(1:m)
    for block in blocks
        setdiff!(masterRows, block)
    end
    subVars = Vector{Vector{Int64}}(undef,nBlocks)
    A0 = Vector{SparseMatrixCSC{Float64,Int64}}(undef,nBlocks)
    for k=1:nBlocks
        subVars[k] = collect(varsInSub[k])
        A0[k] = deepcopy(A[collect(masterRows),subVars[k]])
    end
    b0 = deepcopy(b[collect(masterRows)])
    A0Sense = deepcopy(consType[collect(masterRows)])
    # find matrices for the sub-problems
    ASub = Vector{SparseMatrixCSC{Float64,Int64}}(undef,nBlocks)
    bSub = Vector{Vector{Float64}}(undef,nBlocks)
    cPerSub = Vector{Vector{Float64}}(undef,nBlocks)
    senseSub = Vector{Vector{Sense}}(undef,nBlocks)
    for k=1:nBlocks
        #println("k=$k, blocks[k]=$(blocks[k])")
        ASub[k] = deepcopy(A[blocks[k],subVars[k]])
        bSub[k] = deepcopy(b[blocks[k]])
        cPerSub[k] = deepcopy(c[subVars[k]])
        senseSub[k] = deepcopy(consType[blocks[k]])
    end
    return A0, b0, A0Sense, ASub,bSub, senseSub, subVars, cPerSub
end


end
