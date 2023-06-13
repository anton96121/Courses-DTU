# Benders template (no rays)

#-------------------------------------------------------------------
# Intro definitions
using JuMP
using CPLEX
using LinearAlgebra
using NamedArrays
using LatexPrint
cd("C:\\Users\\anton\\OneDrive - Københavns Universitet\\Uni\\Uni\\10. semester\\decomposition\\project 2")
include("JumpModelToMatrix.jl")
include("DecompMatrix.jl")
include("ColGenGeneric.jl")
#-------------------------------------------------------------------

P = [  0 0 0 0 0 95 0 0 0 95 ]'
D = [  5 35 25 0 50 0 20 30 25 0 ]'


 #-------------------------------------------------------------------
# PARAMETERS
Nodes     = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "j"]
N=length(Nodes)
   node_x = [ 1    3    8    2    4    5    6    7    9    9 ]
   node_y = [ 1    2    1    7    9    5    3    7    9    4 ]
gas_arc = [
              0    1    0    0    0    0    0    0    0    0;
              0    0    0    0    0    1    0    0    0    0;
              0    0    0    0    0    0    0    0    0    1;
              1    0    0    0    1    0    0    0    0    0;
              0    0    0    0    0    1    0    0    0    0;
              1    0    0    1    0    0    0    1    0    0;
              0    1    1    0    0    0    0    0    0    0; 
              0    0    0    0    0    0    0    0    1    0;
              0    0    0    0    1    0    0    0    0    0;
              0    0    0    0    0    1    1    0    1    0
]

c=zeros(Float64,N,N)
for n1=1:N
    for n2=1:N
        c[n1,n2]=floor(sqrt( (node_x[n1]-node_x[n2])*(node_x[n1]-node_x[n2]) +
                          (node_y[n1]-node_y[n2])*(node_y[n1]-node_y[n2]) ) )
    end
end
f=zeros(Float64,N,N)
for n1=1:N
    for n2=1:N
        f[n1,n2]=10*c[n1,n2]*(1-gas_arc[n1,n2])
    end
end

function Abeta(N)
    A = zeros(N,N*N);

    for i=1:N
        for j=1:N
            A[i,(j-1)*N+i] += 1
            A[i,(i-1)*N+j] -= 1
        end
    end

    return A
end

function GAPMIPbranch(node_x, node_y, gas_arc,D,P,N)

    myModel = Model(CPLEX.Optimizer)
    @variable(myModel, 0 <= x[1:N^2]<=100000)
    #@variable(myModel, y[1:N^2], Bin)
    @variable(myModel, 0<=y[1:N^2]<=1, Int)

    # Generate coefficents
    c=zeros(Float64,N,N)
    for n1=1:N
        for n2=1:N
            c[n1,n2]=floor(sqrt( (node_x[n1]-node_x[n2])*(node_x[n1]-node_x[n2]) +
                              (node_y[n1]-node_y[n2])*(node_y[n1]-node_y[n2]) ) )
        end
    end
    f=zeros(Float64,N,N)
    for n1=1:N
        for n2=1:N
            f[n1,n2]=10*c[n1,n2]*(1-gas_arc[n1,n2])
        end
    end
    A = zeros(N,N*N);

    for i=1:N
        for j=1:N
            A[i,(j-1)*N+i] += 1
            A[i,(i-1)*N+j] -= 1
        end
    end

    b0 = zeros(Float64,N*N,1)
    b1 = D-P

    f_vec = vec(f')
    c_vec = vec(c')

    M = 1000
    B0 = M*I(N*N)
    B1 = zeros(Float64,N,N*N)
    A0 = -I(N*N)
    A1 = A

    @objective(myModel, Min, -[c_vec;f_vec]'*[x;y])

    @constraint(myModel,  y[52]>= 1)
    @constraint(myModel,  y[55]>= 1)
    @constraint(myModel,  y[93]>= 1)
    @constraint(myModel,  y[98]>= 1)


    # each job must be served
    @constraint(myModel, con0[i=1:N^2], ([A0[i,:]; B0[i,:]]'*[x;y]) >= b0[i])
    @constraint(myModel, con12[i=1:N^2], y[i] .>= vec(gas_arc')[i])
    @constraint(myModel, con1[i=1:N], ([A1[i,:]; B1[i,:]]'*[x;y]) >= b1[i])
    # We have to construct the blocks in a slightly more complicated way since they now are
    # going to contain constraints of differerent types (both <= and >= once we start branching)
    # We need to make sure that we store "ConstraintRef"s which is the general reference to a constraint.
    blocks = [ConstraintRef[] for i=1:1]
    for i=1:N^2
        push!(blocks[1], con0[i])
        push!(blocks[1], con12[i])
    end

    return myModel, blocks
end

function run()
    P = [  0 0 0 0 0 95 0 0 0 95 ]'
    D = [  5 35 25 0 50 0 20 30 25 0 ]'
    Nodes     = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "j"]
    N=length(Nodes)
    node_x =  [1    3    8    2    4    5    6    7    9    9]
    node_y =  [1    2    1    7    9    5    3    7    9    4]
    gas_arc = [0    1    0    0    0    0    0    0    0    0;
               0    0    0    0    0    1    0    0    0    0;
               0    0    0    0    0    0    0    0    0    1;
               1    0    0    0    1    0    0    0    0    0;
               0    0    0    0    0    1    0    0    0    0;
               1    0    0    1    0    0    0    1    0    0;
               0    1    1    0    0    0    0    0    0    0; 
               0    0    0    0    0    0    0    0    1    0;
               0    0    0    0    1    0    0    0    0    0;
               0    0    0    0    0    1    1    0    1    0]

    gapModel, blocks = GAPMIPbranch(node_x, node_y, gas_arc,D,P,N)
    DW_ColGen.DWColGenEasy(gapModel, blocks)
end

run()

