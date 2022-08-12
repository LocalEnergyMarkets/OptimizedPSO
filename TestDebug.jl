using Pkg
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
include("C:\\Users\\seyednh\\github\\OptimizedPSO\\PSOmodule1.jl")

using BenchmarkTools, DataFrames
const n=Int64(4);
const N=Int64(1*n+1);
using Main.PSO1
##
Ipop=zeros(Float64,n,N);
Vel=zeros(Float64,n,N);
cost=zeros(Float64,N);
xmin=Float64(0);
xmax=Float64(1);
vmax=0.25*(xmax-xmin);

loc_best=0;
Gbest=zeros(Float64,n);
GbestC=100000.0;
Pbest=similar(Ipop);
PbestC=100000*ones(Float64,N);

C1=2
C2=2

i=1
R1=zeros(Float64,n)
R2=zeros(Float64,n)
##
InitPopMatrix!(cost::Vector{Float64},Ipop::Matrix{Float64},Vel::Matrix{Float64},PbestC::Vector{Float64}, Pbest::Matrix{Float64},n::Int64,N::Int64,xmin::Float64,xmax::Float64,vmax::Float64)
GbestC = GlobalBestC!(GbestC,cost)
loc_best = LocBest!(loc_best::Int64, cost::Vector{Float64}, GbestC::Float64)
GlobalBest!(Gbest::Vector{Float64}, loc_best::Int64, Ipop::Matrix{Float64}, n::Int64)
Gbest_init=deepcopy(Gbest)
Ipop_init_D=DataFrame(Ipop, :auto)
Vel_init_D=DataFrame(Vel, :auto)
Cost_init = deepcopy(cost)
Pbest_init=deepcopy(Pbest)
PbestC_init=deepcopy(PbestC)

##

W=rand()*(1-0.4)+0.4

i=1
VelUpdate!(Vel::Matrix{Float64}, n::Int64, i::Int64,W::Float64, C1::Int64, C2::Int64, R1::Vector{Float64}, R2::Vector{Float64}, Pbest::Matrix{Float64}, Gbest::Vector{Float64}, Ipop::Matrix{Float64})
VelMaxCheck!(Vel::Matrix{Float64},vmax::Float64,i::Int64,n::Int64)
VelMinCheck!(Vel::Matrix{Float64},vmax::Float64,i::Int64,n::Int64)
IpopUpdate!(Ipop::Matrix{Float64}, Vel::Matrix{Float64}, i::Int64)
IpopMaxCheck!(Ipop::Matrix{Float64},xmax::Float64,i::Int64, n::Int64)
IpopMinCheck!(Ipop::Matrix{Float64},xmin::Float64,i::Int64, n::Int64)
CostFunc(i::Int64, Ipop::Matrix{Float64}, n::Int64)
PbestUpdate!(PbestC::Vector{Float64}, cost::Vector{Float64}, Pbest::Matrix{Float64}, Ipop::Matrix{Float64}, i::Int64, n::Int64)
GbestUpdate!(GbestC::Float64, Gbest::Vector{Float64}, PbestC::Vector{Float64}, Pbest::Matrix{Float64}, i::Int64, n::Int64)
Gbest_1=deepcopy(Gbest)
Ipop_1_D=DataFrame(Ipop, :auto)
Vel_1_D=DataFrame(Vel, :auto)
Cost_1 = deepcopy(cost)
Pbest_1=deepcopy(Pbest)
PbestC_1=deepcopy(PbestC)




i=2
VelUpdate!(Vel::Matrix{Float64}, n::Int64, i::Int64,W::Float64, C1::Int64, C2::Int64, R1::Vector{Float64}, R2::Vector{Float64}, Pbest::Matrix{Float64}, Gbest::Vector{Float64}, Ipop::Matrix{Float64})
VelMaxCheck!(Vel::Matrix{Float64},vmax::Float64,i::Int64,n::Int64)
VelMinCheck!(Vel::Matrix{Float64},vmax::Float64,i::Int64,n::Int64)
IpopUpdate!(Ipop::Matrix{Float64}, Vel::Matrix{Float64}, i::Int64)
IpopMaxCheck!(Ipop::Matrix{Float64},xmax::Float64,i::Int64, n::Int64)
IpopMinCheck!(Ipop::Matrix{Float64},xmin::Float64,i::Int64, n::Int64)
CostFunc(i::Int64, Ipop::Matrix{Float64}, n::Int64)
PbestUpdate!(PbestC::Vector{Float64}, cost::Vector{Float64}, Pbest::Matrix{Float64}, Ipop::Matrix{Float64}, i::Int64, n::Int64)
GbestUpdate!(GbestC::Float64, Gbest::Vector{Float64}, PbestC::Vector{Float64}, Pbest::Matrix{Float64}, i::Int64, n::Int64)
Gbest_2=deepcopy(Gbest)
Ipop_2_D=DataFrame(Ipop, :auto)
Vel_2_D=DataFrame(Vel, :auto)
Cost_2 = deepcopy(cost)
Pbest_2=deepcopy(Pbest)
PbestC_2=deepcopy(PbestC)



i=3
VelUpdate!(Vel::Matrix{Float64}, n::Int64, i::Int64,W::Float64, C1::Int64, C2::Int64, R1::Vector{Float64}, R2::Vector{Float64}, Pbest::Matrix{Float64}, Gbest::Vector{Float64}, Ipop::Matrix{Float64})
VelMaxCheck!(Vel::Matrix{Float64},vmax::Float64,i::Int64,n::Int64)
VelMinCheck!(Vel::Matrix{Float64},vmax::Float64,i::Int64,n::Int64)
IpopUpdate!(Ipop::Matrix{Float64}, Vel::Matrix{Float64}, i::Int64)
IpopMaxCheck!(Ipop::Matrix{Float64},xmax::Float64,i::Int64, n::Int64)
IpopMinCheck!(Ipop::Matrix{Float64},xmin::Float64,i::Int64, n::Int64)
CostFunc(i::Int64, Ipop::Matrix{Float64}, n::Int64)
PbestUpdate!(PbestC::Vector{Float64}, cost::Vector{Float64}, Pbest::Matrix{Float64}, Ipop::Matrix{Float64}, i::Int64, n::Int64)
GbestUpdate!(GbestC::Float64, Gbest::Vector{Float64}, PbestC::Vector{Float64}, Pbest::Matrix{Float64}, i::Int64, n::Int64)
Gbest_3=deepcopy(Gbest)
Ipop_3_D=DataFrame(Ipop, :auto)
Vel_3_D=DataFrame(Vel, :auto)
Cost_3 = deepcopy(cost)
Pbest_3=deepcopy(Pbest)
PbestC_3=deepcopy(PbestC)




i=4
VelUpdate!(Vel::Matrix{Float64}, n::Int64, i::Int64,W::Float64, C1::Int64, C2::Int64, R1::Vector{Float64}, R2::Vector{Float64}, Pbest::Matrix{Float64}, Gbest::Vector{Float64}, Ipop::Matrix{Float64})
VelMaxCheck!(Vel::Matrix{Float64},vmax::Float64,i::Int64,n::Int64)
VelMinCheck!(Vel::Matrix{Float64},vmax::Float64,i::Int64,n::Int64)
IpopUpdate!(Ipop::Matrix{Float64}, Vel::Matrix{Float64}, i::Int64)
IpopMaxCheck!(Ipop::Matrix{Float64},xmax::Float64,i::Int64, n::Int64)
IpopMinCheck!(Ipop::Matrix{Float64},xmin::Float64,i::Int64, n::Int64)
CostFunc(i::Int64, Ipop::Matrix{Float64}, n::Int64)
PbestUpdate!(PbestC::Vector{Float64}, cost::Vector{Float64}, Pbest::Matrix{Float64}, Ipop::Matrix{Float64}, i::Int64, n::Int64)
GbestUpdate!(GbestC::Float64, Gbest::Vector{Float64}, PbestC::Vector{Float64}, Pbest::Matrix{Float64}, i::Int64, n::Int64)
Gbest_4=deepcopy(Gbest)
Ipop_4_D=DataFrame(Ipop, :auto)
Vel_4_D=DataFrame(Vel, :auto)
Cost_4 = deepcopy(cost)
Pbest_4=deepcopy(Pbest)
PbestC_4=deepcopy(PbestC)


i=5
VelUpdate!(Vel::Matrix{Float64}, n::Int64, i::Int64,W::Float64, C1::Int64, C2::Int64, R1::Vector{Float64}, R2::Vector{Float64}, Pbest::Matrix{Float64}, Gbest::Vector{Float64}, Ipop::Matrix{Float64})
VelMaxCheck!(Vel::Matrix{Float64},vmax::Float64,i::Int64,n::Int64)
VelMinCheck!(Vel::Matrix{Float64},vmax::Float64,i::Int64,n::Int64)
IpopUpdate!(Ipop::Matrix{Float64}, Vel::Matrix{Float64}, i::Int64)
IpopMaxCheck!(Ipop::Matrix{Float64},xmax::Float64,i::Int64, n::Int64)
IpopMinCheck!(Ipop::Matrix{Float64},xmin::Float64,i::Int64, n::Int64)
CostFunc(i::Int64, Ipop::Matrix{Float64}, n::Int64)
PbestUpdate!(PbestC::Vector{Float64}, cost::Vector{Float64}, Pbest::Matrix{Float64}, Ipop::Matrix{Float64}, i::Int64, n::Int64)
GbestUpdate!(GbestC::Float64, Gbest::Vector{Float64}, PbestC::Vector{Float64}, Pbest::Matrix{Float64}, i::Int64, n::Int64)
Gbest_5=deepcopy(Gbest)
Ipop_5_D=DataFrame(Ipop, :auto)
Vel_5_D=DataFrame(Vel, :auto)
Cost_5 = deepcopy(cost)
Pbest_5=deepcopy(Pbest)
PbestC_5=deepcopy(PbestC)
