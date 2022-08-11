using Pkg
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
include("C:\\Users\\seyednh\\github\\OptimizedPSO\\PSOmodule1.jl")

using BenchmarkTools
const n=Int64(4);
const N=Int64(10*n);
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

for itr in 1:100
    W=rand()*(1-0.4)+0.4;
    for i in 1:N
        VelUpdate!(Vel::Matrix{Float64}, n::Int64, i::Int64,W::Float64, C1::Int64, C2::Int64, R1::Vector{Float64}, R2::Vector{Float64}, Pbest::Matrix{Float64}, Gbest::Vector{Float64}, Ipop::Matrix{Float64})
        VelMaxCheck!(Vel::Matrix{Float64},vmax::Float64,i::Int64,n::Int64)
        VelMinCheck!(Vel::Matrix{Float64},vmax::Float64,i::Int64,n::Int64)

        IpopUpdate!(Ipop::Matrix{Float64}, Vel::Matrix{Float64}, i::Int64)
        IpopMaxCheck!(Ipop::Matrix{Float64},xmax::Float64,i::Int64, n::Int64)
        IpopMinCheck!(Ipop::Matrix{Float64},xmin::Float64,i::Int64, n::Int64)

        CostFunc(i::Int64, Ipop::Matrix{Float64}, n::Int64)

        PbestUpdate!(PbestC::Vector{Float64}, cost::Vector{Float64}, Pbest::Matrix{Float64}, Ipop::Matrix{Float64}, i::Int64, n::Int64)

        GbestUpdate!(PbestC::Vector{Float64}, Pbest::Matrix{Float64}, GbestC::Float64, Gbest::Vector{Float64}, i::Int64, n::Int64)
    end
end
