using Pkg
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
include("C:\\Users\\seyednh\\github\\OptimizedPSO\\PSOmodule.jl")

using BenchmarkTools
const n=Int64(1000);
const N=Int64(10*n);
using Main.PSO
##  #TODO: Think about static arrays
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

W=rand()*(1-0.4)+0.4;
C1=2
C2=2

i=1
R1=zeros(Float64,n)
R2=zeros(Float64,n)
## Functions
InitPopMatrix!(Ipop::Matrix{Float64},Vel::Matrix{Float64},n::Int64,N::Int64,xmin::Float64,xmax::Float64,vmax::Float64)
CostAll!(cost::Vector{Float64}, Ipop::Matrix{Float64},N::Int64,n::Int64)
VelAll!(N::Int64,n::Int64, W::Float64, Vel::Matrix{Float64}, C1::Int64, C2::Int64, R1::Vector{Float64}, R2::Vector{Float64}, Pbest::Matrix{Float64}, Gbest::Vector{Float64}, Ipop::Matrix{Float64}, vmax::Float64)
PopAll!(N::Int64,n::Int64, Vel::Matrix{Float64}, Ipop::Matrix{Float64}, xmax::Float64, xmin::Float64)
## PSO algorithm
InitPopMatrix!(Ipop::Matrix{Float64},Vel::Matrix{Float64},n::Int64,N::Int64,xmin::Float64,xmax::Float64,vmax::Float64)
