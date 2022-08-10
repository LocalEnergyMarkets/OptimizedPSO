using Pkg
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
using BenchmarkTools
const n=Int64(1000);
const N=Int64(10*n);
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

#=function FullFunc(Ipop,Vel,cost,xmin,xmax,vmax)

  InitPopMatrix!(Ipop::Matrix{Float64},Vel::Matrix{Float64},n::Int64,N::Int64,xmin::Float64,xmax::Float64,vmax::Float64)
  test1!(cost::Vector{Float64}, Ipop::Matrix{Float64})

  return Ipop, Vel, cost
end

Ipop, Vel, cost=@btime FullFunc(Ipop,Vel,cost,xmin,xmax,vmax)
=#

function InitPopMatrix!(Ipop::Matrix{Float64},Vel::Matrix{Float64},n::Int64,N::Int64,xmin::Float64,xmax::Float64,vmax::Float64)
  for i=1:N
    for j=1:n
      Ipop[j,i]=rand()*(xmax-xmin)+xmin;
      Vel[j,i]=rand()*(2*vmax)-vmax;
    end
  end
  #return Ipop, Vel
end
InitPopMatrix!(Ipop::Matrix{Float64},Vel::Matrix{Float64},n::Int64,N::Int64,xmin::Float64,xmax::Float64,vmax::Float64)
##
function CostFunc(i::Int64, Ipop::Matrix{Float64})
  c = 0
  Ip_i=@view Ipop[:,i]
  for k in 1:n
    c+=exp(-Ip_i[k])*sin(2*pi*Ip_i[k])
  end      # Ipop_i=Ipop[:,i]
  c
end

##
#@btime for i=1:N
#  cost[i] = CostFunc(i::Int64,Ipop::Matrix{Float64})
#end

function test1!(cost::Vector{Float64}, Ipop::Matrix{Float64})
  for i=1:N
    cost[i] = CostFunc(i::Int64,Ipop::Matrix{Float64})
  end
end
test1!(cost::Vector{Float64}, Ipop::Matrix{Float64})
##
W=rand()*(1-0.4)+0.4;
C1=2
C2=2

i=1
R1=rand(n)
R2=rand(n)

function VelUpdate(i::Int64,W::Float64, Vel::Matrix{Float64}, C1::Int64, C2::Int64, R1::Vector{Float64}, R2::Vector{Float64}, Pbest::Matrix{Float64}, Gbest::Vector{Float64}, Ipop::Matrix{Float64})
#  Vel[:,i] = W * Vel[:,i] + C1 * R1.*(Pbest[:,i] - Ipop[:,i]) + C2 * R2 .* (Gbest - Ipop[:,i])
  @. Vel[:,i] = W* (@view Vel[:,i]) + C1 * R1*((@view Pbest[:,i]) - (@view Ipop[:,i])) + C2 * R2 * (Gbest - (@view Ipop[:,i]))
end
VelUpdate(i::Int64,W::Float64, Vel::Matrix{Float64}, C1::Int64, C2::Int64, R1::Vector{Float64}, R2::Vector{Float64}, Pbest::Matrix{Float64}, Gbest::Vector{Float64}, Ipop::Matrix{Float64})
##

function VelMaxCheck!(Vel::Matrix{Float64},vmax::Float64,i::Int64)
  V = @view Vel[:,i]
  for j=1:n
    if V[j]>vmax
      V[j]=vmax
    end
    #V[j] = minimum([V[j],vmax])
  end
end
VelMaxCheck!(Vel::Matrix{Float64},vmax::Float64,i::Int64)

function VelMinCheck!(Vel::Matrix{Float64},vmax::Float64,i::Int64)
  V = @view Vel[:,i]
  for j=1:n
    if V[j]<-vmax
      V[j]=-vmax
    end
    #V[j] = minimum([V[j],vmax])
  end
end
VelMinCheck!(Vel::Matrix{Float64},vmax::Float64,i::Int64)
##
function IpopUpdate!(Ipop::Matrix{Float64}, Vel::Matrix{Float64}, i::Int64)
  @. Ipop[:,i] = (@view Ipop[:,i]) + (@view Vel[:,i])
end
IpopUpdate!(Ipop::Matrix{Float64}, Vel::Matrix{Float64}, i::Int64)

function IpopMaxCheck!(Ipop::Matrix{Float64},xmax::Float64,i::Int64)
  I = @view Ipop[:,i]
  for j=1:n
    if I[j]>xmax
      I[j]=xmax
    end
    #V[j] = minimum([V[j],vmax])
  end
end
IpopMaxCheck!(Ipop::Matrix{Float64},xmax::Float64,i::Int64)

function IpopMinCheck!(Ipop::Matrix{Float64},xmin::Float64,i::Int64)
  I = @view Ipop[:,i]
  for j=1:n
    if I[j]<xmin
      I[j]=xmin
    end
    #V[j] = minimum([V[j],vmax])
  end
end
IpopMinCheck!(Ipop::Matrix{Float64},xmin::Float64,i::Int64)
##
Ip=Ipop[:,i]
function PbestUpdate!(PbestC::Vector{Float64}, cost::Vector{Float64}, Pbest::Matrix{Float64}, Ip::Vector{Float64},i::Int64)
  if PbestC[i]>=cost[i]
    PbestC[i] = cost[i]
    Pbest[:,i] = Ip
  end
  return Pbest
end
PbestUpdate!(PbestC::Vector{Float64}, cost::Vector{Float64}, Pbest::Matrix{Float64}, Ip::Vector{Float64},i::Int64)
##
GbestC = @btime minimum(cost)
loc_best=0
function FindLocBest!(loc_best::Int64, cost::Vector{Float64}, GbestC::Float64)
  loc_best = findfirst(x -> x==GbestC,cost)[1];
  return loc_best
end
loc_best = @btime FindLocBest!(loc_best::Int64, cost::Vector{Float64}, GbestC::Float64)

Gbest=Ipop[:,loc_best]

#=

i=1
PbestC[i]>=cost[i]

CN = 1==1
function PbestUpdate(i::Int64,PbestC::Vector{Float64}, cost::Vector{Float64}, Pbest::Matrix{Float64}, Ipop::Matrix{Float64})
  if PbestC[i]>=cost[i]
    PbestC[i] = cost[i]
    Pbest[:,i] = Ipop[:,i]
  end
end
@btime for i in 1:N
  PbestUpdate(i::Int64,PbestC::Vector{Float64}, cost::Vector{Float64}, Pbest::Matrix{Float64}, Ipop::Matrix{Float64})
end


function evaluator!(loc_best::Int64,Gbest::Vector{Float64},GbestC::Float64,Pbest::Vector{Float64},PbestC::Float64,cost::Vector{Float64})

  GbestC = minimum(cost)
  loc_best=findfirst(x -> x==minimum(cost),cost)[1];
  Gbest=Ipop[:,loc_best];

  PbestC
  Pbest=deepcopy(Ipop);
end
@btime evaluator!(loc_best::Int64,Gbest::Vector{Float64},GbestC::Float64,Pbest::Vector{Float64},PbestC::Float64,cost::Vector{Float64})





function evaluator2!(loc_best::Int64,Gbest::Vector{Float64},GbestC::Float64,Pbest::Vector{Float64},PbestC::Float64,cost::Vector{Float64})
  loc_best=findall(x -> x==minimum(cost),cost)[1];
  Gbest=Ipop[:,loc_best];
  GbestC=cost[loc_best];

  Pbest=deepcopy(Ipop);
  PbestC=deepcopy(cost);
end
@btime evaluator2!(loc_best::Int64,Gbest::Vector{Float64},GbestC::Float64,Pbest::Vector{Float64},PbestC::Float64,cost::Vector{Float64})


##
=#
