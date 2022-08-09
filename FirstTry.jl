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
GbestC=0.0;
Pbest=similar(Ipop);
PbestC=similar(cost);

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

j=1
function VelUpdate(W::Float64, Vel::Matrix{Float64}, C1::Int64, C2::Int64, R1::Vector{Float64}, R2::Vector{Float64}, Pbest::Matrix{Float64}, Gbest::Vector{Float64}, Ipop::Matrix{Float64})
  for j in 1:n
    Vel[j,i] = W*Vel[j,i] + C1 * R1[j]*(Pbest[j,i] - Ipop[j,i]) + C2 * R2[j] .* (Gbest[j] - Ipop[j,i])
  end
end
@btime VelUpdate(W::Float64, Vel::Matrix{Float64}, C1::Int64, C2::Int64, R1::Vector{Float64}, R2::Vector{Float64}, Pbest::Matrix{Float64}, Gbest::Vector{Float64}, Ipop::Matrix{Float64})


@btime Vel[:,i] = w*Vel[:,i] + C1 * R1.*(Pbest[:,i] - Ipop[:,i]) + C2 * R2 .* (Gbest - Ipop[:,i])



##
PbestC=cost
Pbest=Ipop


#aaa=@btime minimum(cost)
#aa = @btime findfirst(x -> x==aaa,cost)[1];

GbestC = minimum(cost)
loc_best=findfirst(x -> x==minimum(cost),cost)[1];
Gbest=Ipop[:,loc_best];

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
