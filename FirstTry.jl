using Pkg
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
using BenchmarkTools
const n=Int64(1000);
const N=Int64(10*n);
#const N=Int64(5);
##
Ipop=zeros(Float64,n,N);
Vel=zeros(Float64,n,N);
cost=zeros(Float64,N);

xmin=Float64(0);
xmax=Float64(1);
vmax=0.25*(xmax-xmin);


function InitPopMatrix(n::Int64,N::Int64,Ipop::Matrix{Float64},Vel::Matrix{Float64},xmin::Float64,xmax::Float64,vmax::Float64)
  for i=1:N
    for j=1:n
      Ipop[j,i]=rand()*(xmax-xmin)+xmin;
      Vel[j,i]=rand()*(2*vmax)-vmax;
    end
  end
  return Ipop, Vel
end
Ipop, Vel = InitPopMatrix(n::Int64,N::Int64,Ipop::Matrix{Float64},Vel::Matrix{Float64},xmin::Float64,xmax::Float64,vmax::Float64)
##
i=5
function CostFunc(i::Int64, Ipop::Matrix{Float64})
  c = 0
  Ip_i=@view Ipop[:,i]
  for k in 1:n
    c+=exp(-Ip_i[k])*sin(2*pi*Ip_i[k])
  end      # Ipop_i=Ipop[:,i]
  c
end
#cost[i] = @btime  CostFunc(i::Int64,Ipop::Matrix{Float64})

##
#@btime for i=1:N
#  cost[i] = CostFunc(i::Int64,Ipop::Matrix{Float64})
#end

function test1!(cost::Vector{Float64}, Ipop::Matrix{Float64})
  for i=1:N
    cost[i] = CostFunc(i::Int64,Ipop::Matrix{Float64})
  end
end
@time test1!(cost::Vector{Float64}, Ipop::Matrix{Float64})

function test2!(cost::Vector{Float64}, Ipop::Matrix{Float64})
  for i=1:N
    cost[i] = CostFunc(i::Int64,Ipop::Matrix{Float64})
  end
end
