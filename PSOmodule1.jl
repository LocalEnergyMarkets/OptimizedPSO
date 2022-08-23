__precompile__()

module PSO1

export InitPopMatrix!, CostFunc, VelUpdate!, VelMaxCheck!, VelMinCheck!, IpopUpdate!, IpopMaxCheck!, IpopMinCheck!, GlobalBestC!, GlobalBest!, LocBest!, PbestUpdate!, GbestUpdate!, PSOalgorithm!

function InitPopMatrix!(cost::Vector{Float64},Ipop::Matrix{Float64},Vel::Matrix{Float64},PbestC::Vector{Float64}, Pbest::Matrix{Float64},n::Int64,N::Int64,xmin::Float64,xmax::Float64,vmax::Float64)
  for i=1:N
    for j=1:n
      Ipop[j,i]=rand()*(xmax-xmin)+xmin;
      Vel[j,i]=rand()*(2*vmax)-vmax;
    end
    cost[i]=CostFunc(i::Int64, Ipop::Matrix{Float64}, n::Int64)
    PbestUpdate!(PbestC::Vector{Float64}, cost::Vector{Float64}, Pbest::Matrix{Float64}, Ipop::Matrix{Float64},i::Int64, n::Int64)
  end
  #return Ipop, Vel
end


function CostFunc(i::Int64, Ipop::Matrix{Float64}, n::Int64)
  c = 0.0
  Ip_i=@view Ipop[:,i]
  for k in 1:n
    c+=exp(-Ip_i[k])*sin(2*pi*Ip_i[k])
  end      # Ipop_i=Ipop[:,i]
  c
end


#=function VelUpdate!(Vel::Matrix{Float64}, n::Int64, i::Int64,W::Float64, C1::Int64, C2::Int64, R1::Vector{Float64}, R2::Vector{Float64}, Pbest::Matrix{Float64}, Gbest::Vector{Float64}, Ipop::Matrix{Float64})
  for j=1:n
    R1[j] = rand()
    R2[j] = rand()
  end
  @. Vel[:,i] = W* (@view Vel[:,i]) + C1 * R1*((@view Pbest[:,i]) - (@view Ipop[:,i])) + C2 * R2 * (Gbest - (@view Ipop[:,i]))
end=#

function VelUpdate!(Vel::Matrix{Float64}, i::Int64, W::Float64, C1::Int64, C2::Int64, Pbest::Matrix{Float64}, Gbest::Vector{Float64}, Ipop::Matrix{Float64})
  #for j=1:n
    #R1[j] = rand()
    #R2[j] = rand()
  #end
  @. Vel[:,i] = W* (@view Vel[:,i]) + C1 * rand()*((@view Pbest[:,i]) - (@view Ipop[:,i])) + C2 * rand() * (Gbest - (@view Ipop[:,i]))
end

function VelMaxCheck!(Vel::Matrix{Float64},vmax::Float64,i::Int64,n::Int64)
  V = @view Vel[:,i]
  for j=1:n
    if V[j]>vmax
      V[j]=vmax
    end
  end
end

function VelMinCheck!(Vel::Matrix{Float64},vmax::Float64,i::Int64,n::Int64)
  V = @view Vel[:,i]
  for j=1:n
    if V[j]<-vmax
      V[j]=-vmax
    end
  end
end


function IpopUpdate!(Ipop::Matrix{Float64}, Vel::Matrix{Float64}, i::Int64)
  @. Ipop[:,i] = (@view Ipop[:,i]) + (@view Vel[:,i])
end

function IpopMaxCheck!(Ipop::Matrix{Float64},xmax::Float64,i::Int64, n::Int64)
  I = @view Ipop[:,i]
  for j=1:n
    if I[j]>xmax
      I[j]=xmax
    end
    #V[j] = minimum([V[j],vmax])
  end
end

function IpopMinCheck!(Ipop::Matrix{Float64},xmin::Float64,i::Int64, n::Int64)
  I = @view Ipop[:,i]
  for j=1:n
    if I[j]<xmin
      I[j]=xmin
    end
    #V[j] = minimum([V[j],vmax])
  end
end

function PbestCalculator!(PbestC::Vector{Float64}, cost::Vector{Float64}, Pbest::Matrix{Float64}, Ip::Vector{Float64},i::Int64)
  if PbestC[i]>=cost[i]
    PbestC[i] = cost[i]
    Pbest[:,i] = Ip
  end
  return Pbest
end

function GlobalBestC!(GbestC,cost)
  GbestC=minimum(cost)
end

function LocBest!(loc_best::Int64, cost::Vector{Float64}, GbestC::Float64)
  loc_best = findfirst(x -> x==GbestC,cost)[1];
  return loc_best
end

function GlobalBest!(Gbest::Vector{Float64}, loc_best::Int64, Ipop::Matrix{Float64}, n::Int64)
  for j=1:n
    Gbest[j]=Ipop[j,loc_best]
  end
end

function PbestUpdate!(PbestC::Vector{Float64}, cost::Vector{Float64}, Pbest::Matrix{Float64}, Ipop::Matrix{Float64}, i::Int64, n::Int64)
  if PbestC[i]>=cost[i]
    PbestC[i] = cost[i]
    for j in 1:n
      Pbest[j,i] = Ipop[j,i]
    end
  end
#  return Pbest
end

function GbestUpdate!(GbestC::Float64, Gbest::Vector{Float64}, PbestC::Vector{Float64}, Pbest::Matrix{Float64}, i::Int64, n::Int64)
  if GbestC>=PbestC[i]
    GbestC = PbestC[i]
    for j in 1:n
      Gbest[j] = Pbest[j,i]
    end
  end
#  return Pbest
end

function GbestCUpdate!(GbestC::Float64, PbestC::Vector{Float64},i::Int64)
  if GbestC > PbestC[i]
    #println("GbestC")
    #println(GbestC)
    #println("PbestC[i]")
    #println(PbestC[i])
    GbestC = PbestC[i]
  else
    GbestC = GbestC
  end
end

function PSOalgorithm!(Gbest::Vector{Float64}, GbestC::Float64, Vel::Matrix{Float64}, n::Int64, N::Int64, C1::Int64, C2::Int64, R1::Vector{Float64}, R2::Vector{Float64}, Pbest::Matrix{Float64}, PbestC::Vector{Float64}, Ipop::Matrix{Float64},vmax::Float64,xmax::Float64,xmin::Float64,cost::Vector{Float64})
  for itr in 1:100
      W=rand()*(1-0.4)+0.4;
      for i in 1:N
          #VelUpdate!(Vel::Matrix{Float64}, n::Int64, i::Int64,W::Float64, C1::Int64, C2::Int64, R1::Vector{Float64}, R2::Vector{Float64}, Pbest::Matrix{Float64}, Gbest::Vector{Float64}, Ipop::Matrix{Float64})
          VelUpdate!(Vel::Matrix{Float64}, i::Int64, W::Float64, C1::Int64, C2::Int64, Pbest::Matrix{Float64}, Gbest::Vector{Float64}, Ipop::Matrix{Float64})
          VelMaxCheck!(Vel::Matrix{Float64},vmax::Float64,i::Int64,n::Int64)
          VelMinCheck!(Vel::Matrix{Float64},vmax::Float64,i::Int64,n::Int64)

          IpopUpdate!(Ipop::Matrix{Float64}, Vel::Matrix{Float64}, i::Int64)
          IpopMaxCheck!(Ipop::Matrix{Float64},xmax::Float64,i::Int64, n::Int64)
          IpopMinCheck!(Ipop::Matrix{Float64},xmin::Float64,i::Int64, n::Int64)

          cost[i]=CostFunc(i::Int64, Ipop::Matrix{Float64}, n::Int64)

          PbestUpdate!(PbestC::Vector{Float64}, cost::Vector{Float64}, Pbest::Matrix{Float64}, Ipop::Matrix{Float64}, i::Int64, n::Int64)

          GbestC = GbestCUpdate!(GbestC::Float64, PbestC::Vector{Float64}, i::Int64)
#println(GbestC)
          GbestUpdate!(GbestC::Float64, Gbest::Vector{Float64}, PbestC::Vector{Float64}, Pbest::Matrix{Float64}, i::Int64, n::Int64)
      end
  end
  return GbestC, Gbest
end



end  # modulePSO
