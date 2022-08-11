__precompile__()

module PSO

export InitPopMatrix!, CostFunc, CostAll!, VelUpdate, VelMaxCheck!, VelMinCheck!, VelAll!, PopAll!

function InitPopMatrix!(Ipop::Matrix{Float64},Vel::Matrix{Float64},n::Int64,N::Int64,xmin::Float64,xmax::Float64,vmax::Float64)
  for i=1:N
    for j=1:n
      Ipop[j,i]=rand()*(xmax-xmin)+xmin;
      Vel[j,i]=rand()*(2*vmax)-vmax;
    end
  end
  #return Ipop, Vel
end


function CostFunc(i::Int64, Ipop::Matrix{Float64}, n::Int64)
  c = 0
  Ip_i=@view Ipop[:,i]
  for k in 1:n
    c+=exp(-Ip_i[k])*sin(2*pi*Ip_i[k])
  end      # Ipop_i=Ipop[:,i]
  c
end

function CostAll!(cost::Vector{Float64}, Ipop::Matrix{Float64}, N::Int64, n::Int64)
  for i=1:N
    cost[i] = CostFunc(i::Int64,Ipop::Matrix{Float64}, n::Int64)
  end
end

function VelUpdate!(Vel::Matrix{Float64}, n::Int64, i::Int64,W::Float64, C1::Int64, C2::Int64, R1::Vector{Float64}, R2::Vector{Float64}, Pbest::Matrix{Float64}, Gbest::Vector{Float64}, Ipop::Matrix{Float64})
  for j=1:n
    R1[j] = rand()
    R2[j] = rand()
  end
  @. Vel[:,i] = W* (@view Vel[:,i]) + C1 * R1*((@view Pbest[:,i]) - (@view Ipop[:,i])) + C2 * R2 * (Gbest - (@view Ipop[:,i]))
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

function VelAll!(N::Int64,n::Int64, W::Float64, Vel::Matrix{Float64}, C1::Int64, C2::Int64, R1::Vector{Float64}, R2::Vector{Float64}, Pbest::Matrix{Float64}, Gbest::Vector{Float64}, Ipop::Matrix{Float64}, vmax::Float64)
  for i=1:N
    VelUpdate!(Vel::Matrix{Float64}, n::Int64, i::Int64,W::Float64, C1::Int64, C2::Int64, R1::Vector{Float64}, R2::Vector{Float64}, Pbest::Matrix{Float64}, Gbest::Vector{Float64}, Ipop::Matrix{Float64})
    VelMaxCheck!(Vel::Matrix{Float64},vmax::Float64,i::Int64,n::Int64)
    VelMinCheck!(Vel::Matrix{Float64},vmax::Float64,i::Int64,n::Int64)
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

function PopAll!(N::Int64, n::Int64, Vel::Matrix{Float64}, Ipop::Matrix{Float64}, xmax::Float64, xmin::Float64)
  for i=1:N
    IpopUpdate!(Ipop::Matrix{Float64}, Vel::Matrix{Float64}, i::Int64)
    IpopMaxCheck!(Ipop::Matrix{Float64},xmax::Float64,i::Int64, n::Int64)
    IpopMinCheck!(Ipop::Matrix{Float64},xmin::Float64,i::Int64, n::Int64)
  end
end


end  # modulePSO
