function discretize(prob::OneDimensionWE, discretization::FiniteDifferences, order::Order2)
    wave! = function (du,u,p,t)
        k = (p[1][1]/p[2][1])^2
        n = Int(floor(length(u)/2))
        @assert p[3][1] == n-1
        for i = 2:n-1
            du[i] = u[n+i]
        end
        for i = 2:n-1
            du[n+i] = k*(u[i+1]+u[i-1]-2*u[i]) + p[4][i] # Convolution and the source term.
        end
    end
    return wave!
end

function discretize(prob::TwoDimensionWE, discretization::FiniteDifferences, order::Order2)
    wave! = function (du,u,p,t)
        kx = (p[1][1]/p[2][1])^2
        ky = (p[1][2]/p[2][2])^2
        nx = Int(p[3][1])+1; ny = Int(p[3][2])+1;
        n = Int((nx+1)*(ny+1))
        # Since this is assignment with a shift, this should be fine.
        for i = 2:nx-1
            for j = 2:ny-1
                k = Int((i-1)*ny + j)
                du[k] = u[n+k]
            end
        end
        # On the other hand, the following code is going to slow down the computation.
        # We need a better ordering and a closed form formula to compute the values of the nearest neighbours.
        for i = 2:nx-1
            for j = 2:ny-1
                k = Int((i-1)*ny + j)
                @show k,i,j
                du[n+k] = kx*(u[k+1]+u[k-1]-2*u[k]) + ky*(u[k-ny]+u[k+ny]-2*u[k]) + p[4][k]
            end
        end
    end
    return wave!
end

function discretize(prob::ThreeDimensionWE, discretization::FiniteDifferences, order::Order2)
    wave! = function (du,u,p,t)

    end
end
