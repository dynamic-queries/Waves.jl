using OrdinaryDiffEq
using Plots

function wave!(du,u,p,t)
    k = (p[1]/p[2])^2
    n = Int(length(u)/2)
    for i = 2:n-1
        du[i] = u[n+i]
    end
    for i = 2:n-1
        du[n+i] = k*(u[i+1]+u[i-1]-2*u[i])
    end
    # Maintain the boundary conditions

end

function IC_BC(res,x)
    u = zeros(res)
    v = zeros(res)
    # u[1:Int(floor(res/3))] = sin.(1:Int(floor(res/3)))
    u[Int(floor(end/2))] = 1.0
    u[1] = 0; u[end] = 0;
    return vcat(u,v)
end

function munge(sol)
    M = zeros(length(sol),length(sol[1]))
    for i = 1:length(sol)
        M[i,:] = sol[i]
    end
    return M
end


begin
    c = 1
    nx = 1000
    x = LinRange(0,30,nx+1)
    dx = x[2]-x[1]
    tspan = (0,8)
    t = 0:0.01:3
    p = [c,dx]
    u0 = IC_BC(nx+1,x)

    # Check the domain
    unit_map = ones(length(x))
    scatter(x,unit_map)

    # Check the initial condition
    plot(x,u0[1:Int(end/2)])

    # Define the ODE problem
    prob = ODEProblem(wave!,u0,tspan,p)
    sol = solve(prob,Tsit5(),saveat=t)

    # Reconstruct the data
    solution = munge(sol)
end

# Visualization
begin
    anim = @animate for i = 1:size(solution,1)
        plot(x,solution[i,1:Int(end/2)],ylim=[-1,1])
    end
    gif(anim,"delta.gif",fps=50)
end
