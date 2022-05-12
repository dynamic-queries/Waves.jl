include("types.jl")
include("utils.jl")
using OrdinaryDiffEq
using Plots

mutable struct Simulation{FT,IT}
    c::FT
    L::FT
    nx::IT
    tspan
    x
    dx
    p
    u0
    solution
    Simulation(c,nx,L,tspan) = new{typeof(c),typeof(nx)}(c,L,nx,tspan)
end

function simulate!(prob::OneDimensionWE, sim::Simulation)

    sim.x = LinRange(0,sim.L,sim.nx+1)
    sim.dx = sim.x[2]-sim.x[1]
    sim.p = [sim.c,sim.dx]
    sim.u0 = IC_BC(sim.nx+1,sim.x)

    # Define the ODE problem
    prob = ODEProblem(wave!,u0,tspan,p)
    sol = solve(prob,Tsit5())

    # Reconstruct the data
    solution = munge(sol)
    sim.solution = solution
end

function visualize(prob::OneDimensionWE, sim::Simulation, filename::String)
    anim = @animate for i = 1:size(sim.solution,1)
        plot(sim.x,sim.solution[i,1:Int(end/2)],ylim=[-1,1])
    end
    gif(anim,filename,fps=50)
end

#TODO : Implement more 1D examples
#TODO : Implement a 2D prototype and test even more examples
#TODO : See if xtalx has implemented the schemes that we spoke of
#TODO : Implement the 2D and 3D version of the Fourier Neural operator.
#TODO : Implement examples of neural networks that could be of use.
#TODO : Implement the 1D FNO
#TODO : Add more examples and establish means to test the fidelity of these models.
