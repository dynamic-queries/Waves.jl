using OrdinaryDiffEq

function initialize(sim::WaveProblem,prob::OneDimensionWE)
    @unpack source_function, forcing_function, icbc_function = sim
    @unpack c,tspan, tdomain = sim
    @unpack domain_min, domain_max, n = sim

    x = Vector(LinRange(domain_min[1],domain_max[1],n[1]+1))
    dx = [x[2] - x[1]]
    source = source_function(x)
    forcing_function = discretize(OneDimensionWE(),FiniteDifferences(),Order2())
    u0 = icbc_function(x)
    p = [c,dx,n,source]
    ode_prob = ODEProblem(forcing_function,u0,tspan,p)
    return WaveSolution(sim,x,tdomain,source,ode_prob)
end

const dimension_map = Dict(1=>OneDimensionWE(),2=>TwoDimensionWE,3=>ThreeDimensionWE)

function solve(problem::WaveProblem)
    domain = dimension_map[problem.ndims]
    sol = initialize(problem,domain)
    sol.solution = OrdinaryDiffEq.solve(sol.odeprob,Tsit5(),saveat=sol.t)
    sol
end
