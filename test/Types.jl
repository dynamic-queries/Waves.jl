using Waves
using OrdinaryDiffEq
using Test

@testset "Basic Type Existence" begin
    # Init all abt types
    one = OneDimensionWE()
    two = TwoDimensionWE()
    three = ThreeDimensionWE()

    # Init all methods
    fd = FiniteDifferences()
    spec = Spectral()

    # Orders
    order = Order2()
end


@testset "WaveProblem and Wave Solution Existence" begin
    ndims = 1
    forcing_function(du,u,p,t) = Vector(u)
    source_function(x) = x
    icbc_function(x) = x

    c = [1.0]
    tspan = (2.0,10.0)
    tdomain = Vector(2:0.1:10)
    domain_min = [0.0]
    domain_max = [1.0]
    n = [100]

    prob = WaveProblem(forcing_function,source_function,icbc_function,
                        c,n,domain_min,domain_max,tspan,tdomain)


    x = Array(LinRange(domain_min[1],domain_max[1],n[1]))
    dx = [x[2]-x[1]]
    source = sin.(x)
    u0 = zero(x)
    odeprob = ODEProblem(forcing_function,u0,tspan,[c[1],dx[1],n[1],source])

    solution = WaveSolution(prob,x,tdomain,source,odeprob)
end
