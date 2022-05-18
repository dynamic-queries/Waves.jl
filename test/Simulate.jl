using Waves
using Waves:discretize,solve
using Test
using Plots
using OrdinaryDiffEq

@testset "1D-WE" begin
    #=
        Defining a Wave Equation Problem involves:
            i.   Define the source function for your domain - source is assumed independent of time.(Suffices for my application)
            ii.  Define the ICBC function which assigns the initial and boundary condition for the field F over your domain.
            iii. Define a forcing function for your discretization method for the space domain.
            iv.  Define the parameters of the simulation
                    - c : wave velocities
                    - n : number of points
                    - Lmin : min extent of the domain
                    - Lmax : max extent of the domain
                    - tspan : time span for time integration
                    - t : time span over which the values have to be intepolated (adaptive integration)
            v.   Pass these parameters to the WaveEquation object
            vi.  Solve the problem.
            vii. Visualization is WIP.
        ----------------------------------------------------------------------------------------------------------------
    =#

    function s(x)
        u = zero(x)
        return u
    end
    #=
    ----------------------------------------------------------------------------------------------------------------
    =#
    function ICBC(x)
        res = length(x)
        u = zeros(res)
        v = zeros(res)
        a = 2.0
        x0 = x[Int(floor(end/2))]
        u .= ((x0.-x)./(sqrt(pi)*a^2)) .* exp.(-((x .- x0)/(sqrt(2)*a)).^2)
        return vcat(u,v)
    end
    #=
    ----------------------------------------------------------------------------------------------------------------
    =#
    forcing_function = discretize(OneDimensionWE(),FiniteDifferences(),Order2())
    #=
    ----------------------------------------------------------------------------------------------------------------
    =#
    c = [1.0]
    n = [1000]
    Lmin = [0.0]
    Lmax = [30.,0]
    tspan = (0,100)
    t = Vector(0:0.01:100)
    #=
    ----------------------------------------------------------------------------------------------------------------
    =#
    prob = WaveProblem(forcing_function,s,ICBC,c,n,Lmin,Lmax,tspan,t)
    sol = solve(prob)
    #=
    ----------------------------------------------------------------------------------------------------------------
    =#
    solution = Waves.munge(sol.solution)
    x = sol.x
    using Plots

    anim = @animate for i = 1:size(solution,1)
        plot(x,solution[i,1:Int(end/2)],ylim=[-1,1])
    end
    gif(anim,"./Waves/examples/1D.gif",fps=50)
end

#=
----------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------
=#

function s(X,Y)
    xmid = X[:,1][Int(floor(end/2))]
    ymid = Y[:,2][Int(floor(end/2))]
    ax = 1
    ay = 1
    return exp.(-(((X .- xmid)./2*ax).^2+((Y .- ymid)./2*ay).^2))
end

function ICBC(X,Y)
    u = zero(X)
    v = zero(Y)
    return vcat(vec(u),vec(v))
end

forcing_function = discretize(TwoDimensionWE(),FiniteDifferences(),Order2())

c = [1.0,1.0]
n = [100,100]
Lmin = [-10.0,-10.0]
Lmax = [10.,10.]
tspan = (0,10)
t = Vector(0:0.01:10)

prob = WaveProblem(forcing_function,s,ICBC,c,n,Lmin,Lmax,tspan,t)
sol = initialize(prob,TwoDimensionWE)
sol = solve(prob)
