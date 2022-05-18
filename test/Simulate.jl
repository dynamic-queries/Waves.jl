using Waves
using Waves:discretize,solve
using Test
using OrdinaryDiffEq

begin
    function s(x)
        # a = 1
        # x0 = x[Int(floor(end/2))]
        u = zero(x)
        # u .= ((-x0.+x)./(sqrt(pi)*a^2)) .* exp.(-((x .- x0)/(sqrt(2)*a)).^2)
        return u
    end
    function ICBC(x)
        res = length(x)
        u = zeros(res)
        v = zeros(res)
        a = 2.0
        u .= cos.(x)
        # x0 = x[Int(floor(end/2))]
        # u .= ((x0.-x)./(sqrt(pi)*a^2)) .* exp.(-((x .- x0)/(sqrt(2)*a)).^2)
        return vcat(u,v)
    end
    forcing_function = discretize(OneDimensionWE(),FiniteDifferences(),Order2())
    c = [1.0]
    n = [1000]
    Lmin = [0.0]
    Lmax = [30.,0]
    tspan = (0,100)
    t = Vector(0:0.01:100)

    prob = WaveProblem(forcing_function,s,ICBC,c,n,Lmin,Lmax,tspan,t)
    sol = solve(prob)

    solution = Waves.munge(sol.solution)
    x = sol.x
    using Plots

    anim = @animate for i = 1:size(solution,1)
        plot(x,solution[i,1:Int(end/2)],ylim=[-1,1])
    end
    gif(anim,"delta.gif",fps=50)
end
