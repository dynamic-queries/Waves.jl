function visualize()

end

function animate(prob::OneDimensionWE, sim::DeterministicSimulation, filename::String)
    anim = @animate for i = 1:size(sim.solution,1)
        plot(sim.x,sim.solution[i,1:Int(end/2)],ylim=[-1,1])
    end
    gif(anim,filename,fps=50)
end
