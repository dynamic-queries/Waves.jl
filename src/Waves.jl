module Waves
# Write your package code here.
    using Parameters: @unpack

    include("Utils.jl")
    include("Types.jl")
    include("SemiDiscretize.jl")
    include("Simulate.jl")

    export OneDimensionWE,TwoDimensionWE,ThreeDimensionWE
    export WaveProblem, WaveSolution
    export FiniteDifferences,Spectral
    export Order2
    export discretize, initialize, solve
end
