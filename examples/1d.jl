using Plots
using OrdinaryDiffEq

# Different types of discretization
abstract type AbstractDiscretization end
struct FiniteDifference <: AbstractDiscretization end
struct PseudoSpectral <: AbstractDiscretization end
struct FiniteElement <: AbstractDiscretization end


# Different types of operators
abstract type AbstractDiffEqOperator end
struct Δ <: AbstractDiffEqOperator
    dims::Int
end
function Δ(dims::Int)
    Δ(dims)
end
function sample!(op::Δ,u::Array)
    if op.dims == 1

    elseif op.dims == 2

    elseif op.dims == 3

    else
        ArgumentError("Laplacian operator must have dimensions ϵ (1,2,3)")
    end
end

# Different types of wave equations
# Modify this to a more generic setup using symbolics
abstract type AbstractWaveEquation end
struct OneDimensionWE <: AbstractWaveEquation end

# PDE Problem interface
# TODO: Include the stochastic wave simulation part from Felix here.
abstract type AbstractPDEProblem end
mutable struct PDEProblem <: AbstractPDEProblem
    # Operators, time and space , an array of custom op objects
    DiffOp
    # Initial values , an array of tensor values depending on the degree of the timeop
    u0
    is_periodic
    # Source terms , an array of functions
    s
    # Discretization knobs
    discretization
    time_integrator
    is_adaptive
    # Discrete versions of ICs and Source
    u0_d
    s_d
end
function PDEProblem(
    eqn::OneDimensionWE,
    u,
    v,
    s,
    discretization::AbstractDiscretization,
    isperiodic::Bool=true,
    isadaptive::Bool=true,
    FT::Type=Float64)

    # Second order integration
    Dop = Δ(1)
    time_integrator = Tsit5() # Nystrom Integrators ?

    # TODO : Create a grid here. Define a grid. Define a length function for the Grid
    # TODO : Define the sample function using SIMD intrinsics.
    grid = Grid(Dop.dims)
    u0 = Array{Array}(Array{FT}(undef,length(grid)),2)
    s_d = Array{FT}(undef,length(grid))
    sample!(u0[1],grid,u)
    sample!(u0[2],grid,v)
    sample!(s_d,grid,s)

    PDEProblem(Dop,[u,v],isperiodic,s,discretization,time_integrator,isadaptive,u0,s_d)
end


# PDE Solution type.
# t,u are the output from the adaptive time stepping scheme and the solution to the IVP at that point.
# TODO: Interpolate the field values at those points.
abstract type AbstractPDESolution end
mutable struct PDESolution <: AbstractPDESolution
    #PDE Problem
    prob::PDEProblem
    # Solution space
    t::Array
    u::Array
    # Results of interpolation - Let's do Hermite interpolation for fun.
    tsaveat::Array
    usaveat::Array
end

function __solve(prob::PDEProblem)

end
