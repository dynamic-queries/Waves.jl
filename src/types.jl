using Base
using OrdinaryDiffEq

# Pre defined models
abstract type AbstractWaveEquation end
struct OneDimensionWE <: AbstractWaveEquation end
struct TwoDimensionWE <: AbstractWaveEquation end
struct ThreeDimensionWE <: AbstractWaveEquation end

# Types of space discretization
abstract type AbstractDiscretization end
struct FiniteDifferences <: AbstractDiscretization end
struct Spectral <: AbstractDiscretization end

# Order of the discretization operator
abstract type AbstractOrder end
struct Order2 <: AbstractOrder end

abstract type AbstractWaveProblem end
abstract type AbstractWaveSolution end

mutable struct WaveProblem{FT,IT} <:AbstractWaveProblem
    ndims::IT
    # Functions for simulation
    forcing_function
    source_function
    icbc_function
    # Simulation arguments
    c::Vector{FT}
    tspan::Tuple

    # Domain parameters
    domain_min::Vector{FT}
    domain_max::Vector{FT}
    n::Vector{IT}

    tdomain::Vector{FT}
    # Constructor
    WaveProblem(f,source,icbc,c,n,Lmin,Lmax,tspan,t) = new{eltype(c),eltype(n)}(length(n),f,source,icbc,c,tspan,Lmin,Lmax,n,t)
end

mutable struct WaveSolution <: AbstractWaveSolution
    prob::WaveProblem
    x::Array
    t::Array
    source::Array
    odeprob::ODEProblem
    solution::ODESolution
    WaveSolution(prob,x,tdomain,source,odeprob) = new(prob,x,tdomain,source,odeprob)
end
