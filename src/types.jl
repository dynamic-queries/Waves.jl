using Base

abstract type AbstractWaveEquation end
struct OneDimensionWE <: AbstractWaveEquation end
struct TwoDimensionWE <: AbstractWaveEquation end
struct ThreeDimensionWE <: AbstractWaveEquation end
