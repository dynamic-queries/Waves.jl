using Waves
using Test

"""
    Ensure that a dummy input does the laplacian computation.
"""
@testset "FD_1D" begin
    relaxer = discretize(OneDimensionWE(),FiniteDifferences(),Order2())
    u = [4.0,6.0,4.0,4.0,6.0,4.0]
    du = zero(u)
    p = [1,1,zero(u)]
    t = 0
    relaxer(du,u,p,t)
    @test du[1] == 0
    @test du[2] == 6.0
    @test du[3] == 0
    @test du[4] == 0
    @test du[5] == -4.0
    @test du[6] == 0
end


"""
    Ensure that a dummy input does the laplacian computation in 2D.
"""
@testset "FD_2D" begin
    relaxer = discretize(TwoDimensionWE(),FiniteDifferences(),Order2())
    u = [1.0,4.0,7.0,2.0,5.0,9.0,3.0,6.0,9.0,1.0,4.0,7.0,2.0,5.0,9.0,3.0,6.0,9.0]
    du = zero(u)
    p = [[1,1],[1,1],[3,3],zero(u)]
    @assert length(u) == p[3][1]*p[3][2]*2
    t = 0
    relaxer(du,u,p,t)

    du = reshape(du,(Int(p[3][1]),Int(p[3][2]),2))
    @test du[1] ≈ 0
    @test du[2] ≈ 0
    @test du[3] ≈ 0
    @test du[4] ≈ 0
    @test du[5] ≈ 5
    @test du[6] ≈ 0
    @test du[7] ≈ 0
    @test du[8] ≈ 0
    @test du[9] ≈ 0

    @test du[10] ≈ 0
    @test du[11] ≈ 0
    @test du[12] ≈ 0
    @test du[13] ≈ 0
    @test du[14] ≈ 1
    @test du[15] ≈ 0
    @test du[16] ≈ 0
    @test du[17] ≈ 0
    @test du[18] ≈ 0
end
