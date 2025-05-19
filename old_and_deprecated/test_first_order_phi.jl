@testitem "phi1 ph" begin
    using StaticArrays

    u = 1.0
    mu = u / 2
    beta = 2.0
    w = im * 4 * pi / beta
    q = SVector(pi / 2, pi)
    k1 = SVector(0.5, 0.5)
    @test isfinite(WeakCouplingParquet.phi1_ph(u, mu, beta, w, q, k1))
end

@testitem "phi1 ph pole" begin
    using StaticArrays

    u = 1.0
    mu = u / 2
    beta = 2.0
    w = im * 0 * pi / beta
    q = SVector(0.0, 0.0)
    k1 = SVector(0.5, 0.5)
    @test isfinite(WeakCouplingParquet.phi1_ph(u, mu, beta, w, q, k1))
end

@testitem "phi1 s" begin
    using StaticArrays

    u = 1.0
    mu = u / 2
    beta = 2.0
    w = im * 4 * pi / beta
    q = SVector(pi / 2, pi)
    k1 = SVector(0.5, 0.5)
    @test isfinite(WeakCouplingParquet.phi1_s(u, mu, beta, w, q, k1))
end

@testitem "phi1 t" begin
    using StaticArrays

    u = 1.0
    mu = u / 2
    beta = 2.0
    w = im * 0 * pi / beta
    q = SVector(0.0, 0.0)
    k1 = SVector(0.5, 0.5)
    @test iszero(WeakCouplingParquet.phi1_t(u, mu, beta, w, q, k1))
end
