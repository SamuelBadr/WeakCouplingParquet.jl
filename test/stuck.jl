@testitem "this previously didn't work" begin
    using StaticArrays

    u = 0.5
    beta = 1.0
    mu = u / 2

    v = 7.0
    vp = 5.0
    w = 0.0

    k = SA[0.0]
    kp = SA[1.5]
    q = SA[0.0]

    abstol = 1e-1
    reltol = 1e-1
    val, err = WeakCouplingParquet.full2_m(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)
    @test isfinite(val)
end
