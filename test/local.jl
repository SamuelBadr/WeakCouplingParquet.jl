@testitem "phi2 local" setup = [Parameters0D] begin
    @test isfinite(WeakCouplingParquet.phi2_d(u, mu, beta, v, vp, w, k, kp, q))
    @test isfinite(WeakCouplingParquet.phi2_m(u, mu, beta, v, vp, w, k, kp, q))
    @test isfinite(WeakCouplingParquet.phi2_s(u, mu, beta, v, vp, w, k, kp, q))
    @test isfinite(WeakCouplingParquet.phi2_t(u, mu, beta, v, vp, w, k, kp, q))
end

@testitem "full2 local" setup = [Parameters0D] begin
    @test isfinite(WeakCouplingParquet.full2_d(u, mu, beta, v, vp, w, k, kp, q))
    @test isfinite(WeakCouplingParquet.full2_m(u, mu, beta, v, vp, w, k, kp, q))
    @test isfinite(WeakCouplingParquet.full2_s(u, mu, beta, v, vp, w, k, kp, q))
    @test isfinite(WeakCouplingParquet.full2_t(u, mu, beta, v, vp, w, k, kp, q))
end

@testitem "gamma2 local" setup = [Parameters0D] begin
    @test isfinite(WeakCouplingParquet.gamma2_d(u, mu, beta, v, vp, w, k, kp, q))
    @test isfinite(WeakCouplingParquet.gamma2_m(u, mu, beta, v, vp, w, k, kp, q))
    @test isfinite(WeakCouplingParquet.gamma2_s(u, mu, beta, v, vp, w, k, kp, q))
    @test isfinite(WeakCouplingParquet.gamma2_t(u, mu, beta, v, vp, w, k, kp, q))
end
