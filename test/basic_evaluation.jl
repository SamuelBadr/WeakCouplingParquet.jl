@testitem "phi not integrated 2D" setup = [Parameters2D] begin
    phi2_d = WeakCouplingParquet.phi2_d(u, mu, beta, v, vp, w, k, kp, q, k1, k2, k3)
    phi2_m = WeakCouplingParquet.phi2_m(u, mu, beta, v, vp, w, k, kp, q, k1, k2, k3)
    phi2_s = WeakCouplingParquet.phi2_s(u, mu, beta, v, vp, w, k, kp, q, k1, k2, k3)
    phi2_t = WeakCouplingParquet.phi2_t(u, mu, beta, v, vp, w, k, kp, q, k1, k2, k3)

    @test isfinite(phi2_d)
    @test isfinite(phi2_m)
    @test isfinite(phi2_s)
    @test isfinite(phi2_t)

    @inferred WeakCouplingParquet.phi2_d(u, mu, beta, v, vp, w, k, kp, q, k1, k2, k3)
    @inferred WeakCouplingParquet.phi2_m(u, mu, beta, v, vp, w, k, kp, q, k1, k2, k3)
    @inferred WeakCouplingParquet.phi2_s(u, mu, beta, v, vp, w, k, kp, q, k1, k2, k3)
    @inferred WeakCouplingParquet.phi2_t(u, mu, beta, v, vp, w, k, kp, q, k1, k2, k3)
end

@testitem "full not integrated 2D" setup = [Parameters2D] begin
    full2_d = WeakCouplingParquet.full2_d(u, mu, beta, v, vp, w, k, kp, q, k1)
    full2_m = WeakCouplingParquet.full2_m(u, mu, beta, v, vp, w, k, kp, q, k1)
    full2_s = WeakCouplingParquet.full2_s(u, mu, beta, v, vp, w, k, kp, q, k1)
    full2_t = WeakCouplingParquet.full2_t(u, mu, beta, v, vp, w, k, kp, q, k1)

    @test isfinite(full2_d)
    @test isfinite(full2_m)
    @test isfinite(full2_s)
    @test isfinite(full2_t)

    @inferred WeakCouplingParquet.full2_d(u, mu, beta, v, vp, w, k, kp, q, k1)
    @inferred WeakCouplingParquet.full2_m(u, mu, beta, v, vp, w, k, kp, q, k1)
    @inferred WeakCouplingParquet.full2_s(u, mu, beta, v, vp, w, k, kp, q, k1)
    @inferred WeakCouplingParquet.full2_t(u, mu, beta, v, vp, w, k, kp, q, k1)
end

@testitem "gamma not integrated 2D" setup = [Parameters2D] begin
    gamma2_d = WeakCouplingParquet.gamma2_d(u, mu, beta, v, vp, w, k, kp, q, k1)
    gamma2_m = WeakCouplingParquet.gamma2_m(u, mu, beta, v, vp, w, k, kp, q, k1)
    gamma2_s = WeakCouplingParquet.gamma2_s(u, mu, beta, v, vp, w, k, kp, q, k1)
    gamma2_t = WeakCouplingParquet.gamma2_t(u, mu, beta, v, vp, w, k, kp, q, k1)

    @test isfinite(gamma2_d)
    @test isfinite(gamma2_m)
    @test isfinite(gamma2_s)
    @test isfinite(gamma2_t)

    @inferred WeakCouplingParquet.gamma2_d(u, mu, beta, v, vp, w, k, kp, q, k1)
    @inferred WeakCouplingParquet.gamma2_m(u, mu, beta, v, vp, w, k, kp, q, k1)
    @inferred WeakCouplingParquet.gamma2_s(u, mu, beta, v, vp, w, k, kp, q, k1)
    @inferred WeakCouplingParquet.gamma2_t(u, mu, beta, v, vp, w, k, kp, q, k1)
end

@testitem "phi 0D" setup = [Parameters0D] begin
    phi2_d = WeakCouplingParquet.phi2_d(u, mu, beta, v, vp, w)
    phi2_m = WeakCouplingParquet.phi2_m(u, mu, beta, v, vp, w)
    phi2_s = WeakCouplingParquet.phi2_s(u, mu, beta, v, vp, w)
    phi2_t = WeakCouplingParquet.phi2_t(u, mu, beta, v, vp, w)

    @test isfinite(phi2_d)
    @test isfinite(phi2_m)
    @test isfinite(phi2_s)
    @test isfinite(phi2_t)

    @inferred WeakCouplingParquet.phi2_d(u, mu, beta, v, vp, w)
    @inferred WeakCouplingParquet.phi2_m(u, mu, beta, v, vp, w)
    @inferred WeakCouplingParquet.phi2_s(u, mu, beta, v, vp, w)
    @inferred WeakCouplingParquet.phi2_t(u, mu, beta, v, vp, w)
end

@testitem "full 0D" setup = [Parameters0D] begin
    full2_d = WeakCouplingParquet.full2_d(u, mu, beta, v, vp, w)
    full2_m = WeakCouplingParquet.full2_m(u, mu, beta, v, vp, w)
    full2_s = WeakCouplingParquet.full2_s(u, mu, beta, v, vp, w)
    full2_t = WeakCouplingParquet.full2_t(u, mu, beta, v, vp, w)

    @test isfinite(full2_d)
    @test isfinite(full2_m)
    @test isfinite(full2_s)
    @test isfinite(full2_t)

    @inferred WeakCouplingParquet.full2_d(u, mu, beta, v, vp, w)
    @inferred WeakCouplingParquet.full2_m(u, mu, beta, v, vp, w)
    @inferred WeakCouplingParquet.full2_s(u, mu, beta, v, vp, w)
    @inferred WeakCouplingParquet.full2_t(u, mu, beta, v, vp, w)
end

@testitem "gamma 0D" setup = [Parameters0D] begin
    gamma2_d = WeakCouplingParquet.gamma2_d(u, mu, beta, v, vp, w)
    gamma2_m = WeakCouplingParquet.gamma2_m(u, mu, beta, v, vp, w)
    gamma2_s = WeakCouplingParquet.gamma2_s(u, mu, beta, v, vp, w)
    gamma2_t = WeakCouplingParquet.gamma2_t(u, mu, beta, v, vp, w)

    @test isfinite(gamma2_d)
    @test isfinite(gamma2_m)
    @test isfinite(gamma2_s)
    @test isfinite(gamma2_t)

    @inferred WeakCouplingParquet.gamma2_d(u, mu, beta, v, vp, w)
    @inferred WeakCouplingParquet.gamma2_m(u, mu, beta, v, vp, w)
    @inferred WeakCouplingParquet.gamma2_s(u, mu, beta, v, vp, w)
    @inferred WeakCouplingParquet.gamma2_t(u, mu, beta, v, vp, w)
end
