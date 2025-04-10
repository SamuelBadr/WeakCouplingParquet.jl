using Test
import WeakCouplingParquet as WCP
using StaticArrays
using Integrals

@testset "evaluation" begin
    beta = 10.0
    u = 3.0
    mu = u / 2

    n = 1
    np = 11
    m = -4

    v = im * n * pi / beta
    vp = im * np * pi / beta
    w = im * m * pi / beta

    k = SA[1.2, 1.9]
    kp = SA[0.5, 1.4]
    q = SA[0.3, 0.8]

    k1 = SA[0.1, 0.9]
    k2 = SA[0.2, 1.0]
    k3 = SA[0.3, 1.1]

    @testset "phi" begin
        phi2_d = WCP.phi2_d(u, mu, beta, v, vp, w, k, kp, q, k1, k2, k3)
        phi2_m = WCP.phi2_m(u, mu, beta, v, vp, w, k, kp, q, k1, k2, k3)
        phi2_s = WCP.phi2_s(u, mu, beta, v, vp, w, k, kp, q, k1, k2, k3)
        phi2_t = WCP.phi2_t(u, mu, beta, v, vp, w, k, kp, q, k1, k2, k3)

        @test isfinite(phi2_d)
        @test isfinite(phi2_m)
        @test isfinite(phi2_s)
        @test isfinite(phi2_t)

        @inferred WCP.phi2_d(u, mu, beta, v, vp, w, k, kp, q, k1, k2, k3)
        @inferred WCP.phi2_m(u, mu, beta, v, vp, w, k, kp, q, k1, k2, k3)
        @inferred WCP.phi2_s(u, mu, beta, v, vp, w, k, kp, q, k1, k2, k3)
        @inferred WCP.phi2_t(u, mu, beta, v, vp, w, k, kp, q, k1, k2, k3)

        @test @allocated(WCP.phi2_d(u, mu, beta, v, vp, w, k, kp, q, k1, k2, k3)) == 0
        @test @allocated(WCP.phi2_m(u, mu, beta, v, vp, w, k, kp, q, k1, k2, k3)) == 0
        @test @allocated(WCP.phi2_s(u, mu, beta, v, vp, w, k, kp, q, k1, k2, k3)) == 0
        @test @allocated(WCP.phi2_t(u, mu, beta, v, vp, w, k, kp, q, k1, k2, k3)) == 0
    end

    @testset "full" begin
        full2_d = WCP.full2_d(u, mu, beta, v, vp, w, k, kp, q, k1)
        full2_m = WCP.full2_m(u, mu, beta, v, vp, w, k, kp, q, k1)
        full2_s = WCP.full2_s(u, mu, beta, v, vp, w, k, kp, q, k1)
        full2_t = WCP.full2_t(u, mu, beta, v, vp, w, k, kp, q, k1)

        @test isfinite(full2_d)
        @test isfinite(full2_m)
        @test isfinite(full2_s)
        @test isfinite(full2_t)

        @inferred WCP.full2_d(u, mu, beta, v, vp, w, k, kp, q, k1)
        @inferred WCP.full2_m(u, mu, beta, v, vp, w, k, kp, q, k1)
        @inferred WCP.full2_s(u, mu, beta, v, vp, w, k, kp, q, k1)
        @inferred WCP.full2_t(u, mu, beta, v, vp, w, k, kp, q, k1)

        @test @allocated(WCP.full2_d(u, mu, beta, v, vp, w, k, kp, q, k1)) == 0
        @test @allocated(WCP.full2_m(u, mu, beta, v, vp, w, k, kp, q, k1)) == 0
        @test @allocated(WCP.full2_s(u, mu, beta, v, vp, w, k, kp, q, k1)) == 0
        @test @allocated(WCP.full2_t(u, mu, beta, v, vp, w, k, kp, q, k1)) == 0
    end

    @testset "gamma" begin
        gamma2_d = WCP.gamma2_d(u, mu, beta, v, vp, w, k, kp, q, k1)
        gamma2_m = WCP.gamma2_m(u, mu, beta, v, vp, w, k, kp, q, k1)
        gamma2_s = WCP.gamma2_s(u, mu, beta, v, vp, w, k, kp, q, k1)
        gamma2_t = WCP.gamma2_t(u, mu, beta, v, vp, w, k, kp, q, k1)

        @test isfinite(gamma2_d)
        @test isfinite(gamma2_m)
        @test isfinite(gamma2_s)
        @test isfinite(gamma2_t)

        @inferred WCP.gamma2_d(u, mu, beta, v, vp, w, k, kp, q, k1)
        @inferred WCP.gamma2_m(u, mu, beta, v, vp, w, k, kp, q, k1)
        @inferred WCP.gamma2_s(u, mu, beta, v, vp, w, k, kp, q, k1)
        @inferred WCP.gamma2_t(u, mu, beta, v, vp, w, k, kp, q, k1)

        @test @allocated(WCP.gamma2_d(u, mu, beta, v, vp, w, k, kp, q, k1)) == 0
        @test @allocated(WCP.gamma2_m(u, mu, beta, v, vp, w, k, kp, q, k1)) == 0
        @test @allocated(WCP.gamma2_s(u, mu, beta, v, vp, w, k, kp, q, k1)) == 0
        @test @allocated(WCP.gamma2_t(u, mu, beta, v, vp, w, k, kp, q, k1)) == 0
    end
end

@testset "1d - integrated" begin
    beta = 1.0
    u = 10.0
    mu = u / 2

    n = 1
    np = 11
    m = -4

    v = im * n * pi / beta
    vp = im * np * pi / beta
    w = im * m * pi / beta

    k = SA[1.2]
    kp = SA[0.5]
    q = SA[0.3]

    abstol = 1e-5
    reltol = 1e-5

    @testset "phi" begin
        val_d, err_d = WCP.phi2_d(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)
        val_m, err_m = WCP.phi2_m(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)
        val_s, err_s = WCP.phi2_s(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)
        val_t, err_t = WCP.phi2_t(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)

        for part in (real, imag)
            @test abs(part(err_d)) < max(abstol, reltol * abs(part(val_d)))
            @test abs(part(err_m)) < max(abstol, reltol * abs(part(val_m)))
            @test abs(part(err_s)) < max(abstol, reltol * abs(part(val_s)))
            @test abs(part(err_t)) < max(abstol, reltol * abs(part(val_t)))
        end
        # @show val_d val_m val_s val_t
    end

    @testset "full" begin
        val_d, err_d = WCP.full2_d(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)
        val_m, err_m = WCP.full2_m(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)
        val_s, err_s = WCP.full2_s(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)
        val_t, err_t = WCP.full2_t(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)

        for part in (real, imag)
            @test abs(part(err_d)) < max(abstol, reltol * abs(part(val_d)))
            @test abs(part(err_m)) < max(abstol, reltol * abs(part(val_m)))
            @test abs(part(err_s)) < max(abstol, reltol * abs(part(val_s)))
            @test abs(part(err_t)) < max(abstol, reltol * abs(part(val_t)))
        end
        # @show val_d val_m val_s val_t
    end

    @testset "gamma" begin
        val_d, err_d = WCP.gamma2_d(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)
        val_m, err_m = WCP.gamma2_m(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)
        val_s, err_s = WCP.gamma2_s(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)
        val_t, err_t = WCP.gamma2_t(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)

        for part in (real, imag)
            @test abs(part(err_d)) < max(abstol, reltol * abs(part(val_d)))
            @test abs(part(err_m)) < max(abstol, reltol * abs(part(val_m)))
            @test abs(part(err_s)) < max(abstol, reltol * abs(part(val_s)))
            @test abs(part(err_t)) < max(abstol, reltol * abs(part(val_t)))
        end
        # @show val_d val_m val_s val_t
    end
end


# this is currently pretty slow. not terribly useful

# @testset "2d - integrated" begin
#     beta = 1.0
#     u = 1.0
#     mu = u / 2# + im * 0.01

#     n = 1
#     np = 11
#     m = -4

#     v = im * n * pi / beta
#     vp = im * np * pi / beta
#     w = im * m * pi / beta

#     k = SA[1.2, 0.9]
#     kp = SA[0.5, 0.9]
#     q = SA[0.3, 0.9]

#     reltol = 1e-2

#     @time phi2_d = WCP.phi2_d(u, mu, beta, v, vp, w, k, kp, q, CubatureJLh(); reltol)
#     @time phi2_m = WCP.phi2_m(u, mu, beta, v, vp, w, k, kp, q, CubatureJLh(); reltol)
#     @time phi2_s = WCP.phi2_s(u, mu, beta, v, vp, w, k, kp, q, CubatureJLh(); reltol)
#     @time phi2_t = WCP.phi2_t(u, mu, beta, v, vp, w, k, kp, q, CubatureJLh(); reltol)

#     @show phi2_d phi2_m phi2_s phi2_t
# end
