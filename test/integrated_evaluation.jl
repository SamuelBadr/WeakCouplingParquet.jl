@testitem "phi integrated 1D" setup = [Parameters1D] begin
    abstol = 1e-5
    reltol = 1e-5

    val_d, err_d = WeakCouplingParquet.phi2_d(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)
    val_m, err_m = WeakCouplingParquet.phi2_m(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)
    val_s, err_s = WeakCouplingParquet.phi2_s(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)
    val_t, err_t = WeakCouplingParquet.phi2_t(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)

    for part in (real, imag)
        @test abs(part(err_d)) < max(abstol, reltol * abs(part(val_d)))
        @test abs(part(err_m)) < max(abstol, reltol * abs(part(val_m)))
        @test abs(part(err_s)) < max(abstol, reltol * abs(part(val_s)))
        @test abs(part(err_t)) < max(abstol, reltol * abs(part(val_t)))
    end
    # @show val_d val_m val_s val_t
end

@testitem "full integrated 1D" setup = [Parameters1D] begin
    abstol = 1e-5
    reltol = 1e-5

    val_d, err_d = WeakCouplingParquet.full2_d(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)
    val_m, err_m = WeakCouplingParquet.full2_m(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)
    val_s, err_s = WeakCouplingParquet.full2_s(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)
    val_t, err_t = WeakCouplingParquet.full2_t(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)

    for part in (real, imag)
        @test abs(part(err_d)) < max(abstol, reltol * abs(part(val_d)))
        @test abs(part(err_m)) < max(abstol, reltol * abs(part(val_m)))
        @test abs(part(err_s)) < max(abstol, reltol * abs(part(val_s)))
        @test abs(part(err_t)) < max(abstol, reltol * abs(part(val_t)))
    end
    # @show val_d val_m val_s val_t
end

@testitem "gamma integrated 1D" setup = [Parameters1D] begin
    abstol = 1e-5
    reltol = 1e-5

    val_d, err_d = WeakCouplingParquet.gamma2_d(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)
    val_m, err_m = WeakCouplingParquet.gamma2_m(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)
    val_s, err_s = WeakCouplingParquet.gamma2_s(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)
    val_t, err_t = WeakCouplingParquet.gamma2_t(u, mu, beta, v, vp, w, k, kp, q; abstol, reltol)

    for part in (real, imag)
        @test abs(part(err_d)) < max(abstol, reltol * abs(part(val_d)))
        @test abs(part(err_m)) < max(abstol, reltol * abs(part(val_m)))
        @test abs(part(err_s)) < max(abstol, reltol * abs(part(val_s)))
        @test abs(part(err_t)) < max(abstol, reltol * abs(part(val_t)))
    end
    # @show val_d val_m val_s val_t
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

#     @time phi2_d = WeakCouplingParquet.phi2_d(u, mu, beta, v, vp, w, k, kp, q, CubatureJLh(); reltol)
#     @time phi2_m = WeakCouplingParquet.phi2_m(u, mu, beta, v, vp, w, k, kp, q, CubatureJLh(); reltol)
#     @time phi2_s = WeakCouplingParquet.phi2_s(u, mu, beta, v, vp, w, k, kp, q, CubatureJLh(); reltol)
#     @time phi2_t = WeakCouplingParquet.phi2_t(u, mu, beta, v, vp, w, k, kp, q, CubatureJLh(); reltol)

#     @show phi2_d phi2_m phi2_s phi2_t
# end
