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

@testitem "gamma integrated 2D" setup = [Parameters2D] begin
    abstol = 1e-2
    reltol = 1e-2

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
end

@testitem "no nan value in gamma when v == vp and k != kp" begin
    using StaticArrays

    u = 0.3
    beta = 0.5
    val_d, err_d = WeakCouplingParquet.gamma2_d(u, u / 2, beta, 17 * im * pi / beta, 17 * im * pi / beta, -12 * im * pi / beta, SA[pi/2, pi/2], SA[0.0, 0.0], SA[0.0, 0.0])

    @test !isnan(val_d)
    @test !isnan(err_d)
end

@testitem "no nan value in phi when v == vp" begin
    using StaticArrays

    u = 0.3
    beta = 0.5
    val_d = WeakCouplingParquet.phi2_m(u, u / 2, beta, -13 * im * pi / beta, -13 * im * pi / beta, 30 * im * pi / beta, SA[0π/4, 2π/4], SA[2π/4, 0π/4], SA[4π/4, 6π/4], SA[1.9486832980505138, 0.05131670194948623], SA[1.0, 1.0], SA[1.0, 1.0])

    @test !isnan(val_d)
end

@testitem "no nan value in integrated phi when v == vp" begin
    using StaticArrays

    u = 0.3
    beta = 0.5
    abstol = 1e-5
    reltol = 1e-5
    val_d, err_d = WeakCouplingParquet.phi2_m(u, u / 2, beta, -13 * im * pi / beta, -13 * im * pi / beta, 30 * im * pi / beta, SA[0π/4, 2π/4], SA[2π/4, 0π/4], SA[4π/4, 6π/4]; abstol, reltol)

    @test !isnan(val_d)
end
