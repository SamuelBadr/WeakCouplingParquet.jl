@testitem "phi integrated 1D" setup = [Parameters1D] begin
    val_d = WeakCouplingParquet.phi2_d(u, mu, beta, v, vp, w, k, kp, q)
    val_m = WeakCouplingParquet.phi2_m(u, mu, beta, v, vp, w, k, kp, q)
    val_s = WeakCouplingParquet.phi2_s(u, mu, beta, v, vp, w, k, kp, q)
    val_t = WeakCouplingParquet.phi2_t(u, mu, beta, v, vp, w, k, kp, q)

    val_d_better = WeakCouplingParquet.phi2_d(u, mu, beta, v, vp, w, k, kp, q; n=41)
    val_m_better = WeakCouplingParquet.phi2_m(u, mu, beta, v, vp, w, k, kp, q; n=41)
    val_s_better = WeakCouplingParquet.phi2_s(u, mu, beta, v, vp, w, k, kp, q; n=41)
    val_t_better = WeakCouplingParquet.phi2_t(u, mu, beta, v, vp, w, k, kp, q; n=41)

    for part in (real, imag)
        @test abs(part(val_d_better - val_d)) / abs(val_d_better) < 1e-8
        @test abs(part(val_m_better - val_m)) / abs(val_m_better) < 1e-8
        @test abs(part(val_s_better - val_s)) / abs(val_s_better) < 1e-8
        @test abs(part(val_t_better - val_t)) / abs(val_t_better) < 1e-8
    end
end

@testitem "full integrated 1D" setup = [Parameters1D] begin
    val_d = WeakCouplingParquet.full2_d(u, mu, beta, v, vp, w, k, kp, q)
    val_m = WeakCouplingParquet.full2_m(u, mu, beta, v, vp, w, k, kp, q)
    val_s = WeakCouplingParquet.full2_s(u, mu, beta, v, vp, w, k, kp, q)
    val_t = WeakCouplingParquet.full2_t(u, mu, beta, v, vp, w, k, kp, q)

    val_d_better = WeakCouplingParquet.full2_d(u, mu, beta, v, vp, w, k, kp, q; n=41)
    val_m_better = WeakCouplingParquet.full2_m(u, mu, beta, v, vp, w, k, kp, q; n=41)
    val_s_better = WeakCouplingParquet.full2_s(u, mu, beta, v, vp, w, k, kp, q; n=41)
    val_t_better = WeakCouplingParquet.full2_t(u, mu, beta, v, vp, w, k, kp, q; n=41)

    for part in (real, imag)
        @test abs(part(val_d_better - val_d)) / abs(val_d_better) < 1e-8
        @test abs(part(val_m_better - val_m)) / abs(val_m_better) < 1e-8
        @test abs(part(val_s_better - val_s)) / abs(val_s_better) < 1e-8
        @test abs(part(val_t_better - val_t)) / abs(val_t_better) < 1e-6
    end
end

@testitem "gamma integrated 1D" setup = [Parameters1D] begin
    val_d = WeakCouplingParquet.gamma2_d(u, mu, beta, v, vp, w, k, kp, q)
    val_m = WeakCouplingParquet.gamma2_m(u, mu, beta, v, vp, w, k, kp, q)
    val_s = WeakCouplingParquet.gamma2_s(u, mu, beta, v, vp, w, k, kp, q)
    val_t = WeakCouplingParquet.gamma2_t(u, mu, beta, v, vp, w, k, kp, q)

    val_d_better = WeakCouplingParquet.gamma2_d(u, mu, beta, v, vp, w, k, kp, q; n=41)
    val_m_better = WeakCouplingParquet.gamma2_m(u, mu, beta, v, vp, w, k, kp, q; n=41)
    val_s_better = WeakCouplingParquet.gamma2_s(u, mu, beta, v, vp, w, k, kp, q; n=41)
    val_t_better = WeakCouplingParquet.gamma2_t(u, mu, beta, v, vp, w, k, kp, q; n=41)

    for part in (real, imag)
        @test abs(part(val_d_better - val_d)) / abs(val_d_better) < 1e-8
        @test abs(part(val_m_better - val_m)) / abs(val_m_better) < 1e-8
        @test abs(part(val_s_better - val_s)) / abs(val_s_better) < 1e-8
        @test abs(part(val_t_better - val_t)) / abs(val_t_better) < 1e-6
    end
end


# this is currently pretty slow. not terribly useful

# @testitem "gamma integrated 2D" setup = [Parameters2D] begin
#     abstol = 1e-2
#     reltol = 1e-2

#     val_d, err_d = WeakCouplingParquet.gamma2_d(u, mu, beta, v, vp, w, k, kp, q)
#     val_m, err_m = WeakCouplingParquet.gamma2_m(u, mu, beta, v, vp, w, k, kp, q)
#     val_s, err_s = WeakCouplingParquet.gamma2_s(u, mu, beta, v, vp, w, k, kp, q)
#     val_t, err_t = WeakCouplingParquet.gamma2_t(u, mu, beta, v, vp, w, k, kp, q)

#     for part in (real, imag)
#         @test abs(part(err_d)) < max(abstol, reltol * abs(part(val_d)))
#         @test abs(part(err_m)) < max(abstol, reltol * abs(part(val_m)))
#         @test abs(part(err_s)) < max(abstol, reltol * abs(part(val_s)))
#         @test abs(part(err_t)) < max(abstol, reltol * abs(part(val_t)))
#     end
# end

# @testitem "no nan value in gamma when v == vp and k != kp" begin
#     using StaticArrays

#     u = 0.3
#     beta = 0.5
#     val_d, err_d = WeakCouplingParquet.gamma2_d(u, u / 2, beta, 17 * im * pi / beta, 17 * im * pi / beta, -12 * im * pi / beta, SA[pi/2, pi/2], SA[0.0, 0.0], SA[0.0, 0.0])

#     @test !isnan(val_d)
#     @test !isnan(err_d)
# end

# @testitem "no nan value in phi when v == vp" begin
#     using StaticArrays

#     u = 0.3
#     beta = 0.5
#     val_d = WeakCouplingParquet.phi2_m(u, u / 2, beta, -13 * im * pi / beta, -13 * im * pi / beta, 30 * im * pi / beta, SA[0π/4, 2π/4], SA[2π/4, 0π/4], SA[4π/4, 6π/4], SA[1.9486832980505138, 0.05131670194948623], SA[1.0, 1.0], SA[1.0, 1.0])

#     @test !isnan(val_d)
# end

# @testitem "no nan value in integrated phi when v == vp" begin
#     using StaticArrays

#     u = 0.3
#     beta = 0.5
#     abstol = 1e-5
#     reltol = 1e-5
#     val_d, err_d = WeakCouplingParquet.phi2_m(u, u / 2, beta, -13 * im * pi / beta, -13 * im * pi / beta, 30 * im * pi / beta, SA[0π/4, 2π/4], SA[2π/4, 0π/4], SA[4π/4, 6π/4])

#     @test !isnan(val_d)
# end
