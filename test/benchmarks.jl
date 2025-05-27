@testitem "phi integrated 1D benchmark" setup = [Parameters2D] begin
    abstol = 1e-3
    reltol = 1e-3

    @time WeakCouplingParquet.phi2_d(u, mu, beta, v, vp, w, k, kp, q)
    @time WeakCouplingParquet.phi2_m(u, mu, beta, v, vp, w, k, kp, q)
    @time WeakCouplingParquet.phi2_s(u, mu, beta, v, vp, w, k, kp, q)
    @time WeakCouplingParquet.phi2_t(u, mu, beta, v, vp, w, k, kp, q)
end
