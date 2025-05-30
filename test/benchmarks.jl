# @testitem "phi integrated 1D benchmark" setup = [Parameters1D] begin
#     time_d = @elapsed WeakCouplingParquet.phi2_d(u, mu, beta, v, vp, w, k, kp, q)
#     time_m = @elapsed WeakCouplingParquet.phi2_m(u, mu, beta, v, vp, w, k, kp, q)
#     time_s = @elapsed WeakCouplingParquet.phi2_s(u, mu, beta, v, vp, w, k, kp, q)
#     time_t = @elapsed WeakCouplingParquet.phi2_t(u, mu, beta, v, vp, w, k, kp, q)

#     @test time_d < 0.1
#     @test time_m < 0.1
#     @test time_s < 0.1
#     @test time_t < 0.1
# end
