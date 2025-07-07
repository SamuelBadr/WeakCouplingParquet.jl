# @testitem "Bethe-Salpeter Equation D = 0" begin
using StaticArrays
using WeakCouplingParquet
using CairoMakie

beta = 20.0
u = 0.5
mu = 0.2

nfermi_box = 32
nfermi_sum = 640

vs = (-(nfermi_box - 1):2:+(nfermi_box - 1)) .* im * pi / beta
vs_sum = (-(nfermi_sum - 1):2:+(nfermi_sum - 1)) .* im * pi / beta
ws = [10] .* im * pi / beta
k = SA{Float64}[]

bse_res = map(Iterators.product(vs, vs, ws)) do (v, vp, w)
    res = 0.0 + 0.0im
    for v1 in vs_sum
        res += WeakCouplingParquet.gamma2_m(u, mu, beta, v, v1, w, k, k, k) *
               WeakCouplingParquet.chi0_m(beta, mu, v1, w, k, k) *
               WeakCouplingParquet.full2_m(u, mu, beta, v1, vp, w, k, k, k)
    end
    res / beta
end

phi = map(Iterators.product(vs, vs, ws)) do (v, vp, w)
    res = WeakCouplingParquet.phi2_m(u, mu, beta, v, vp, w, k, k, k)
end

##

fig = Figure(size=(400, 800))
res = real(bse_res)[:, :, 1]
ana = real(phi)[:, :, 1]
diff = res - ana
ext = extrema([res; ana])
colorrange = ext

ax1 = Axis(fig[1, 1]; aspect=1, title="BSE Result")
hm1 = heatmap!(ax1, imag.(vs), imag.(vs), res; colorrange)
Colorbar(fig[1, 1][1, 2], hm1)

ax2 = Axis(fig[2, 1]; aspect=1, title="Analytic Φ")
hm2 = heatmap!(ax2, imag.(vs), imag.(vs), ana; colorrange)
Colorbar(fig[2, 1][1, 2], hm2)

ax3 = Axis(fig[3, 1]; aspect=1, title="Difference")
hm3 = heatmap!(ax3, imag.(vs), imag.(vs), abs.(diff))
Colorbar(fig[3, 1][1, 2], hm3)

fig
# end
