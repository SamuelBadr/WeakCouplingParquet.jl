module WeakCouplingParquet

using StaticArrays
# using Cubature
# using Integrals
using CommonSubexpressions
using Scratch

@inline disp(k::AbstractVector) = -2 * sum(cos, k)
@inline disp(k::SVector{0}) = 0.0

@inline fermidist(x::Number, beta::Number) = fermidist(complex(x), beta)
@inline fermidist(x::Complex, beta::Number) = @fastmath 1 / (exp(beta * x) + 1)

@inline function dfermidist(x::Number, beta::Number)
    f = fermidist(x, beta)
    beta * f * (f - 1)
end

@inline and(x, b) = x && b

# const julia_function_dir = @get_scratch!("julia_functions")
# const julia_function_dir = joinpath(@__DIR__, "functions")
# include(joinpath(julia_function_dir, "_includes.jl"))
include("functions/_includes.jl")


include("integration.jl")

function greensfunction(mu, v, k)
    1 / (v - disp(k) + mu)
end

function chi0_ph(beta, mu, v, vp, w, k, kp, q)
    if v ≈ vp && k ≈ kp
        beta * chi0_ph(beta, mu, v, w, k, q)
    else
        complex(0.0)
    end
end

function chi0_ph(beta, mu, v, w, k, q)
    greensfunction(mu, v, k) * greensfunction(mu, v + w, k + q)
end

const chi0_d = chi0_ph
const chi0_m = chi0_ph

function chi0_s(beta, mu, v, vp, w, k, kp, q)
    if v ≈ vp && k ≈ kp
        beta * chi0_s(beta, mu, v, w, k, q)
    else
        complex(0.0)
    end
end

function chi0_s(beta, mu, v, w, k, q)
    -1 / 2 * greensfunction(mu, v, k) * greensfunction(mu, v + w, k + q)
end

function chi0_t(beta, mu, v, vp, w, k, kp, q)
    if v ≈ vp && k ≈ kp
        beta * chi0_t(beta, mu, v, w, k, q)
    else
        complex(0.0)
    end
end

function chi0_t(beta, mu, v, w, k, q)
    1 / 2 * greensfunction(mu, v, k) * greensfunction(mu, v + w, k + q)
end

for fun in (:phi2_d, :phi2_m, :phi2_s, :phi2_t)
    # @eval function $(fun)(dx, x, p)
    #     (; u, mu, beta, v, vp, w, k, kp, q) = p
    #     dim = length(k)
    #     inds = SVector(ntuple(identity, dim)...)
    #     inds1 = inds .+ 0dim
    #     inds2 = inds .+ 1dim
    #     inds3 = inds .+ 2dim
    #     @inbounds Threads.@threads for i in axes(x, 2)
    #         k1 = x[inds1, i]
    #         k2 = x[inds2, i]
    #         k3 = x[inds3, i]
    #         res = $(fun)(u, mu, beta, v, vp, w, k, kp, q, k1, k2, k3)
    #         if isnan(res)
    #             error("$((u, mu, beta, v, vp, w, k, kp, q, k1, k2, k3))")
    #         end
    #         dx[1, i], dx[2, i] = real(res), imag(res)
    #     end
    # end

    # @eval function $(fun)(u, mu, beta, v, vp, w, k::SVector{dim}, kp::SVector{dim}, q::SVector{dim}, alg=CubatureJLh(); kwargs...) where {dim}
    #     prototype = @SMatrix zeros(2, 0)
    #     intfun = BatchIntegralFunction($fun, prototype)

    #     corner1 = @SVector fill(0.0, 3dim)
    #     corner2 = @SVector fill(2.0, 3dim)
    #     domain = (corner1, corner2)
    #     p = (; u, mu, beta, v, vp, w, k, kp, q)
    #     prob = IntegralProblem(intfun, domain, p)

    #     sol = solve(prob, alg; kwargs...)
    #     result = sol.u[1] + im * sol.u[2]
    #     residuum = sol.resid[1] + im * sol.resid[2]
    #     return result, residuum
    # end

    @eval function $(fun)(u, mu, beta, v, vp, w, k::SVector{dim}, kp::SVector{dim}, q::SVector{dim}; kwargs...) where {dim}
        inds1 = SVector{dim,Int}(1:dim...) .+ 0dim
        inds2 = SVector{dim,Int}(1:dim...) .+ 1dim
        inds3 = SVector{dim,Int}(1:dim...) .+ 2dim

        func(k1k2k3) = $(fun)(u, mu, beta, v, vp, w, k, kp, q, k1k2k3[inds1], k1k2k3[inds2], k1k2k3[inds3])

        return integrate(Val(3dim), ComplexF64, func; kwargs...)
    end
end

for fun in (:full2_d, :full2_m, :full2_s, :full2_t, :gamma2_d, :gamma2_m, :gamma2_s, :gamma2_t)
    # @eval function $(fun)(dx, x, p)
    #     (; u, mu, beta, v, vp, w, k, kp, q) = p
    #     dim = length(k)
    #     inds = SVector(ntuple(identity, dim)...)
    #     @inbounds Threads.@threads for i in axes(x, 2)
    #         k1 = x[inds, i]
    #         res = $(fun)(u, mu, beta, v, vp, w, k, kp, q, k1)
    #         dx[1, i], dx[2, i] = real(res), imag(res)
    #     end
    # end

    # @eval function $(fun)(u, mu, beta, v, vp, w, k::SVector{dim}, kp::SVector{dim}, q::SVector{dim}, alg=CubatureJLh(); kwargs...) where {dim}
    #     prototype = @SMatrix zeros(2, 0)
    #     intfun = BatchIntegralFunction($fun, prototype)

    #     corner1 = @SVector fill(0.0, dim)
    #     corner2 = @SVector fill(2pi, dim)
    #     domain = (corner1, corner2)
    #     p = (; u, mu, beta, v, vp, w, k, kp, q)
    #     prob = IntegralProblem(intfun, domain, p)

    #     sol = solve(prob, alg; kwargs...)
    #     result = sol.u[1] + im * sol.u[2]
    #     residuum = sol.resid[1] + im * sol.resid[2]
    #     return result, residuum
    # end

    @eval function $(fun)(u, mu, beta, v, vp, w, k::SVector{dim}, kp::SVector{dim}, q::SVector{dim}; kwargs...) where {dim}
        func(k1) = $(fun)(u, mu, beta, v, vp, w, k, kp, q, k1)

        return integrate(Val(dim), ComplexF64, func; kwargs...)
    end
end

end # module WeakCouplingParquet
