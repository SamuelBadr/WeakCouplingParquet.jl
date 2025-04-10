module WeakCouplingParquet

using MathLink
using JuliaFormatter
using StaticArrays
using HCubature
using Cubature
using Integrals
using FastGaussQuadrature
using MCIntegration
using Cuba

disp(k::AbstractVector) = -2 * sum(cospi, k)

fermidist(x::Number, beta::Number) = 1 / (exp(beta * x) + 1)
bosedist(x::Number, beta::Number) = 1 / (exp(beta * x) - 1)

function make_julia_function(wolfram_path, out_path, function_name)
    wolfram_expr_string = read(wolfram_path, String)
    wolfram_expr = MathLink.parseexpr(wolfram_expr_string)
    expr = W2JuliaExpr(wolfram_expr)
    expr = replace(string(expr), "+ -1 *" => "-", "-1 * " => "-", "+ -1" => "- ")
    signature = if occursin("k2", expr) && occursin("k3", expr)
        "(u, mu, beta, v, vp, w, k::SVector, kp::SVector, q::SVector, k1::SVector, k2::SVector, k3::SVector)"
    else
        "(u, mu, beta, v, vp, w, k::SVector, kp::SVector, q::SVector, k1::SVector)"
    end
    function_string = """
    function $(function_name)$(signature)
        $expr
    end
    """
    function_string = format_text(function_string)
    write(out_path, function_string)
    return out_path
end

wolfram_expr_dir = joinpath(@__DIR__, "wolfram_expressions")
julia_expr_dir = joinpath(@__DIR__, "julia_expressions")
for in_file in readdir(wolfram_expr_dir)
    fun_name = splitext(in_file)[1]
    out_file = "$(fun_name).jl"
    wolfram_path = joinpath(wolfram_expr_dir, in_file)
    julia_path = joinpath(julia_expr_dir, out_file)
    make_julia_function(wolfram_path, julia_path, fun_name)
    include(julia_path)
end

for fun in (:phi2_d, :phi2_m, :phi2_s, :phi2_t)
    @eval function $(fun)(dx, x, p)
        (; u, mu, beta, v, vp, w, k, kp, q) = p
        dim = length(k)
        inds = SVector(ntuple(identity, dim)...)
        inds1 = inds .+ 0dim
        inds2 = inds .+ 1dim
        inds3 = inds .+ 2dim
        # @show size(x, 2)
        @inbounds Threads.@threads for i in axes(x, 2)
            k1 = x[inds1, i]
            k2 = x[inds2, i]
            k3 = x[inds3, i]

            while true
                res = $(fun)(u, mu, beta, v, vp, w, k, kp, q, k1, k2, k3)
                if isfinite(res)
                    dx[1, i], dx[2, i] = real(res), imag(res)
                    break
                end
                # if the integrand blows up, slightly
                # perturb kp until it doesn't
                kp = kp + 0.0000001 * @SVector randn(dim)
            end
        end
    end

    @eval function $(fun)(u, mu, beta, v, vp, w, k::SVector{dim}, kp::SVector{dim}, q::SVector{dim}, alg=CubatureJLh(); kwargs...) where dim
        prototype = @SMatrix zeros(2, 0)
        intfun = BatchIntegralFunction($fun, prototype)

        corner1 = @SVector fill(0.0, 3dim)
        corner2 = @SVector fill(2.0, 3dim)
        domain = (corner1, corner2)
        p = (; u, mu, beta, v, vp, w, k, kp, q)
        prob = IntegralProblem(intfun, domain, p)

        sol = solve(prob, alg; kwargs...)
        result = sol.u[1] + im * sol.u[2]
        residuum = sol.resid[1] + im * sol.resid[2]
        return result, residuum
    end
end

for fun in (:full2_d, :full2_m, :full2_s, :full2_t, :gamma2_d, :gamma2_m, :gamma2_s, :gamma2_t)
    @eval function $(fun)(dx, x, p)
        (; u, mu, beta, v, vp, w, k, kp, q) = p
        dim = length(k)
        inds = SVector(ntuple(identity, dim)...)
        # @show size(x, 2)
        @inbounds Threads.@threads for i in axes(x, 2)
            k1 = x[inds, i]
            while true
                res = $(fun)(u, mu, beta, v, vp, w, k, kp, q, k1)
                if isfinite(res)
                    dx[1, i], dx[2, i] = real(res), imag(res)
                    break
                end
                # if the integrand blows up, slightly
                # perturb kp until it doesn't
                kp = kp + 0.0000001 * @SVector randn(dim)
            end
        end
    end

    @eval function $(fun)(u, mu, beta, v, vp, w, k::SVector{dim}, kp::SVector{dim}, q::SVector{dim}, alg=CubatureJLh(); kwargs...) where dim
        prototype = @SMatrix zeros(2, 0)
        intfun = BatchIntegralFunction($fun, prototype)

        corner1 = @SVector fill(0.0, dim)
        corner2 = @SVector fill(2.0, dim)
        domain = (corner1, corner2)
        p = (; u, mu, beta, v, vp, w, k, kp, q)
        prob = IntegralProblem(intfun, domain, p)

        sol = solve(prob, alg; kwargs...)
        result = sol.u[1] + im * sol.u[2]
        residuum = sol.resid[1] + im * sol.resid[2]
        return result, residuum
    end
end

end # module WeakCouplingParquet
