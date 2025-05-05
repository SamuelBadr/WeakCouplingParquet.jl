using MathLink
using JuliaFormatter
using Scratch

function make_julia_function(wolfram_path, out_path, function_name)
    wolfram_expr_string = read(wolfram_path, String)
    wolfram_expr = MathLink.parseexpr(wolfram_expr_string)
    @show wolfram_path function_name
    expr = W2JuliaExpr(wolfram_expr)
    expr = replace(string(expr), "+ -1 *" => "-", "-1 * " => "-", "+ -1" => "- ")
    if occursin("k2", expr) && occursin("k3", expr)
        signature = "(u, mu, beta, v, vp, w, k::SVector, kp::SVector, q::SVector, k1::SVector, k2::SVector, k3::SVector)"
        normalization_constant = "cnorm = 1 / (2pi)^length(k)"
    elseif occursin("k1", expr)
        signature = "(u, mu, beta, v, vp, w, k::SVector, kp::SVector, q::SVector, k1::SVector)"
        normalization_constant = "cnorm = 1 / (2pi)^length(k)"
    else
        signature = "(u, mu, beta, v, vp, w)"
        normalization_constant = ""
    end
    function_string = """
    @cse function $(function_name)$(signature)
        $normalization_constant
        $expr
    end
    """
    function_string = format_text(function_string)
    write(out_path, function_string)
    return out_path
end

function make_julia_functions()
    wolfram_expr_dir = joinpath(@__DIR__, "wolfram_expressions")
    julia_function_dir = get_scratch!(Base.UUID("62b8c3c3-4004-406c-9797-34cb1d81b562"), "julia_functions")
    includes = String[]
    for in_file in readdir(wolfram_expr_dir)
        fun_name = splitext(in_file)[1]
        out_file = "$(fun_name).jl"
        wolfram_path = joinpath(wolfram_expr_dir, in_file)
        julia_path = joinpath(julia_function_dir, out_file)
        make_julia_function(wolfram_path, julia_path, fun_name)
        push!(includes, "include(\"$out_file\")")
    end
    write(joinpath(julia_function_dir, "_includes.jl"), join(includes, "\n"))
end

make_julia_functions()
