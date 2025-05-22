using QuadGK

function integrate(::Val{D}, ::Type{T}, f; n=21, rule=gauss(n, 0.0, 2pi)) where {D,T}
    xs, ws = rule

    total = zero(T)
    for (x, w) in zip(Iterators.product(ntuple(_ -> xs, D)...), Iterators.product(ntuple(_ -> ws, D)...))
        xvec = SVector(x)
        wprod = prod(w)
        total += wprod * f(xvec)
    end
    return total
end
