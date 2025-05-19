@inline nf(x::Number, beta::Number) = 1 / (exp(beta * x) + 1)

@inline function dnf(x::Number, beta::Number)
    f = nf(x, beta)
    beta * f * (f - 1)
end

@inline function ddnf(x::Number, beta::Number)
    f = nf(x, beta)
    beta^2 * f * (f - 1) * (2f - 1)
end

@inline function dddnf(x::Number, beta::Number)
    f = nf(x, beta)
    beta^3 * f * (f - 1) * (6f * (f - 1) + 1)
end


"""
    matsubara_sum_fermi(beta, a)

Compute the Matsubara sum
```math
    1/beta * sum(1 / (z - a), z)
```
at temperature 1/`beta` assuming `z` is fermionic.
"""
function matsubara_sum_fermi(beta::Real, a::T) where {T<:Number}
    nf(a, beta)
end

"""
    matsubara_sum_fermi(beta, a, b)

Compute the Matsubara sum
```math
    1/beta * sum(1 / ((z - a) (z - b)), z)
```
at temperature 1/`beta` assuming `z` is fermionic.
"""
function matsubara_sum_fermi(beta::Real, a::T, b::T) where {T<:Number}
    if a == b
        dnf(a, beta)
    else
        nf(a, beta) / (a - b) +
        nf(b, beta) / (b - a)
    end
end

matsubara_sum_fermi(beta::Real, a, b) = matsubara_sum_fermi(beta, promote(a, b)...)

"""
    matsubara_sum_fermi(beta, a, b, c)

Compute the Matsubara sum
```math
    1/beta * sum(1 / ((z - a) (z - b) (z - c)), z)
```
at temperature 1/`beta` assuming `z` is fermionic.
"""
function matsubara_sum_fermi(beta::Real, a::T, b::T, c::T) where {T<:Number}
    a, b, c = equals_first(a, b, c)
    if a == b == c
        1 / 2 * ddnf(a, beta)
    elseif a == b
        -(nf(a, beta) - nf(c, beta)) / (a - c)^2 +
        dnf(a, beta) / (a - c)
    else
        nf(a, beta) / ((a - b) * (a - c)) +
        nf(b, beta) / ((b - a) * (b - c)) +
        nf(c, beta) / ((c - a) * (c - b))
    end
end

matsubara_sum_fermi(beta::Real, a, b, c) = matsubara_sum_fermi(beta, promote(a, b, c)...)

"""
    matsubara_sum_fermi(beta, a, b, c, d)

Compute the Matsubara sum
```math
    1/beta * sum(1 / ((z - a) (z - b) (z - c) (z - d)), z)
```
at temperature 1/`beta` assuming `z` is fermionic.
"""
function matsubara_sum_fermi(beta::Real, a::T, b::T, c::T, d::T) where {T<:Number}
    a, b, c, d = equals_first(a, b, c, d)
    if a == b == c == d
        1 / 6 * dddnf(a, beta)
    elseif a == b == c
        (nf(a, beta) - nf(d, beta)) / (a - d)^3 -
        dnf(a, beta) / (a - d)^2 +
        1 / 2 * ddnf(a, beta) / (a - d)
    elseif a == b
        if c == d
            -2 * (nf(a, beta) - nf(c, beta)) / (a - c)^3 +
            (dnf(a, beta) + dnf(c, beta)) / (a - c)^2
        else
            (-2a + c + d) * nf(a, beta) / ((a - c)^2 * (a - d)^2) +
            nf(c, beta) / ((a - c)^2 * (c - d)) +
            nf(d, beta) / ((a - d)^2 * (d - c)) +
            dnf(a, beta) / ((a - c) * (a - d))
        end
    else
        nf(a, beta) / ((a - b) * (a - c) * (a - d)) +
        nf(b, beta) / ((b - a) * (b - c) * (b - d)) +
        nf(c, beta) / ((c - a) * (c - b) * (c - d)) +
        nf(d, beta) / ((d - a) * (d - b) * (d - c))
    end
end

matsubara_sum_fermi(beta::Real, a, b, c, d) = matsubara_sum_fermi(beta, promote(a, b, c, d)...)

"""
    equals_first(a_1, a_2, …, a_n)

Reorder its arguments so that all equal values appear first, with the most-frequent values
coming before less-frequent ones. Return a tuple of the reordered values.
"""
function equals_first(a)
    (a,)
end

function equals_first(a, b)
    (a, b)
end

function equals_first(a, b, c)
    if a == b
        (a, b, c)
    elseif a == c
        (a, c, b)
    else
        (b, c, a)
    end
end

function equals_first(a, b, c, d)
    if a == b
        if a == c
            (a, b, c, d)
        else
            (a, b, d, c)
        end
    elseif a == c
        (a, c, d, b)
    elseif a == d
        (a, d, c, b)
    elseif b == c
        (b, c, d, a)
    elseif b == d
        (b, d, c, a)
    else
        (c, d, b, a)
    end
end
