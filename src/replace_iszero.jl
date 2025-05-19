using MacroTools

"""
@replace_iszero func_def

Macro that walks through a function definition `func_def`, inserts a local absolute
 tolerance `atol = eps()` at the top of the function body, and replaces every
 `iszero(x)` call by `isapprox(x, 0; atol=atol)`.

# Arguments

- `func_def`: a Julia `function` definition expression.

# Example

```julia
@replace_iszero function foo(x)
    a = iszero(x)
    b = iszero(x + 2)
    return a || b
end
```

After macro expansion, calls become:

```julia
function foo(x)
    atol = eps()
    a = isapprox(x, 0; atol=atol)
    b = isapprox(x + 2, 0; atol=atol)
    return a || b
end
```
"""
macro replace_iszero(func_def)
    # Ensure we have a function definition
    if !(func_def isa Expr && func_def.head == :function)
        error("@replace_iszero must be applied to a function definition")
    end

    # Destructure the function expression: signature and body block
    sig = func_def.args[1]
    body = func_def.args[2:end]

    # Create a new body with `atol = eps()` inserted
    new_body = Expr(:block, :(atol = eps()), body...)

    # Walk and transform the new body, replacing iszero(x) with isapprox(x, 0; atol=atol)
    transformed_body = MacroTools.postwalk(new_body) do node
        if node isa Expr && node.head == :call && node.args[1] == :iszero
            arg = node.args[2]
            return :(isapprox($arg, 0; atol))
        else
            return node
        end
    end

    # Reconstruct the function definition
    new_func = Expr(:function, sig, transformed_body)

    # Return with hygiene
    return esc(new_func)
end
