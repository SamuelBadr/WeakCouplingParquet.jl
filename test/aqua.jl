@testitem "Aqua.jl" begin
    using Aqua
    Aqua.test_all(WeakCouplingParquet)
end
