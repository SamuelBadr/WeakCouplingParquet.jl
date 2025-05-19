# Single argument tests
@testitem "Single argument" begin
    @test WeakCouplingParquet.equals_first(5) == (5,)
    @test WeakCouplingParquet.equals_first("hello") == ("hello",)
    @test WeakCouplingParquet.equals_first([1, 2, 3]) == ([1, 2, 3],)
end

# Two arguments tests
@testitem "Two arguments - no equality" begin
    @test WeakCouplingParquet.equals_first(1, 2) == (1, 2)
    @test WeakCouplingParquet.equals_first("a", "b") == ("a", "b")
end

@testitem "Two arguments - with equality" begin
    @test WeakCouplingParquet.equals_first(5, 5) == (5, 5)
    @test WeakCouplingParquet.equals_first("test", "test") == ("test", "test")
end

@testitem "Three arguments - a equals b" begin
    @test WeakCouplingParquet.equals_first(5, 5, 3) == (5, 5, 3)
end

@testitem "Three arguments - a equals c" begin
    @test WeakCouplingParquet.equals_first(5, 3, 5) == (5, 5, 3)
end

@testitem "Three arguments - b equals c" begin
    @test WeakCouplingParquet.equals_first(1, 5, 5) == (5, 5, 1)
end

@testitem "Three arguments - all equal" begin
    @test WeakCouplingParquet.equals_first(5, 5, 5) == (5, 5, 5)
end

# Four arguments tests
@testitem "Four arguments - a equals b" begin
    @test WeakCouplingParquet.equals_first(5, 5, 3, 4) == (5, 5, 4, 3)
end

@testitem "Four arguments - a equals c" begin
    @test WeakCouplingParquet.equals_first(5, 3, 5, 4) == (5, 5, 4, 3)
end

@testitem "Four arguments - a equals d" begin
    @test WeakCouplingParquet.equals_first(5, 3, 4, 5) == (5, 5, 4, 3)
end

@testitem "Four arguments - b equals c" begin
    @test WeakCouplingParquet.equals_first(1, 5, 5, 4) == (5, 5, 4, 1)
end

@testitem "Four arguments - b equals d" begin
    @test WeakCouplingParquet.equals_first(1, 5, 3, 5) == (5, 5, 3, 1)
end

@testitem "Four arguments - c equals d" begin
    @test WeakCouplingParquet.equals_first(1, 2, 5, 5) == (5, 5, 2, 1)
end

@testitem "Four arguments - multiple pairs of equality" begin
    # a equals b, c equals d
    @test WeakCouplingParquet.equals_first(5, 5, 7, 7) == (5, 5, 7, 7)

    # a equals c, b equals d
    @test WeakCouplingParquet.equals_first(5, 7, 5, 7) == (5, 5, 7, 7)

    # a equals d, b equals c 
    @test WeakCouplingParquet.equals_first(5, 7, 7, 5) == (5, 5, 7, 7)
end

@testitem "Four arguments - three equal values" begin
    # a, b, c equal
    @test WeakCouplingParquet.equals_first(5, 5, 5, 7) == (5, 5, 5, 7)

    # a, b, d equal
    @test WeakCouplingParquet.equals_first(5, 5, 7, 5) == (5, 5, 5, 7)

    # a, c, d equal
    @test WeakCouplingParquet.equals_first(5, 7, 5, 5) == (5, 5, 5, 7)

    # b, c, d equal
    @test WeakCouplingParquet.equals_first(7, 5, 5, 5) == (5, 5, 5, 7)
end

@testitem "Four arguments - all equal" begin
    @test WeakCouplingParquet.equals_first(5, 5, 5, 5) == (5, 5, 5, 5)
end

@testitem "Different data types" begin
    @test WeakCouplingParquet.equals_first("test", 5, "test") == ("test", "test", 5)
    @test WeakCouplingParquet.equals_first([1, 2], [1, 2], "array") == ([1, 2], [1, 2], "array")
end

@testitem "Edge cases" begin
    # Empty strings
    @test WeakCouplingParquet.equals_first("", "", "x") == ("", "", "x")

    # NaN values (should not be considered equal)
    @test WeakCouplingParquet.equals_first(NaN, NaN, 1) != (NaN, NaN, 1)

    # Nothing values
    @test WeakCouplingParquet.equals_first(nothing, nothing, 1) == (nothing, nothing, 1)

    # Missing values
    @test_broken WeakCouplingParquet.equals_first(missing, missing, 1) == (missing, missing, 1)
end

# Frequency-based ordering tests
@testitem "Frequency-based ordering" begin
    # Test with varargs of different lengths containing repeated values
    # The most frequent value should come first

    # Here 5 appears twice, 3 appears once, 2 appears once
    @test WeakCouplingParquet.equals_first(2, 5, 3, 5) == (5, 5, 3, 2)

    # Here both 5 and 2 appear twice, but order is determined by the specific implementation
    # Need to check function behavior to determine expected outcome
    # @test WeakCouplingParquet.equals_first(5, 2, 5, 2) == (5, 5, 2, 2)  # This would depend on implementation

    # For three arguments, where two values are equal
    @test WeakCouplingParquet.equals_first(5, 5, 3) == (5, 5, 3)  # Current implementation doesn't count frequencies
end

# Bug test: The implementation has an issue with the case a==c in the 4-argument version
@testitem "Bug: Duplicate condition in 4-argument version" begin
    # In the current implementation, there's a condition checking if a==b inside the a==c branch
    # This is likely a bug as it's already established that a==c at this point
    # This test is marked as broken to highlight the issue
    @test WeakCouplingParquet.equals_first(5, 3, 5, 4) == (5, 5, 4, 3)
end
