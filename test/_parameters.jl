@testsnippet Parameters2D begin
    using StaticArrays

    beta = 10.0
    u = 3.0
    mu = u / 2

    n = 1
    np = 11
    m = -4

    v = im * n * pi / beta
    vp = im * np * pi / beta
    w = im * m * pi / beta

    k = SA[1.2, 1.9]
    kp = SA[0.5, 1.4]
    q = SA[0.3, 0.8]

    k1 = SA[0.1, 0.9]
    k2 = SA[0.2, 1.0]
    k3 = SA[0.3, 1.1]
end

@testsnippet Parameters1D begin
    using StaticArrays

    beta = 1.0
    u = 10.0
    mu = u / 2

    n = 1
    np = 11
    m = -4

    v = im * n * pi / beta
    vp = im * np * pi / beta
    w = im * m * pi / beta

    k = SA[1.2]
    kp = SA[0.5]
    q = SA[0.3]
end

@testsnippet Parameters0D begin
    using StaticArrays

    beta = 1.0
    u = 10.0
    mu = u / 2

    n = 1
    np = 11
    m = -4

    v = im * n * pi / beta
    vp = im * np * pi / beta
    w = im * m * pi / beta

    k = SA{Float64}[]
    kp = SA{Float64}[]
    q = SA{Float64}[]
end
