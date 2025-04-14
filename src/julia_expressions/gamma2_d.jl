function gamma2_d(u, mu, beta, v, vp, w, k::SVector, kp::SVector, q::SVector, k1::SVector)
    cnorm = 1 / (2pi)^length(k)
    u +
    (
        cnorm *
        u ^ 2 *
        (
            -fermidist(-1mu + disp(k1), beta) +
            fermidist(mu + v + vp + w - disp(k - k1 + kp + q), beta)
        )
    ) * (2mu + v + vp + w - disp(k1) - disp(k - k1 + kp + q)) ^ -1 - (
        2 *
        cnorm *
        ifelse(
            and(iszero(-1v + vp), iszero(-1k + kp)),
            u ^ 2 * dfermidist(-1mu + disp(k1), beta),
            (
                u ^ 2 * (
                    -fermidist(-1mu + disp(k1), beta) +
                    fermidist(-1mu + v - vp + disp(-1k + k1 + kp), beta)
                )
            ) * (v - vp - disp(k1) + disp(-1k + k1 + kp)) ^ -1,
        )
    )
end
