function full2_t(u, mu, beta, v, vp, w, k::SVector, kp::SVector, q::SVector, k1::SVector)
    cnorm = 1 / (2pi)^length(k)
    cnorm * (
        ifelse(
            and(iszero(-1v + vp), iszero(-1k + kp)),
            -(u ^ 2 * dfermidist(-1mu + disp(k1), beta)),
            (
                u ^ 2 * (
                    fermidist(-1mu + disp(k1), beta) -
                    fermidist(-1mu + v - vp + disp(-1k + k1 + kp), beta)
                )
            ) * (v - vp - disp(k1) + disp(-1k + k1 + kp)) ^ -1,
        ) + ifelse(
            and(iszero(-1v - vp + w), iszero(-1k - kp + q)),
            u ^ 2 * dfermidist(-1mu + disp(k1), beta),
            (
                u ^ 2 * (
                    -fermidist(-1mu + disp(k1), beta) +
                    fermidist(-1mu + v + vp - w + disp(-1k + k1 - kp + q), beta)
                )
            ) * (v + vp - w - disp(k1) + disp(-1k + k1 - kp + q)) ^ -1,
        )
    )
end
