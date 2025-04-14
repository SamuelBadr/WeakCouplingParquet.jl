function phi2_t(
    u,
    mu,
    beta,
    v,
    vp,
    w,
    k::SVector,
    kp::SVector,
    q::SVector,
    k1::SVector,
    k2::SVector,
    k3::SVector,
)
    cnorm = 1 / (2pi)^length(k)
    cnorm * (
        (
            cnorm ^ 2 *
            fermidist(-1mu + disp(k1), beta) *
            (
                ifelse(
                    and(iszero(mu + vp - disp(k1)), iszero(-1k1 + kp)),
                    -(u ^ 2 * dfermidist(-1mu + disp(k3), beta)),
                    (
                        u ^ 2 * (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                        )
                    ) * (-1mu - vp + disp(k1) - disp(k3) + disp(-1k1 + k3 + kp)) ^ -1,
                ) + ifelse(
                    and(iszero(mu - vp + w - disp(k1)), iszero(-1k1 - kp + q)),
                    u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                    (
                        u ^ 2 * (
                            -fermidist(-1mu + disp(k3), beta) +
                            fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                        )
                    ) *
                    (-1mu + vp - w + disp(k1) - disp(k3) + disp(-1k1 + k3 - kp + q)) ^ -1,
                )
            ) *
            (
                ifelse(
                    and(iszero(mu - v + w - disp(k1)), iszero(-1k - k1 + q)),
                    u ^ 2 * dfermidist(-1mu + disp(k2), beta),
                    (
                        u ^ 2 * (
                            -fermidist(-1mu + disp(k2), beta) +
                            fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
                        )
                    ) * (-1mu + v - w + disp(k1) - disp(k2) + disp(-1k - k1 + k2 + q)) ^ -1,
                ) + ifelse(
                    and(iszero(-1mu - v + disp(k1)), iszero(-1k + k1)),
                    -(u ^ 2 * dfermidist(-1mu + disp(k2), beta)),
                    (
                        u ^ 2 * (
                            fermidist(-1mu + disp(k2), beta) -
                            fermidist(-1mu + disp(-1k + k1 + k2), beta)
                        )
                    ) * (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) ^ -1,
                )
            )
        ) * (2 * (2mu + w - disp(k1) - disp(-1k1 + q))) ^ -1 - (
            (
                cnorm ^ 2 *
                fermidist(mu - disp(-1k1 + q), beta) *
                (
                    ifelse(
                        and(iszero(mu - v + w - disp(-1k1 + q)), iszero(-1k + k1)),
                        -(u ^ 2 * dfermidist(-1mu + disp(k2), beta)),
                        (
                            u ^ 2 * (
                                fermidist(-1mu + disp(k2), beta) -
                                fermidist(-1mu + disp(-1k + k1 + k2), beta)
                            )
                        ) *
                        (-1mu + v - w - disp(k2) + disp(-1k + k1 + k2) + disp(-1k1 + q)) ^
                        -1,
                    ) + ifelse(
                        and(iszero(-1mu - v + disp(-1k1 + q)), iszero(-1k - k1 + q)),
                        u ^ 2 * dfermidist(-1mu + disp(k2), beta),
                        (
                            u ^ 2 * (
                                -fermidist(-1mu + disp(k2), beta) +
                                fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
                            )
                        ) *
                        (mu + v - disp(k2) - disp(-1k1 + q) + disp(-1k - k1 + k2 + q)) ^ -1,
                    )
                ) *
                (
                    ifelse(
                        and(iszero(-1mu - vp + disp(-1k1 + q)), iszero(-1k1 - kp + q)),
                        u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                        (
                            u ^ 2 * (
                                -fermidist(-1mu + disp(k3), beta) +
                                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                            )
                        ) *
                        (mu + vp - disp(k3) - disp(-1k1 + q) + disp(-1k1 + k3 - kp + q)) ^
                        -1,
                    ) + ifelse(
                        and(iszero(-1mu + vp - w + disp(-1k1 + q)), iszero(-1k1 + kp)),
                        -(u ^ 2 * dfermidist(-1mu + disp(k3), beta)),
                        (
                            u ^ 2 * (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                            )
                        ) *
                        (mu - vp + w - disp(k3) + disp(-1k1 + k3 + kp) - disp(-1k1 + q)) ^
                        -1,
                    )
                )
            ) * (2 * (2mu + w - disp(k1) - disp(-1k1 + q))) ^ -1
        )
    )
end
