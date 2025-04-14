function phi2_m(
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
            u *
            (
                -1v + vp + disp(k2) - disp(k3) + disp(k + k1 - k2 + q) -
                disp(k1 - k3 + kp + q) + cnorm * u * fermidist(-1mu + disp(k2), beta) -
                (cnorm * u * fermidist(mu - disp(k + k1 - k2 + q), beta))
            ) *
            (
                cnorm * u ^ 2 * fermidist(-1mu + disp(k3), beta) -
                (cnorm * u ^ 2 * fermidist(mu - disp(k1 - k3 + kp + q), beta))
            ) *
            fermidist(-2mu + disp(k3) + disp(k1 - k3 + kp + q), beta)
        ) *
        (
            (mu + vp + w + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q)) *
            (mu + vp - disp(k3) + disp(k1 + q) - disp(k1 - k3 + kp + q)) *
            (
                -1v + vp + disp(k2) - disp(k3) + disp(k + k1 - k2 + q) -
                disp(k1 - k3 + kp + q)
            ) *
            (1 - (2 * fermidist(-2mu + disp(k3) + disp(k1 - k3 + kp + q), beta)))
        ) ^ -1 +
        (
            u *
            fermidist(-1mu + disp(k1), beta) *
            (
                -1mu - v - w - disp(k1) +
                disp(k2) +
                disp(k + k1 - k2 + q) +
                cnorm * u * fermidist(-1mu + disp(k2), beta) -
                (cnorm * u * fermidist(mu - disp(k + k1 - k2 + q), beta))
            ) *
            (
                -(mu * u) - (u * vp) - (u * w) - (u * disp(k1)) +
                u * disp(k3) +
                u * disp(k1 - k3 + kp + q) +
                cnorm * u ^ 2 * fermidist(-1mu + disp(k3), beta) -
                (cnorm * u ^ 2 * fermidist(mu - disp(k1 - k3 + kp + q), beta)) +
                cnorm *
                mu *
                ifelse(
                    and(iszero(w), iszero(q)),
                    u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                    (
                        u ^ 2 * (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(-1mu + disp(k3 + q), beta)
                        )
                    ) * (w + disp(k3) - disp(k3 + q)) ^ -1,
                ) +
                cnorm *
                vp *
                ifelse(
                    and(iszero(w), iszero(q)),
                    u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                    (
                        u ^ 2 * (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(-1mu + disp(k3 + q), beta)
                        )
                    ) * (w + disp(k3) - disp(k3 + q)) ^ -1,
                ) +
                cnorm *
                w *
                ifelse(
                    and(iszero(w), iszero(q)),
                    u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                    (
                        u ^ 2 * (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(-1mu + disp(k3 + q), beta)
                        )
                    ) * (w + disp(k3) - disp(k3 + q)) ^ -1,
                ) +
                cnorm *
                disp(k1) *
                ifelse(
                    and(iszero(w), iszero(q)),
                    u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                    (
                        u ^ 2 * (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(-1mu + disp(k3 + q), beta)
                        )
                    ) * (w + disp(k3) - disp(k3 + q)) ^ -1,
                ) - (
                    cnorm *
                    disp(k3) *
                    ifelse(
                        and(iszero(w), iszero(q)),
                        u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                        (
                            u ^ 2 * (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(k3 + q), beta)
                            )
                        ) * (w + disp(k3) - disp(k3 + q)) ^ -1,
                    )
                ) - (
                    cnorm *
                    disp(k1 - k3 + kp + q) *
                    ifelse(
                        and(iszero(w), iszero(q)),
                        u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                        (
                            u ^ 2 * (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(k3 + q), beta)
                            )
                        ) * (w + disp(k3) - disp(k3 + q)) ^ -1,
                    )
                )
            )
        ) *
        (
            (w + disp(k1) - disp(k1 + q)) *
            (mu + v + w + disp(k1) - disp(k2) - disp(k + k1 - k2 + q)) *
            (mu + vp + w + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q))
        ) ^ -1 - (
            (
                u *
                fermidist(-1mu + disp(k1 + q), beta) *
                (
                    -1mu - v + disp(k2) - disp(k1 + q) +
                    disp(k + k1 - k2 + q) +
                    cnorm * u * fermidist(-1mu + disp(k2), beta) -
                    (cnorm * u * fermidist(mu - disp(k + k1 - k2 + q), beta))
                ) *
                (
                    -(mu * u) - (u * vp) + u * disp(k3) - (u * disp(k1 + q)) +
                    u * disp(k1 - k3 + kp + q) +
                    cnorm * u ^ 2 * fermidist(-1mu + disp(k3), beta) -
                    (cnorm * u ^ 2 * fermidist(mu - disp(k1 - k3 + kp + q), beta)) +
                    cnorm *
                    mu *
                    ifelse(
                        and(iszero(w), iszero(q)),
                        u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                        (
                            u ^ 2 * (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(k3 + q), beta)
                            )
                        ) * (w + disp(k3) - disp(k3 + q)) ^ -1,
                    ) +
                    cnorm *
                    vp *
                    ifelse(
                        and(iszero(w), iszero(q)),
                        u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                        (
                            u ^ 2 * (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(k3 + q), beta)
                            )
                        ) * (w + disp(k3) - disp(k3 + q)) ^ -1,
                    ) - (
                        cnorm *
                        disp(k3) *
                        ifelse(
                            and(iszero(w), iszero(q)),
                            u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                            (
                                u ^ 2 * (
                                    fermidist(-1mu + disp(k3), beta) -
                                    fermidist(-1mu + disp(k3 + q), beta)
                                )
                            ) * (w + disp(k3) - disp(k3 + q)) ^ -1,
                        )
                    ) +
                    cnorm *
                    disp(k1 + q) *
                    ifelse(
                        and(iszero(w), iszero(q)),
                        u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                        (
                            u ^ 2 * (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(k3 + q), beta)
                            )
                        ) * (w + disp(k3) - disp(k3 + q)) ^ -1,
                    ) - (
                        cnorm *
                        disp(k1 - k3 + kp + q) *
                        ifelse(
                            and(iszero(w), iszero(q)),
                            u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                            (
                                u ^ 2 * (
                                    fermidist(-1mu + disp(k3), beta) -
                                    fermidist(-1mu + disp(k3 + q), beta)
                                )
                            ) * (w + disp(k3) - disp(k3 + q)) ^ -1,
                        )
                    )
                )
            ) *
            (
                (w + disp(k1) - disp(k1 + q)) *
                (mu + v - disp(k2) + disp(k1 + q) - disp(k + k1 - k2 + q)) *
                (mu + vp - disp(k3) + disp(k1 + q) - disp(k1 - k3 + kp + q))
            ) ^ -1
        ) +
        (
            u *
            (
                cnorm * u * fermidist(-1mu + disp(k2), beta) -
                (cnorm * u * fermidist(mu - disp(k + k1 - k2 + q), beta))
            ) *
            fermidist(-2mu + disp(k2) + disp(k + k1 - k2 + q), beta) *
            (
                u * v - (u * vp) - (u * disp(k2)) + u * disp(k3) -
                (u * disp(k + k1 - k2 + q)) +
                u * disp(k1 - k3 + kp + q) +
                cnorm * u ^ 2 * fermidist(-1mu + disp(k3), beta) -
                (cnorm * u ^ 2 * fermidist(mu - disp(k1 - k3 + kp + q), beta)) - (
                    cnorm *
                    v *
                    ifelse(
                        and(iszero(w), iszero(q)),
                        u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                        (
                            u ^ 2 * (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(k3 + q), beta)
                            )
                        ) * (w + disp(k3) - disp(k3 + q)) ^ -1,
                    )
                ) +
                cnorm *
                vp *
                ifelse(
                    and(iszero(w), iszero(q)),
                    u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                    (
                        u ^ 2 * (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(-1mu + disp(k3 + q), beta)
                        )
                    ) * (w + disp(k3) - disp(k3 + q)) ^ -1,
                ) +
                cnorm *
                disp(k2) *
                ifelse(
                    and(iszero(w), iszero(q)),
                    u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                    (
                        u ^ 2 * (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(-1mu + disp(k3 + q), beta)
                        )
                    ) * (w + disp(k3) - disp(k3 + q)) ^ -1,
                ) - (
                    cnorm *
                    disp(k3) *
                    ifelse(
                        and(iszero(w), iszero(q)),
                        u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                        (
                            u ^ 2 * (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(k3 + q), beta)
                            )
                        ) * (w + disp(k3) - disp(k3 + q)) ^ -1,
                    )
                ) +
                cnorm *
                disp(k + k1 - k2 + q) *
                ifelse(
                    and(iszero(w), iszero(q)),
                    u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                    (
                        u ^ 2 * (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(-1mu + disp(k3 + q), beta)
                        )
                    ) * (w + disp(k3) - disp(k3 + q)) ^ -1,
                ) - (
                    cnorm *
                    disp(k1 - k3 + kp + q) *
                    ifelse(
                        and(iszero(w), iszero(q)),
                        u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                        (
                            u ^ 2 * (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(k3 + q), beta)
                            )
                        ) * (w + disp(k3) - disp(k3 + q)) ^ -1,
                    )
                )
            )
        ) *
        (
            (mu + v + w + disp(k1) - disp(k2) - disp(k + k1 - k2 + q)) *
            (mu + v - disp(k2) + disp(k1 + q) - disp(k + k1 - k2 + q)) *
            (
                v - vp - disp(k2) + disp(k3) - disp(k + k1 - k2 + q) +
                disp(k1 - k3 + kp + q)
            ) *
            (1 - (2 * fermidist(-2mu + disp(k2) + disp(k + k1 - k2 + q), beta)))
        ) ^ -1
    )
end
