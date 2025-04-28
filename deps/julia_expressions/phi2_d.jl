function phi2_d(
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
            fermidist(-1mu + disp(k1), beta) *
            (
                mu * u + u * vp + u * w - (u * disp(k1 - k3 + kp + q)) -
                (cnorm * u ^ 2 * fermidist(-1mu + disp(k3), beta)) +
                cnorm * u ^ 2 * fermidist(mu - disp(k1 - k3 + kp + q), beta) +
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
                ) - (
                    2 *
                    cnorm *
                    mu *
                    ifelse(
                        and(iszero(mu + vp - disp(k1)), iszero(-1k1 + kp)),
                        u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                        (
                            u ^ 2 * (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                            )
                        ) * (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) ^ -1,
                    )
                ) - (
                    2 *
                    cnorm *
                    vp *
                    ifelse(
                        and(iszero(mu + vp - disp(k1)), iszero(-1k1 + kp)),
                        u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                        (
                            u ^ 2 * (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                            )
                        ) * (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) ^ -1,
                    )
                ) - (
                    2 *
                    cnorm *
                    w *
                    ifelse(
                        and(iszero(mu + vp - disp(k1)), iszero(-1k1 + kp)),
                        u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                        (
                            u ^ 2 * (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                            )
                        ) * (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) ^ -1,
                    )
                ) +
                2 *
                cnorm *
                disp(k1 - k3 + kp + q) *
                ifelse(
                    and(iszero(mu + vp - disp(k1)), iszero(-1k1 + kp)),
                    u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                    (
                        u ^ 2 * (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                        )
                    ) * (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) ^ -1,
                ) +
                disp(k1) * (
                    u +
                    cnorm * ifelse(
                        and(iszero(w), iszero(q)),
                        u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                        (
                            u ^ 2 * (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(k3 + q), beta)
                            )
                        ) * (w + disp(k3) - disp(k3 + q)) ^ -1,
                    ) - (
                        2 *
                        cnorm *
                        ifelse(
                            and(iszero(mu + vp - disp(k1)), iszero(-1k1 + kp)),
                            u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                            (
                                u ^ 2 * (
                                    fermidist(-1mu + disp(k3), beta) -
                                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                                )
                            ) * (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) ^ -1,
                        )
                    )
                ) - (
                    disp(k3) * (
                        u +
                        cnorm * ifelse(
                            and(iszero(w), iszero(q)),
                            u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                            (
                                u ^ 2 * (
                                    fermidist(-1mu + disp(k3), beta) -
                                    fermidist(-1mu + disp(k3 + q), beta)
                                )
                            ) * (w + disp(k3) - disp(k3 + q)) ^ -1,
                        ) - (
                            2 *
                            cnorm *
                            ifelse(
                                and(iszero(mu + vp - disp(k1)), iszero(-1k1 + kp)),
                                u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                                (
                                    u ^ 2 * (
                                        fermidist(-1mu + disp(k3), beta) -
                                        fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                                    )
                                ) *
                                (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) ^ -1,
                            )
                        )
                    )
                )
            ) *
            (
                mu * u + u * v + u * w - (u * disp(k + k1 - k2 + q)) -
                (cnorm * u ^ 2 * fermidist(-1mu + disp(k2), beta)) +
                cnorm * u ^ 2 * fermidist(mu - disp(k + k1 - k2 + q), beta) - (
                    2 *
                    cnorm *
                    mu *
                    ifelse(
                        and(iszero(-1mu - v + disp(k1)), iszero(-1k + k1)),
                        u ^ 2 * dfermidist(-1mu + disp(k2), beta),
                        (
                            u ^ 2 * (
                                -fermidist(-1mu + disp(k2), beta) +
                                fermidist(-1mu + disp(-1k + k1 + k2), beta)
                            )
                        ) * (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) ^ -1,
                    )
                ) - (
                    2 *
                    cnorm *
                    v *
                    ifelse(
                        and(iszero(-1mu - v + disp(k1)), iszero(-1k + k1)),
                        u ^ 2 * dfermidist(-1mu + disp(k2), beta),
                        (
                            u ^ 2 * (
                                -fermidist(-1mu + disp(k2), beta) +
                                fermidist(-1mu + disp(-1k + k1 + k2), beta)
                            )
                        ) * (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) ^ -1,
                    )
                ) - (
                    2 *
                    cnorm *
                    w *
                    ifelse(
                        and(iszero(-1mu - v + disp(k1)), iszero(-1k + k1)),
                        u ^ 2 * dfermidist(-1mu + disp(k2), beta),
                        (
                            u ^ 2 * (
                                -fermidist(-1mu + disp(k2), beta) +
                                fermidist(-1mu + disp(-1k + k1 + k2), beta)
                            )
                        ) * (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) ^ -1,
                    )
                ) +
                2 *
                cnorm *
                disp(k + k1 - k2 + q) *
                ifelse(
                    and(iszero(-1mu - v + disp(k1)), iszero(-1k + k1)),
                    u ^ 2 * dfermidist(-1mu + disp(k2), beta),
                    (
                        u ^ 2 * (
                            -fermidist(-1mu + disp(k2), beta) +
                            fermidist(-1mu + disp(-1k + k1 + k2), beta)
                        )
                    ) * (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) ^ -1,
                ) +
                disp(k1) * (
                    u - (
                        2 *
                        cnorm *
                        ifelse(
                            and(iszero(-1mu - v + disp(k1)), iszero(-1k + k1)),
                            u ^ 2 * dfermidist(-1mu + disp(k2), beta),
                            (
                                u ^ 2 * (
                                    -fermidist(-1mu + disp(k2), beta) +
                                    fermidist(-1mu + disp(-1k + k1 + k2), beta)
                                )
                            ) * (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) ^ -1,
                        )
                    )
                ) - (
                    disp(k2) * (
                        u - (
                            2 *
                            cnorm *
                            ifelse(
                                and(iszero(-1mu - v + disp(k1)), iszero(-1k + k1)),
                                u ^ 2 * dfermidist(-1mu + disp(k2), beta),
                                (
                                    u ^ 2 * (
                                        -fermidist(-1mu + disp(k2), beta) +
                                        fermidist(-1mu + disp(-1k + k1 + k2), beta)
                                    )
                                ) *
                                (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) ^ -1,
                            )
                        )
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
                fermidist(-1mu + disp(k1 + q), beta) *
                (
                    mu * u + u * vp - (u * disp(k1 - k3 + kp + q)) -
                    (cnorm * u ^ 2 * fermidist(-1mu + disp(k3), beta)) +
                    cnorm * u ^ 2 * fermidist(mu - disp(k1 - k3 + kp + q), beta) +
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
                    ) - (
                        2 *
                        cnorm *
                        mu *
                        ifelse(
                            and(iszero(mu + vp + w - disp(k1 + q)), iszero(-1k1 + kp)),
                            u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                            (
                                u ^ 2 * (
                                    fermidist(-1mu + disp(k3), beta) -
                                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                                )
                            ) *
                            (mu + vp + w + disp(k3) - disp(-1k1 + k3 + kp) - disp(k1 + q)) ^
                            -1,
                        )
                    ) - (
                        2 *
                        cnorm *
                        vp *
                        ifelse(
                            and(iszero(mu + vp + w - disp(k1 + q)), iszero(-1k1 + kp)),
                            u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                            (
                                u ^ 2 * (
                                    fermidist(-1mu + disp(k3), beta) -
                                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                                )
                            ) *
                            (mu + vp + w + disp(k3) - disp(-1k1 + k3 + kp) - disp(k1 + q)) ^
                            -1,
                        )
                    ) +
                    2 *
                    cnorm *
                    disp(k1 - k3 + kp + q) *
                    ifelse(
                        and(iszero(mu + vp + w - disp(k1 + q)), iszero(-1k1 + kp)),
                        u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                        (
                            u ^ 2 * (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                            )
                        ) *
                        (mu + vp + w + disp(k3) - disp(-1k1 + k3 + kp) - disp(k1 + q)) ^ -1,
                    ) - (
                        disp(k3) * (
                            u +
                            cnorm * ifelse(
                                and(iszero(w), iszero(q)),
                                u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                                (
                                    u ^ 2 * (
                                        fermidist(-1mu + disp(k3), beta) -
                                        fermidist(-1mu + disp(k3 + q), beta)
                                    )
                                ) * (w + disp(k3) - disp(k3 + q)) ^ -1,
                            ) - (
                                2 *
                                cnorm *
                                ifelse(
                                    and(
                                        iszero(mu + vp + w - disp(k1 + q)),
                                        iszero(-1k1 + kp),
                                    ),
                                    u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                                    (
                                        u ^ 2 * (
                                            fermidist(-1mu + disp(k3), beta) -
                                            fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                                        )
                                    ) *
                                    (
                                        mu + vp + w + disp(k3) - disp(-1k1 + k3 + kp) -
                                        disp(k1 + q)
                                    ) ^ -1,
                                )
                            )
                        )
                    ) +
                    disp(k1 + q) * (
                        u +
                        cnorm * ifelse(
                            and(iszero(w), iszero(q)),
                            u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                            (
                                u ^ 2 * (
                                    fermidist(-1mu + disp(k3), beta) -
                                    fermidist(-1mu + disp(k3 + q), beta)
                                )
                            ) * (w + disp(k3) - disp(k3 + q)) ^ -1,
                        ) - (
                            2 *
                            cnorm *
                            ifelse(
                                and(iszero(mu + vp + w - disp(k1 + q)), iszero(-1k1 + kp)),
                                u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                                (
                                    u ^ 2 * (
                                        fermidist(-1mu + disp(k3), beta) -
                                        fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                                    )
                                ) *
                                (
                                    mu + vp + w + disp(k3) - disp(-1k1 + k3 + kp) -
                                    disp(k1 + q)
                                ) ^ -1,
                            )
                        )
                    )
                ) *
                (
                    mu * u + u * v - (u * disp(k + k1 - k2 + q)) -
                    (cnorm * u ^ 2 * fermidist(-1mu + disp(k2), beta)) +
                    cnorm * u ^ 2 * fermidist(mu - disp(k + k1 - k2 + q), beta) - (
                        2 *
                        cnorm *
                        mu *
                        ifelse(
                            and(iszero(-1mu - v - w + disp(k1 + q)), iszero(-1k + k1)),
                            u ^ 2 * dfermidist(-1mu + disp(k2), beta),
                            (
                                u ^ 2 * (
                                    -fermidist(-1mu + disp(k2), beta) +
                                    fermidist(-1mu + disp(-1k + k1 + k2), beta)
                                )
                            ) *
                            (mu + v + w - disp(k2) + disp(-1k + k1 + k2) - disp(k1 + q)) ^
                            -1,
                        )
                    ) - (
                        2 *
                        cnorm *
                        v *
                        ifelse(
                            and(iszero(-1mu - v - w + disp(k1 + q)), iszero(-1k + k1)),
                            u ^ 2 * dfermidist(-1mu + disp(k2), beta),
                            (
                                u ^ 2 * (
                                    -fermidist(-1mu + disp(k2), beta) +
                                    fermidist(-1mu + disp(-1k + k1 + k2), beta)
                                )
                            ) *
                            (mu + v + w - disp(k2) + disp(-1k + k1 + k2) - disp(k1 + q)) ^
                            -1,
                        )
                    ) +
                    2 *
                    cnorm *
                    disp(k + k1 - k2 + q) *
                    ifelse(
                        and(iszero(-1mu - v - w + disp(k1 + q)), iszero(-1k + k1)),
                        u ^ 2 * dfermidist(-1mu + disp(k2), beta),
                        (
                            u ^ 2 * (
                                -fermidist(-1mu + disp(k2), beta) +
                                fermidist(-1mu + disp(-1k + k1 + k2), beta)
                            )
                        ) *
                        (mu + v + w - disp(k2) + disp(-1k + k1 + k2) - disp(k1 + q)) ^ -1,
                    ) - (
                        disp(k2) * (
                            u - (
                                2 *
                                cnorm *
                                ifelse(
                                    and(
                                        iszero(-1mu - v - w + disp(k1 + q)),
                                        iszero(-1k + k1),
                                    ),
                                    u ^ 2 * dfermidist(-1mu + disp(k2), beta),
                                    (
                                        u ^ 2 * (
                                            -fermidist(-1mu + disp(k2), beta) +
                                            fermidist(-1mu + disp(-1k + k1 + k2), beta)
                                        )
                                    ) *
                                    (
                                        mu + v + w - disp(k2) + disp(-1k + k1 + k2) -
                                        disp(k1 + q)
                                    ) ^ -1,
                                )
                            )
                        )
                    ) +
                    disp(k1 + q) * (
                        u - (
                            2 *
                            cnorm *
                            ifelse(
                                and(iszero(-1mu - v - w + disp(k1 + q)), iszero(-1k + k1)),
                                u ^ 2 * dfermidist(-1mu + disp(k2), beta),
                                (
                                    u ^ 2 * (
                                        -fermidist(-1mu + disp(k2), beta) +
                                        fermidist(-1mu + disp(-1k + k1 + k2), beta)
                                    )
                                ) *
                                (
                                    mu + v + w - disp(k2) + disp(-1k + k1 + k2) -
                                    disp(k1 + q)
                                ) ^ -1,
                            )
                        )
                    )
                )
            ) *
            (
                (w + disp(k1) - disp(k1 + q)) *
                (mu + v - disp(k2) + disp(k1 + q) - disp(k + k1 - k2 + q)) *
                (mu + vp - disp(k3) + disp(k1 + q) - disp(k1 - k3 + kp + q))
            ) ^ -1
        ) - (
            (
                cnorm *
                u ^ 2 *
                (
                    fermidist(-1mu + disp(k2), beta) -
                    fermidist(mu - disp(k + k1 - k2 + q), beta)
                ) *
                fermidist(-2mu + disp(k2) + disp(k + k1 - k2 + q), beta) *
                (
                    u * v - (u * vp) - (u * disp(k + k1 - k2 + q)) +
                    u * disp(k1 - k3 + kp + q) +
                    cnorm * u ^ 2 * fermidist(-1mu + disp(k3), beta) -
                    (cnorm * u ^ 2 * fermidist(mu - disp(k1 - k3 + kp + q), beta)) +
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
                    ) - (
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
                        )
                    ) - (
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
                        )
                    ) +
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
                    ) - (
                        2 *
                        cnorm *
                        v *
                        ifelse(
                            and(
                                iszero(2mu + v + vp + w - disp(k2) - disp(k + k1 - k2 + q)),
                                iszero(-1k1 + kp),
                            ),
                            u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                            (
                                u ^ 2 * (
                                    fermidist(-1mu + disp(k3), beta) -
                                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                                )
                            ) *
                            (
                                2mu + v + vp + w - disp(k2) + disp(k3) -
                                disp(-1k1 + k3 + kp) - disp(k + k1 - k2 + q)
                            ) ^ -1,
                        )
                    ) +
                    2 *
                    cnorm *
                    vp *
                    ifelse(
                        and(
                            iszero(2mu + v + vp + w - disp(k2) - disp(k + k1 - k2 + q)),
                            iszero(-1k1 + kp),
                        ),
                        u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                        (
                            u ^ 2 * (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                            )
                        ) *
                        (
                            2mu + v + vp + w - disp(k2) + disp(k3) - disp(-1k1 + k3 + kp) -
                            disp(k + k1 - k2 + q)
                        ) ^ -1,
                    ) +
                    2 *
                    cnorm *
                    disp(k + k1 - k2 + q) *
                    ifelse(
                        and(
                            iszero(2mu + v + vp + w - disp(k2) - disp(k + k1 - k2 + q)),
                            iszero(-1k1 + kp),
                        ),
                        u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                        (
                            u ^ 2 * (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                            )
                        ) *
                        (
                            2mu + v + vp + w - disp(k2) + disp(k3) - disp(-1k1 + k3 + kp) -
                            disp(k + k1 - k2 + q)
                        ) ^ -1,
                    ) - (
                        2 *
                        cnorm *
                        disp(k1 - k3 + kp + q) *
                        ifelse(
                            and(
                                iszero(2mu + v + vp + w - disp(k2) - disp(k + k1 - k2 + q)),
                                iszero(-1k1 + kp),
                            ),
                            u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                            (
                                u ^ 2 * (
                                    fermidist(-1mu + disp(k3), beta) -
                                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                                )
                            ) *
                            (
                                2mu + v + vp + w - disp(k2) + disp(k3) -
                                disp(-1k1 + k3 + kp) - disp(k + k1 - k2 + q)
                            ) ^ -1,
                        )
                    ) - (
                        disp(k2) * (
                            u +
                            cnorm * ifelse(
                                and(iszero(w), iszero(q)),
                                u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                                (
                                    u ^ 2 * (
                                        fermidist(-1mu + disp(k3), beta) -
                                        fermidist(-1mu + disp(k3 + q), beta)
                                    )
                                ) * (w + disp(k3) - disp(k3 + q)) ^ -1,
                            ) - (
                                2 *
                                cnorm *
                                ifelse(
                                    and(
                                        iszero(
                                            2mu + v + vp + w - disp(k2) -
                                            disp(k + k1 - k2 + q),
                                        ),
                                        iszero(-1k1 + kp),
                                    ),
                                    u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                                    (
                                        u ^ 2 * (
                                            fermidist(-1mu + disp(k3), beta) -
                                            fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                                        )
                                    ) *
                                    (
                                        2mu + v + vp + w - disp(k2) + disp(k3) -
                                        disp(-1k1 + k3 + kp) - disp(k + k1 - k2 + q)
                                    ) ^ -1,
                                )
                            )
                        )
                    ) +
                    disp(k3) * (
                        u +
                        cnorm * ifelse(
                            and(iszero(w), iszero(q)),
                            u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                            (
                                u ^ 2 * (
                                    fermidist(-1mu + disp(k3), beta) -
                                    fermidist(-1mu + disp(k3 + q), beta)
                                )
                            ) * (w + disp(k3) - disp(k3 + q)) ^ -1,
                        ) - (
                            2 *
                            cnorm *
                            ifelse(
                                and(
                                    iszero(
                                        2mu + v + vp + w - disp(k2) - disp(k + k1 - k2 + q),
                                    ),
                                    iszero(-1k1 + kp),
                                ),
                                u ^ 2 * dfermidist(-1mu + disp(k3), beta),
                                (
                                    u ^ 2 * (
                                        fermidist(-1mu + disp(k3), beta) -
                                        fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                                    )
                                ) *
                                (
                                    2mu + v + vp + w - disp(k2) + disp(k3) -
                                    disp(-1k1 + k3 + kp) - disp(k + k1 - k2 + q)
                                ) ^ -1,
                            )
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
                (-1 + 2 * fermidist(-2mu + disp(k2) + disp(k + k1 - k2 + q), beta))
            ) ^ -1
        ) +
        (
            cnorm *
            u ^ 2 *
            (
                fermidist(-1mu + disp(k3), beta) -
                fermidist(mu - disp(k1 - k3 + kp + q), beta)
            ) *
            fermidist(-2mu + disp(k3) + disp(k1 - k3 + kp + q), beta) *
            (
                u * v - (u * vp) - (u * disp(k + k1 - k2 + q)) +
                u * disp(k1 - k3 + kp + q) -
                (cnorm * u ^ 2 * fermidist(-1mu + disp(k2), beta)) +
                cnorm * u ^ 2 * fermidist(mu - disp(k + k1 - k2 + q), beta) - (
                    2 *
                    cnorm *
                    v *
                    ifelse(
                        and(
                            iszero(-2mu - v - vp - w + disp(k3) + disp(k1 - k3 + kp + q)),
                            iszero(-1k + k1),
                        ),
                        u ^ 2 * dfermidist(-1mu + disp(k2), beta),
                        (
                            u ^ 2 * (
                                -fermidist(-1mu + disp(k2), beta) +
                                fermidist(-1mu + disp(-1k + k1 + k2), beta)
                            )
                        ) *
                        (
                            2mu + v + vp + w - disp(k2) + disp(-1k + k1 + k2) - disp(k3) -
                            disp(k1 - k3 + kp + q)
                        ) ^ -1,
                    )
                ) +
                2 *
                cnorm *
                vp *
                ifelse(
                    and(
                        iszero(-2mu - v - vp - w + disp(k3) + disp(k1 - k3 + kp + q)),
                        iszero(-1k + k1),
                    ),
                    u ^ 2 * dfermidist(-1mu + disp(k2), beta),
                    (
                        u ^ 2 * (
                            -fermidist(-1mu + disp(k2), beta) +
                            fermidist(-1mu + disp(-1k + k1 + k2), beta)
                        )
                    ) *
                    (
                        2mu + v + vp + w - disp(k2) + disp(-1k + k1 + k2) - disp(k3) -
                        disp(k1 - k3 + kp + q)
                    ) ^ -1,
                ) +
                2 *
                cnorm *
                disp(k + k1 - k2 + q) *
                ifelse(
                    and(
                        iszero(-2mu - v - vp - w + disp(k3) + disp(k1 - k3 + kp + q)),
                        iszero(-1k + k1),
                    ),
                    u ^ 2 * dfermidist(-1mu + disp(k2), beta),
                    (
                        u ^ 2 * (
                            -fermidist(-1mu + disp(k2), beta) +
                            fermidist(-1mu + disp(-1k + k1 + k2), beta)
                        )
                    ) *
                    (
                        2mu + v + vp + w - disp(k2) + disp(-1k + k1 + k2) - disp(k3) -
                        disp(k1 - k3 + kp + q)
                    ) ^ -1,
                ) - (
                    2 *
                    cnorm *
                    disp(k1 - k3 + kp + q) *
                    ifelse(
                        and(
                            iszero(-2mu - v - vp - w + disp(k3) + disp(k1 - k3 + kp + q)),
                            iszero(-1k + k1),
                        ),
                        u ^ 2 * dfermidist(-1mu + disp(k2), beta),
                        (
                            u ^ 2 * (
                                -fermidist(-1mu + disp(k2), beta) +
                                fermidist(-1mu + disp(-1k + k1 + k2), beta)
                            )
                        ) *
                        (
                            2mu + v + vp + w - disp(k2) + disp(-1k + k1 + k2) - disp(k3) -
                            disp(k1 - k3 + kp + q)
                        ) ^ -1,
                    )
                ) - (
                    disp(k2) * (
                        u - (
                            2 *
                            cnorm *
                            ifelse(
                                and(
                                    iszero(
                                        -2mu - v - vp - w +
                                        disp(k3) +
                                        disp(k1 - k3 + kp + q),
                                    ),
                                    iszero(-1k + k1),
                                ),
                                u ^ 2 * dfermidist(-1mu + disp(k2), beta),
                                (
                                    u ^ 2 * (
                                        -fermidist(-1mu + disp(k2), beta) +
                                        fermidist(-1mu + disp(-1k + k1 + k2), beta)
                                    )
                                ) *
                                (
                                    2mu + v + vp + w - disp(k2) + disp(-1k + k1 + k2) -
                                    disp(k3) - disp(k1 - k3 + kp + q)
                                ) ^ -1,
                            )
                        )
                    )
                ) +
                disp(k3) * (
                    u - (
                        2 *
                        cnorm *
                        ifelse(
                            and(
                                iszero(
                                    -2mu - v - vp - w + disp(k3) + disp(k1 - k3 + kp + q),
                                ),
                                iszero(-1k + k1),
                            ),
                            u ^ 2 * dfermidist(-1mu + disp(k2), beta),
                            (
                                u ^ 2 * (
                                    -fermidist(-1mu + disp(k2), beta) +
                                    fermidist(-1mu + disp(-1k + k1 + k2), beta)
                                )
                            ) *
                            (
                                2mu + v + vp + w - disp(k2) + disp(-1k + k1 + k2) -
                                disp(k3) - disp(k1 - k3 + kp + q)
                            ) ^ -1,
                        )
                    )
                )
            )
        ) *
        (
            (mu + vp + w + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q)) *
            (mu + vp - disp(k3) + disp(k1 + q) - disp(k1 - k3 + kp + q)) *
            (
                -1v + vp + disp(k2) - disp(k3) + disp(k + k1 - k2 + q) -
                disp(k1 - k3 + kp + q)
            ) *
            (-1 + 2 * fermidist(-2mu + disp(k3) + disp(k1 - k3 + kp + q), beta))
        ) ^ -1
    )
end
