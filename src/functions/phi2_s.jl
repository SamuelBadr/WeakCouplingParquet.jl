@cse function phi2_s(
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
    atol = 1e-6 #eps()
    cnorm = 1 / (2pi)^length(k)
    ifelse(
        and(
            isless(
                abs2(
                    v - vp - disp(k2) + disp(-1k + k1 + k2) - disp(k3) +
                    disp(-1k1 + k3 + kp),
                ),
                atol,
            ),
            isless(abs2(2mu + w - disp(k1) - disp(-1k1 + q)), atol),
        ),
        (
            cnorm ^ 3 *
            u ^ 4 *
            (
                -(
                    (
                        dfermidist(-1mu + disp(k1), beta) *
                        (mu + v - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) +
                        2 * fermidist(-1mu + disp(k1), beta)
                    ) *
                    (
                        fermidist(-1mu + disp(k2), beta) -
                        fermidist(-1mu + disp(k2) + disp(k3) - disp(-1k1 + k3 + kp), beta)
                    ) *
                    (
                        fermidist(-1mu + disp(k3), beta) -
                        fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                    )
                ) - (
                    (
                        (mu + v - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) * (
                            -(
                                dfermidist(
                                    -1mu + disp(k2) + disp(k3) - disp(-1k1 + k3 + kp),
                                    beta,
                                ) * fermidist(v + disp(k3) - disp(-1k1 + k3 + kp), beta)
                            ) +
                            dfermidist(v + disp(k3) - disp(-1k1 + k3 + kp), beta) * (
                                fermidist(-1mu + disp(k2), beta) - fermidist(
                                    -1mu + disp(k2) + disp(k3) - disp(-1k1 + k3 + kp),
                                    beta,
                                )
                            )
                        ) - (
                            2 *
                            (-1 + fermidist(-1mu + disp(k2), beta)) *
                            fermidist(
                                -1mu + disp(k2) + disp(k3) - disp(-1k1 + k3 + kp),
                                beta,
                            )
                        )
                    ) * (
                        fermidist(-1mu + disp(k3), beta) -
                        fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                    )
                ) +
                dfermidist(-1mu + disp(k2) + disp(k3) - disp(-1k1 + k3 + kp), beta) *
                (mu + v - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) *
                fermidist(-1mu + disp(k3), beta) *
                (-1 + fermidist(-1mu + disp(-1k1 + k3 + kp), beta))
            )
        ) * (2 * (mu + v - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) ^ 3) ^ -1,
        ifelse(
            isless(abs2(2mu + w - disp(k1) - disp(-1k1 + q)), atol),
            (-2 ^ -1) * (
                (
                    cnorm ^ 3 *
                    u ^ 4 *
                    (
                        (
                            v - vp - disp(k2) + disp(-1k + k1 + k2) - disp(k3) +
                            disp(-1k1 + k3 + kp)
                        ) *
                        (
                            (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) *
                            fermidist(-1mu + disp(k1), beta) +
                            (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) * (
                                dfermidist(-1mu + disp(k1), beta) *
                                (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) +
                                fermidist(-1mu + disp(k1), beta)
                            )
                        ) *
                        (
                            fermidist(-1mu + disp(k2), beta) -
                            fermidist(-1mu + disp(-1k + k1 + k2), beta)
                        ) *
                        (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                        ) +
                        (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) ^ 2 *
                        (-1 + fermidist(-1mu + disp(k2), beta)) *
                        fermidist(-1mu + disp(-1k + k1 + k2), beta) *
                        (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                        ) +
                        (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) ^ 2 *
                        (
                            fermidist(-1mu + disp(k2), beta) -
                            fermidist(-1mu + disp(-1k + k1 + k2), beta)
                        ) *
                        fermidist(-1mu + disp(k3), beta) *
                        (-1 + fermidist(-1mu + disp(-1k1 + k3 + kp), beta))
                    )
                ) *
                (
                    (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) ^ 2 *
                    (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) ^ 2 *
                    (
                        v - vp - disp(k2) + disp(-1k + k1 + k2) - disp(k3) +
                        disp(-1k1 + k3 + kp)
                    )
                ) ^ -1
            ),
            ifelse(
                isless(
                    abs2(
                        v - vp - disp(k2) + disp(-1k + k1 + k2) - disp(k3) +
                        disp(-1k1 + k3 + kp),
                    ),
                    atol,
                ),
                (
                    cnorm ^ 3 *
                    u ^ 4 *
                    (
                        (
                            fermidist(-1mu + disp(k1), beta) *
                            (
                                fermidist(-1mu + disp(k2), beta) - fermidist(
                                    -1mu + disp(k2) + disp(k3) - disp(-1k1 + k3 + kp),
                                    beta,
                                )
                            ) *
                            (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                            )
                        ) *
                        (
                            (mu + v - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) ^ 2 *
                            (2mu + w - disp(k1) - disp(-1k1 + q))
                        ) ^ -1 +
                        (
                            (
                                -(
                                    (
                                        mu - v + w - disp(k3) + disp(-1k1 + k3 + kp) -
                                        disp(-1k1 + q)
                                    ) *
                                    (-1 + fermidist(-1mu + disp(k2), beta)) *
                                    fermidist(
                                        -1mu + disp(k2) + disp(k3) - disp(-1k1 + k3 + kp),
                                        beta,
                                    )
                                ) +
                                (mu + v - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) * (
                                    (
                                        mu - v + w - disp(k3) + disp(-1k1 + k3 + kp) -
                                        disp(-1k1 + q)
                                    ) * (
                                        -(
                                            dfermidist(
                                                -1mu + disp(k2) + disp(k3) -
                                                disp(-1k1 + k3 + kp),
                                                beta,
                                            ) * fermidist(
                                                v + disp(k3) - disp(-1k1 + k3 + kp),
                                                beta,
                                            )
                                        ) +
                                        dfermidist(
                                            v + disp(k3) - disp(-1k1 + k3 + kp),
                                            beta,
                                        ) * (
                                            fermidist(-1mu + disp(k2), beta) - fermidist(
                                                -1mu + disp(k2) + disp(k3) -
                                                disp(-1k1 + k3 + kp),
                                                beta,
                                            )
                                        )
                                    ) +
                                    (-1 + fermidist(-1mu + disp(k2), beta)) * fermidist(
                                        -1mu + disp(k2) + disp(k3) - disp(-1k1 + k3 + kp),
                                        beta,
                                    )
                                )
                            ) * (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                            )
                        ) *
                        (
                            (mu + v - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) ^ 2 *
                            (
                                mu - v + w - disp(k3) + disp(-1k1 + k3 + kp) -
                                disp(-1k1 + q)
                            ) ^ 2
                        ) ^ -1 - (
                            (
                                dfermidist(
                                    -1mu + disp(k2) + disp(k3) - disp(-1k1 + k3 + kp),
                                    beta,
                                ) *
                                fermidist(-1mu + disp(k3), beta) *
                                (-1 + fermidist(-1mu + disp(-1k1 + k3 + kp), beta))
                            ) *
                            (
                                (mu + v - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) * (
                                    mu - v + w - disp(k3) + disp(-1k1 + k3 + kp) -
                                    disp(-1k1 + q)
                                )
                            ) ^ -1
                        ) - (
                            (
                                (
                                    fermidist(-1mu + disp(k2), beta) - fermidist(
                                        -1mu + disp(k2) + disp(k3) - disp(-1k1 + k3 + kp),
                                        beta,
                                    )
                                ) *
                                (
                                    fermidist(-1mu + disp(k3), beta) -
                                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                                ) *
                                fermidist(mu + w - disp(-1k1 + q), beta)
                            ) *
                            (
                                (2mu + w - disp(k1) - disp(-1k1 + q)) *
                                (
                                    mu - v + w - disp(k3) + disp(-1k1 + k3 + kp) -
                                    disp(-1k1 + q)
                                ) ^ 2
                            ) ^ -1
                        )
                    )
                ) * 2 ^ -1,
                (
                    cnorm ^ 3 *
                    u ^ 4 *
                    (
                        (
                            fermidist(-1mu + disp(k1), beta) *
                            (
                                fermidist(-1mu + disp(k2), beta) -
                                fermidist(-1mu + disp(-1k + k1 + k2), beta)
                            ) *
                            (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                            )
                        ) *
                        (
                            (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) *
                            (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) *
                            (2mu + w - disp(k1) - disp(-1k1 + q))
                        ) ^ -1 +
                        (
                            (-1 + fermidist(-1mu + disp(k2), beta)) *
                            fermidist(-1mu + disp(-1k + k1 + k2), beta) *
                            (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                            )
                        ) *
                        (
                            (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) *
                            (
                                v - vp - disp(k2) + disp(-1k + k1 + k2) - disp(k3) +
                                disp(-1k1 + k3 + kp)
                            ) *
                            (mu - v + w + disp(k2) - disp(-1k + k1 + k2) - disp(-1k1 + q))
                        ) ^ -1 - (
                            (
                                (
                                    fermidist(-1mu + disp(k2), beta) -
                                    fermidist(-1mu + disp(-1k + k1 + k2), beta)
                                ) *
                                fermidist(-1mu + disp(k3), beta) *
                                (-1 + fermidist(-1mu + disp(-1k1 + k3 + kp), beta))
                            ) *
                            (
                                (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) *
                                (
                                    -1v + vp + disp(k2) - disp(-1k + k1 + k2) + disp(k3) -
                                    disp(-1k1 + k3 + kp)
                                ) *
                                (
                                    mu - vp + w - disp(k3) + disp(-1k1 + k3 + kp) -
                                    disp(-1k1 + q)
                                )
                            ) ^ -1
                        ) - (
                            (
                                (
                                    fermidist(-1mu + disp(k2), beta) -
                                    fermidist(-1mu + disp(-1k + k1 + k2), beta)
                                ) *
                                (
                                    fermidist(-1mu + disp(k3), beta) -
                                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                                ) *
                                fermidist(mu + w - disp(-1k1 + q), beta)
                            ) *
                            (
                                (2mu + w - disp(k1) - disp(-1k1 + q)) *
                                (
                                    mu - v + w + disp(k2) - disp(-1k + k1 + k2) -
                                    disp(-1k1 + q)
                                ) *
                                (
                                    mu - vp + w - disp(k3) + disp(-1k1 + k3 + kp) -
                                    disp(-1k1 + q)
                                )
                            ) ^ -1
                        )
                    )
                ) * 2 ^ -1,
            ),
        ),
    ) +
    ifelse(
        and(
            isless(abs2(2mu + w - disp(k1) - disp(-1k1 + q)), atol),
            isless(
                abs2(
                    v + vp - w - disp(k2) + disp(k3) - disp(-1k1 + k3 + kp) +
                    disp(-1k - k1 + k2 + q),
                ),
                atol,
            ),
        ),
        (
            cnorm ^ 3 *
            u ^ 4 *
            (
                (v + vp) *
                (
                    (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) *
                    fermidist(-1mu + disp(k1), beta) +
                    (mu - v - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) * (
                        dfermidist(-1mu + disp(k1), beta) *
                        (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) +
                        fermidist(-1mu + disp(k1), beta)
                    )
                ) *
                (
                    fermidist(-1mu + disp(k3), beta) -
                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                ) +
                (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) ^ 2 *
                fermidist(-1mu + disp(k3), beta) *
                (-1 + fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) - (
                    (-1mu + v + disp(k1) - disp(k3) + disp(-1k1 + k3 + kp)) ^ 2 *
                    fermidist(-1mu + disp(k3), beta) *
                    (-1 + fermidist(-1mu + disp(-1k1 + k3 + kp), beta))
                )
            ) *
            (
                fermidist(-1mu + disp(k2), beta) -
                fermidist(-1mu + disp(k2) - disp(k3) + disp(-1k1 + k3 + kp), beta)
            )
        ) *
        (
            2 *
            (v + vp) *
            (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) ^ 2 *
            (-1mu + v + disp(k1) - disp(k3) + disp(-1k1 + k3 + kp)) ^ 2
        ) ^ -1,
        ifelse(
            isless(
                abs2(
                    v + vp - w - disp(k2) + disp(k3) - disp(-1k1 + k3 + kp) +
                    disp(-1k - k1 + k2 + q),
                ),
                atol,
            ),
            (
                cnorm ^ 3 *
                u ^ 4 *
                (
                    -(
                        (
                            dfermidist(
                                -1mu + disp(k2) - disp(k3) + disp(-1k1 + k3 + kp),
                                beta,
                            ) *
                            fermidist(-1mu + disp(k3), beta) *
                            (-1 + fermidist(-1mu + disp(-1k1 + k3 + kp), beta))
                        ) *
                        (
                            (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) *
                            (mu + v - disp(k3) + disp(-1k1 + k3 + kp) - disp(-1k1 + q))
                        ) ^ -1
                    ) +
                    (
                        (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                        ) * (
                            (mu + v - disp(k3) + disp(-1k1 + k3 + kp) - disp(-1k1 + q)) *
                            (
                                -(
                                    (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) *
                                    (
                                        dfermidist(
                                            -1mu + disp(k2) - disp(k3) +
                                            disp(-1k1 + k3 + kp),
                                            beta,
                                        ) * fermidist(
                                            vp + disp(k3) - disp(-1k1 + k3 + kp),
                                            beta,
                                        ) +
                                        dfermidist(
                                            vp + disp(k3) - disp(-1k1 + k3 + kp),
                                            beta,
                                        ) * (
                                            fermidist(-1mu + disp(k2), beta) - fermidist(
                                                -1mu + disp(k2) - disp(k3) +
                                                disp(-1k1 + k3 + kp),
                                                beta,
                                            )
                                        )
                                    )
                                ) - (
                                    fermidist(-1mu + disp(k2), beta) * (
                                        -1 + fermidist(
                                            -1mu + disp(k2) - disp(k3) +
                                            disp(-1k1 + k3 + kp),
                                            beta,
                                        )
                                    )
                                )
                            ) +
                            (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) *
                            fermidist(-1mu + disp(k2), beta) *
                            (
                                -1 + fermidist(
                                    -1mu + disp(k2) - disp(k3) + disp(-1k1 + k3 + kp),
                                    beta,
                                )
                            )
                        )
                    ) *
                    (
                        (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) ^ 2 *
                        (mu + v - disp(k3) + disp(-1k1 + k3 + kp) - disp(-1k1 + q)) ^ 2
                    ) ^ -1 - (
                        (
                            fermidist(-1mu + disp(k1), beta) *
                            (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                            ) *
                            (
                                fermidist(-1mu + disp(k2), beta) - fermidist(
                                    -1mu + disp(k2) - disp(k3) + disp(-1k1 + k3 + kp),
                                    beta,
                                )
                            )
                        ) *
                        (
                            (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) ^ 2 *
                            (2mu + v + vp - disp(k1) - disp(-1k1 + q))
                        ) ^ -1
                    ) +
                    (
                        (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                        ) *
                        (
                            fermidist(-1mu + disp(k2), beta) - fermidist(
                                -1mu + disp(k2) - disp(k3) + disp(-1k1 + k3 + kp),
                                beta,
                            )
                        ) *
                        fermidist(mu + v + vp - disp(-1k1 + q), beta)
                    ) *
                    (
                        (2mu + v + vp - disp(k1) - disp(-1k1 + q)) *
                        (mu + v - disp(k3) + disp(-1k1 + k3 + kp) - disp(-1k1 + q)) ^ 2
                    ) ^ -1
                )
            ) * 2 ^ -1,
            ifelse(
                isless(abs2(2mu + w - disp(k1) - disp(-1k1 + q)), atol),
                (
                    cnorm ^ 3 *
                    u ^ 4 *
                    (
                        (
                            v + vp - disp(k2) + disp(k3) - disp(-1k1 + k3 + kp) +
                            disp(-1k - k1 + k2 + q)
                        ) *
                        (
                            (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) *
                            fermidist(-1mu + disp(k1), beta) +
                            (mu - v - disp(k1) + disp(k2) - disp(-1k - k1 + k2 + q)) * (
                                dfermidist(-1mu + disp(k1), beta) *
                                (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) +
                                fermidist(-1mu + disp(k1), beta)
                            )
                        ) *
                        (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                        ) *
                        (
                            fermidist(-1mu + disp(k2), beta) -
                            fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
                        ) - (
                            (-1mu + v + disp(k1) - disp(k2) + disp(-1k - k1 + k2 + q)) ^ 2 *
                            fermidist(-1mu + disp(k3), beta) *
                            (-1 + fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) *
                            (
                                fermidist(-1mu + disp(k2), beta) -
                                fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
                            )
                        ) +
                        (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) ^ 2 *
                        fermidist(-1mu + disp(k2), beta) *
                        (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                        ) *
                        (-1 + fermidist(-1mu + disp(-1k - k1 + k2 + q), beta))
                    )
                ) *
                (
                    2 *
                    (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) ^ 2 *
                    (-1mu + v + disp(k1) - disp(k2) + disp(-1k - k1 + k2 + q)) ^ 2 *
                    (
                        v + vp - disp(k2) + disp(k3) - disp(-1k1 + k3 + kp) +
                        disp(-1k - k1 + k2 + q)
                    )
                ) ^ -1,
                (
                    cnorm ^ 3 *
                    u ^ 4 *
                    (
                        -(
                            (
                                fermidist(-1mu + disp(k1), beta) *
                                (
                                    fermidist(-1mu + disp(k3), beta) -
                                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                                ) *
                                (
                                    fermidist(-1mu + disp(k2), beta) -
                                    fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
                                )
                            ) *
                            (
                                (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) *
                                (2mu + w - disp(k1) - disp(-1k1 + q)) *
                                (mu - v + w - disp(k1) + disp(k2) - disp(-1k - k1 + k2 + q))
                            ) ^ -1
                        ) +
                        (
                            fermidist(-1mu + disp(k3), beta) *
                            (-1 + fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) *
                            (
                                fermidist(-1mu + disp(k2), beta) -
                                fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
                            )
                        ) *
                        (
                            (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) *
                            (
                                mu - vp + w - disp(k3) + disp(-1k1 + k3 + kp) -
                                disp(-1k1 + q)
                            ) *
                            (
                                v + vp - w - disp(k2) + disp(k3) - disp(-1k1 + k3 + kp) +
                                disp(-1k - k1 + k2 + q)
                            )
                        ) ^ -1 +
                        (
                            (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                            ) *
                            fermidist(mu + w - disp(-1k1 + q), beta) *
                            (
                                fermidist(-1mu + disp(k2), beta) -
                                fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
                            )
                        ) *
                        (
                            (2mu + w - disp(k1) - disp(-1k1 + q)) *
                            (
                                mu - vp + w - disp(k3) + disp(-1k1 + k3 + kp) -
                                disp(-1k1 + q)
                            ) *
                            (mu + v - disp(k2) - disp(-1k1 + q) + disp(-1k - k1 + k2 + q))
                        ) ^ -1 - (
                            (
                                fermidist(-1mu + disp(k2), beta) *
                                (
                                    fermidist(-1mu + disp(k3), beta) -
                                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                                ) *
                                (-1 + fermidist(-1mu + disp(-1k - k1 + k2 + q), beta))
                            ) *
                            (
                                (
                                    mu - v + w - disp(k1) + disp(k2) -
                                    disp(-1k - k1 + k2 + q)
                                ) *
                                (
                                    v + vp - w - disp(k2) + disp(k3) -
                                    disp(-1k1 + k3 + kp) + disp(-1k - k1 + k2 + q)
                                ) *
                                (
                                    mu + v - disp(k2) - disp(-1k1 + q) +
                                    disp(-1k - k1 + k2 + q)
                                )
                            ) ^ -1
                        )
                    )
                ) * 2 ^ -1,
            ),
        ),
    ) +
    ifelse(
        and(
            isless(abs2(2mu + w - disp(k1) - disp(-1k1 + q)), atol),
            isless(abs2(2mu + w - disp(k3) - disp(-1k3 + q)), atol),
        ),
        2 *
        cnorm ^ 3 *
        u ^ 3 *
        dfermidist(-1mu + disp(k1), beta) *
        dfermidist(-1mu + disp(k3), beta),
        ifelse(
            isless(abs2(2mu + w - disp(k3) - disp(-1k3 + q)), atol),
            (
                2 *
                cnorm ^ 3 *
                u ^ 3 *
                dfermidist(-1mu + disp(k3), beta) *
                (fermidist(-1mu + disp(k1), beta) - fermidist(mu - disp(-1k1 + q), beta))
            ) * (-2mu + disp(k1) + disp(-1k1 + q)) ^ -1,
            ifelse(
                isless(abs2(2mu + w - disp(k1) - disp(-1k1 + q)), atol),
                (
                    2 *
                    cnorm ^ 3 *
                    u ^ 3 *
                    dfermidist(-1mu + disp(k1), beta) *
                    (
                        fermidist(-1mu + disp(k3), beta) -
                        fermidist(mu - disp(-1k3 + q), beta)
                    )
                ) * (-2mu + disp(k3) + disp(-1k3 + q)) ^ -1,
                (
                    2 *
                    cnorm ^ 3 *
                    u ^ 3 *
                    (
                        fermidist(-1mu + disp(k1), beta) -
                        fermidist(mu + w - disp(-1k1 + q), beta)
                    ) *
                    (
                        fermidist(-1mu + disp(k3), beta) -
                        fermidist(mu - disp(-1k3 + q), beta)
                    )
                ) *
                (
                    (2mu + w - disp(k1) - disp(-1k1 + q)) *
                    (2mu + w - disp(k3) - disp(-1k3 + q))
                ) ^ -1,
            ),
        ),
    ) +
    ifelse(
        and(
            isless(abs2(2mu + w - disp(k1) - disp(-1k1 + q)), atol),
            isless(abs2(2mu + w - disp(k3) - disp(-1k3 + q)), atol),
        ),
        (
            cnorm ^ 3 *
            u ^ 4 *
            dfermidist(-1mu + disp(k3), beta) *
            (
                fermidist(-1mu + disp(k2), beta) -
                fermidist(-1mu + disp(-1k + k1 + k2), beta)
            ) *
            (
                dfermidist(-1mu + disp(k1), beta) *
                (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) +
                fermidist(-1mu + disp(k1), beta) -
                fermidist(v - disp(k2) + disp(-1k + k1 + k2), beta)
            )
        ) * ((mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) ^ 2) ^ -1,
        ifelse(
            isless(abs2(2mu + w - disp(k3) - disp(-1k3 + q)), atol),
            -(
                (
                    cnorm ^ 3 *
                    u ^ 4 *
                    dfermidist(-1mu + disp(k3), beta) *
                    (
                        -(
                            (mu - v + disp(k2) - disp(-1k + k1 + k2) - disp(-1k1 + q)) *
                            fermidist(-1mu + disp(k1), beta) *
                            (
                                fermidist(-1mu + disp(k2), beta) -
                                fermidist(-1mu + disp(-1k + k1 + k2), beta)
                            )
                        ) +
                        (2mu - disp(k1) - disp(-1k1 + q)) *
                        (-1 + fermidist(-1mu + disp(k2), beta)) *
                        fermidist(-1mu + disp(-1k + k1 + k2), beta) - (
                            (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) *
                            (
                                fermidist(-1mu + disp(k2), beta) -
                                fermidist(-1mu + disp(-1k + k1 + k2), beta)
                            ) *
                            fermidist(mu - disp(-1k1 + q), beta)
                        )
                    )
                ) *
                (
                    (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) *
                    (mu - v + disp(k2) - disp(-1k + k1 + k2) - disp(-1k1 + q)) *
                    (-2mu + disp(k1) + disp(-1k1 + q))
                ) ^ -1
            ),
            ifelse(
                isless(abs2(2mu + w - disp(k1) - disp(-1k1 + q)), atol),
                (
                    cnorm ^ 3 *
                    u ^ 4 *
                    (
                        fermidist(-1mu + disp(k2), beta) -
                        fermidist(-1mu + disp(-1k + k1 + k2), beta)
                    ) *
                    (
                        dfermidist(-1mu + disp(k1), beta) *
                        (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) +
                        fermidist(-1mu + disp(k1), beta) -
                        fermidist(v - disp(k2) + disp(-1k + k1 + k2), beta)
                    ) *
                    (
                        fermidist(-1mu + disp(k3), beta) -
                        fermidist(mu - disp(-1k3 + q), beta)
                    )
                ) *
                (
                    (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) ^ 2 *
                    (-2mu + disp(k3) + disp(-1k3 + q))
                ) ^ -1,
                (
                    cnorm ^ 3 *
                    u ^ 4 *
                    (
                        (mu - v + w + disp(k2) - disp(-1k + k1 + k2) - disp(-1k1 + q)) *
                        fermidist(-1mu + disp(k1), beta) *
                        (
                            fermidist(-1mu + disp(k2), beta) -
                            fermidist(-1mu + disp(-1k + k1 + k2), beta)
                        ) - (
                            (2mu + w - disp(k1) - disp(-1k1 + q)) *
                            (-1 + fermidist(-1mu + disp(k2), beta)) *
                            fermidist(-1mu + disp(-1k + k1 + k2), beta)
                        ) +
                        (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) *
                        (
                            fermidist(-1mu + disp(k2), beta) -
                            fermidist(-1mu + disp(-1k + k1 + k2), beta)
                        ) *
                        fermidist(mu + w - disp(-1k1 + q), beta)
                    ) *
                    (
                        fermidist(-1mu + disp(k3), beta) -
                        fermidist(mu - disp(-1k3 + q), beta)
                    )
                ) *
                (
                    (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) *
                    (2mu + w - disp(k1) - disp(-1k1 + q)) *
                    (mu - v + w + disp(k2) - disp(-1k + k1 + k2) - disp(-1k1 + q)) *
                    (2mu + w - disp(k3) - disp(-1k3 + q))
                ) ^ -1,
            ),
        ),
    ) +
    ifelse(
        and(
            isless(abs2(2mu + w - disp(k1) - disp(-1k1 + q)), atol),
            isless(abs2(2mu + w - disp(k3) - disp(-1k3 + q)), atol),
        ),
        -(
            (
                cnorm ^ 3 *
                u ^ 4 *
                dfermidist(-1mu + disp(k3), beta) *
                (
                    dfermidist(-1mu + disp(k1), beta) *
                    (mu - v - disp(k1) + disp(k2) - disp(-1k - k1 + k2 + q)) +
                    fermidist(-1mu + disp(k1), beta) -
                    fermidist(-1v + disp(k2) - disp(-1k - k1 + k2 + q), beta)
                ) *
                (
                    fermidist(-1mu + disp(k2), beta) -
                    fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
                )
            ) * ((-1mu + v + disp(k1) - disp(k2) + disp(-1k - k1 + k2 + q)) ^ 2) ^ -1
        ),
        ifelse(
            isless(abs2(2mu + w - disp(k3) - disp(-1k3 + q)), atol),
            -(
                (
                    cnorm ^ 3 *
                    u ^ 4 *
                    dfermidist(-1mu + disp(k3), beta) *
                    (
                        (mu + v - disp(k2) - disp(-1k1 + q) + disp(-1k - k1 + k2 + q)) *
                        fermidist(-1mu + disp(k1), beta) *
                        (
                            fermidist(-1mu + disp(k2), beta) -
                            fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
                        ) +
                        (mu - v - disp(k1) + disp(k2) - disp(-1k - k1 + k2 + q)) *
                        fermidist(mu - disp(-1k1 + q), beta) *
                        (
                            fermidist(-1mu + disp(k2), beta) -
                            fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
                        ) +
                        (2mu - disp(k1) - disp(-1k1 + q)) *
                        fermidist(-1mu + disp(k2), beta) *
                        (-1 + fermidist(-1mu + disp(-1k - k1 + k2 + q), beta))
                    )
                ) *
                (
                    (-2mu + disp(k1) + disp(-1k1 + q)) *
                    (mu - v - disp(k1) + disp(k2) - disp(-1k - k1 + k2 + q)) *
                    (mu + v - disp(k2) - disp(-1k1 + q) + disp(-1k - k1 + k2 + q))
                ) ^ -1
            ),
            ifelse(
                isless(abs2(2mu + w - disp(k1) - disp(-1k1 + q)), atol),
                -(
                    (
                        cnorm ^ 3 *
                        u ^ 4 *
                        (
                            dfermidist(-1mu + disp(k1), beta) *
                            (mu - v - disp(k1) + disp(k2) - disp(-1k - k1 + k2 + q)) +
                            fermidist(-1mu + disp(k1), beta) -
                            fermidist(-1v + disp(k2) - disp(-1k - k1 + k2 + q), beta)
                        ) *
                        (
                            fermidist(-1mu + disp(k2), beta) -
                            fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
                        ) *
                        (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(mu - disp(-1k3 + q), beta)
                        )
                    ) *
                    (
                        (-1mu + v + disp(k1) - disp(k2) + disp(-1k - k1 + k2 + q)) ^ 2 *
                        (-2mu + disp(k3) + disp(-1k3 + q))
                    ) ^ -1
                ),
                -(
                    (
                        cnorm ^ 3 *
                        u ^ 4 *
                        (
                            (mu + v - disp(k2) - disp(-1k1 + q) + disp(-1k - k1 + k2 + q)) *
                            fermidist(-1mu + disp(k1), beta) *
                            (
                                fermidist(-1mu + disp(k2), beta) -
                                fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
                            ) +
                            (mu - v + w - disp(k1) + disp(k2) - disp(-1k - k1 + k2 + q)) *
                            fermidist(mu + w - disp(-1k1 + q), beta) *
                            (
                                fermidist(-1mu + disp(k2), beta) -
                                fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
                            ) +
                            (2mu + w - disp(k1) - disp(-1k1 + q)) *
                            fermidist(-1mu + disp(k2), beta) *
                            (-1 + fermidist(-1mu + disp(-1k - k1 + k2 + q), beta))
                        ) *
                        (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(mu - disp(-1k3 + q), beta)
                        )
                    ) *
                    (
                        (2mu + w - disp(k1) - disp(-1k1 + q)) *
                        (mu - v + w - disp(k1) + disp(k2) - disp(-1k - k1 + k2 + q)) *
                        (mu + v - disp(k2) - disp(-1k1 + q) + disp(-1k - k1 + k2 + q)) *
                        (2mu + w - disp(k3) - disp(-1k3 + q))
                    ) ^ -1
                ),
            ),
        ),
    ) +
    ifelse(
        and(
            isless(abs2(2mu + w - disp(k1) - disp(-1k1 + q)), atol),
            isless(
                abs2(
                    v - vp - disp(k2) + disp(k3) + disp(-1k - k1 + k2 + q) -
                    disp(-1k1 + k3 - kp + q),
                ),
                atol,
            ),
        ),
        (
            cnorm ^ 3 *
            u ^ 4 *
            (
                dfermidist(-1mu + disp(k2) - disp(k3) + disp(-1k1 + k3 - kp + q), beta) *
                (mu - v - disp(k1) + disp(k3) - disp(-1k1 + k3 - kp + q)) *
                fermidist(-1mu + disp(k3), beta) *
                (-1 + fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) - (
                    (
                        fermidist(-1mu + disp(k3), beta) -
                        fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                    ) * (
                        -(
                            (mu - v - disp(k1) + disp(k3) - disp(-1k1 + k3 - kp + q)) * (
                                dfermidist(
                                    -1mu + disp(k2) - disp(k3) + disp(-1k1 + k3 - kp + q),
                                    beta,
                                ) *
                                fermidist(-1v + disp(k3) - disp(-1k1 + k3 - kp + q), beta) +
                                dfermidist(
                                    -1v + disp(k3) - disp(-1k1 + k3 - kp + q),
                                    beta,
                                ) * (
                                    fermidist(-1mu + disp(k2), beta) - fermidist(
                                        -1mu + disp(k2) - disp(k3) +
                                        disp(-1k1 + k3 - kp + q),
                                        beta,
                                    )
                                )
                            )
                        ) - (
                            2 *
                            fermidist(-1mu + disp(k2), beta) *
                            (
                                -1 + fermidist(
                                    -1mu + disp(k2) - disp(k3) + disp(-1k1 + k3 - kp + q),
                                    beta,
                                )
                            )
                        )
                    )
                ) +
                (
                    dfermidist(-1mu + disp(k1), beta) *
                    (mu - v - disp(k1) + disp(k3) - disp(-1k1 + k3 - kp + q)) +
                    2 * fermidist(-1mu + disp(k1), beta)
                ) *
                (
                    fermidist(-1mu + disp(k3), beta) -
                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                ) *
                (
                    fermidist(-1mu + disp(k2), beta) -
                    fermidist(-1mu + disp(k2) - disp(k3) + disp(-1k1 + k3 - kp + q), beta)
                )
            )
        ) * (2 * (mu - v - disp(k1) + disp(k3) - disp(-1k1 + k3 - kp + q)) ^ 3) ^ -1,
        ifelse(
            isless(
                abs2(
                    v - vp - disp(k2) + disp(k3) + disp(-1k - k1 + k2 + q) -
                    disp(-1k1 + k3 - kp + q),
                ),
                atol,
            ),
            (
                cnorm ^ 3 *
                u ^ 4 *
                (
                    -(
                        (
                            dfermidist(
                                -1mu + disp(k2) - disp(k3) + disp(-1k1 + k3 - kp + q),
                                beta,
                            ) *
                            fermidist(-1mu + disp(k3), beta) *
                            (-1 + fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta))
                        ) *
                        (
                            (mu - v + w - disp(k1) + disp(k3) - disp(-1k1 + k3 - kp + q)) *
                            (mu + v - disp(k3) - disp(-1k1 + q) + disp(-1k1 + k3 - kp + q))
                        ) ^ -1
                    ) +
                    (
                        (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                        ) * (
                            (
                                mu + v - disp(k3) - disp(-1k1 + q) +
                                disp(-1k1 + k3 - kp + q)
                            ) * (
                                -(
                                    (
                                        mu - v + w - disp(k1) + disp(k3) -
                                        disp(-1k1 + k3 - kp + q)
                                    ) * (
                                        dfermidist(
                                            -1mu + disp(k2) - disp(k3) +
                                            disp(-1k1 + k3 - kp + q),
                                            beta,
                                        ) * fermidist(
                                            -1v + w + disp(k3) - disp(-1k1 + k3 - kp + q),
                                            beta,
                                        ) +
                                        dfermidist(
                                            -1v + w + disp(k3) - disp(-1k1 + k3 - kp + q),
                                            beta,
                                        ) * (
                                            fermidist(-1mu + disp(k2), beta) - fermidist(
                                                -1mu + disp(k2) - disp(k3) +
                                                disp(-1k1 + k3 - kp + q),
                                                beta,
                                            )
                                        )
                                    )
                                ) - (
                                    fermidist(-1mu + disp(k2), beta) * (
                                        -1 + fermidist(
                                            -1mu + disp(k2) - disp(k3) +
                                            disp(-1k1 + k3 - kp + q),
                                            beta,
                                        )
                                    )
                                )
                            ) +
                            (mu - v + w - disp(k1) + disp(k3) - disp(-1k1 + k3 - kp + q)) *
                            fermidist(-1mu + disp(k2), beta) *
                            (
                                -1 + fermidist(
                                    -1mu + disp(k2) - disp(k3) + disp(-1k1 + k3 - kp + q),
                                    beta,
                                )
                            )
                        )
                    ) *
                    (
                        (mu - v + w - disp(k1) + disp(k3) - disp(-1k1 + k3 - kp + q)) ^ 2 *
                        (mu + v - disp(k3) - disp(-1k1 + q) + disp(-1k1 + k3 - kp + q)) ^ 2
                    ) ^ -1 - (
                        (
                            fermidist(-1mu + disp(k1), beta) *
                            (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                            ) *
                            (
                                fermidist(-1mu + disp(k2), beta) - fermidist(
                                    -1mu + disp(k2) - disp(k3) + disp(-1k1 + k3 - kp + q),
                                    beta,
                                )
                            )
                        ) *
                        (
                            (2mu + w - disp(k1) - disp(-1k1 + q)) *
                            (mu - v + w - disp(k1) + disp(k3) - disp(-1k1 + k3 - kp + q)) ^
                            2
                        ) ^ -1
                    ) +
                    (
                        fermidist(mu + w - disp(-1k1 + q), beta) *
                        (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                        ) *
                        (
                            fermidist(-1mu + disp(k2), beta) - fermidist(
                                -1mu + disp(k2) - disp(k3) + disp(-1k1 + k3 - kp + q),
                                beta,
                            )
                        )
                    ) *
                    (
                        (2mu + w - disp(k1) - disp(-1k1 + q)) *
                        (mu + v - disp(k3) - disp(-1k1 + q) + disp(-1k1 + k3 - kp + q)) ^ 2
                    ) ^ -1
                )
            ) * 2 ^ -1,
            ifelse(
                isless(abs2(2mu + w - disp(k1) - disp(-1k1 + q)), atol),
                (
                    cnorm ^ 3 *
                    u ^ 4 *
                    (
                        (
                            (
                                (
                                    -1mu + vp + disp(k1) - disp(k3) +
                                    disp(-1k1 + k3 - kp + q)
                                ) * (
                                    dfermidist(-1mu + disp(k1), beta) * (
                                        -1mu + v + disp(k1) - disp(k2) +
                                        disp(-1k - k1 + k2 + q)
                                    ) - fermidist(-1mu + disp(k1), beta)
                                ) - (
                                    (
                                        -1mu + v + disp(k1) - disp(k2) +
                                        disp(-1k - k1 + k2 + q)
                                    ) * fermidist(-1mu + disp(k1), beta)
                                )
                            ) *
                            (
                                fermidist(-1mu + disp(k2), beta) -
                                fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
                            ) *
                            (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                            )
                        ) *
                        (
                            (-1mu + v + disp(k1) - disp(k2) + disp(-1k - k1 + k2 + q)) ^ 2 *
                            (-1mu + vp + disp(k1) - disp(k3) + disp(-1k1 + k3 - kp + q)) ^
                            2
                        ) ^ -1 +
                        (
                            fermidist(-1mu + disp(k2), beta) *
                            (-1 + fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) *
                            (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                            )
                        ) *
                        (
                            (-1mu + v + disp(k1) - disp(k2) + disp(-1k - k1 + k2 + q)) ^ 2 *
                            (
                                v - vp - disp(k2) + disp(k3) + disp(-1k - k1 + k2 + q) -
                                disp(-1k1 + k3 - kp + q)
                            )
                        ) ^ -1 +
                        (
                            fermidist(-1mu + disp(k3), beta) *
                            (
                                fermidist(-1mu + disp(k2), beta) -
                                fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
                            ) *
                            (-1 + fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta))
                        ) *
                        (
                            (-1mu + vp + disp(k1) - disp(k3) + disp(-1k1 + k3 - kp + q)) ^
                            2 * (
                                -1v + vp + disp(k2) - disp(k3) - disp(-1k - k1 + k2 + q) +
                                disp(-1k1 + k3 - kp + q)
                            )
                        ) ^ -1
                    )
                ) * 2 ^ -1,
                (
                    cnorm ^ 3 *
                    u ^ 4 *
                    (
                        -(
                            (
                                fermidist(-1mu + disp(k1), beta) *
                                (
                                    fermidist(-1mu + disp(k2), beta) -
                                    fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
                                ) *
                                (
                                    fermidist(-1mu + disp(k3), beta) -
                                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                                )
                            ) *
                            (
                                (2mu + w - disp(k1) - disp(-1k1 + q)) *
                                (
                                    mu - v + w - disp(k1) + disp(k2) -
                                    disp(-1k - k1 + k2 + q)
                                ) *
                                (
                                    mu - vp + w - disp(k1) + disp(k3) -
                                    disp(-1k1 + k3 - kp + q)
                                )
                            ) ^ -1
                        ) +
                        (
                            fermidist(mu + w - disp(-1k1 + q), beta) *
                            (
                                fermidist(-1mu + disp(k2), beta) -
                                fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
                            ) *
                            (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                            )
                        ) *
                        (
                            (2mu + w - disp(k1) - disp(-1k1 + q)) *
                            (mu + v - disp(k2) - disp(-1k1 + q) + disp(-1k - k1 + k2 + q)) *
                            (mu + vp - disp(k3) - disp(-1k1 + q) + disp(-1k1 + k3 - kp + q))
                        ) ^ -1 - (
                            (
                                fermidist(-1mu + disp(k2), beta) *
                                (-1 + fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) *
                                (
                                    fermidist(-1mu + disp(k3), beta) -
                                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                                )
                            ) *
                            (
                                (
                                    mu - v + w - disp(k1) + disp(k2) -
                                    disp(-1k - k1 + k2 + q)
                                ) *
                                (
                                    mu + v - disp(k2) - disp(-1k1 + q) +
                                    disp(-1k - k1 + k2 + q)
                                ) *
                                (
                                    v - vp - disp(k2) + disp(k3) + disp(-1k - k1 + k2 + q) -
                                    disp(-1k1 + k3 - kp + q)
                                )
                            ) ^ -1
                        ) - (
                            (
                                fermidist(-1mu + disp(k3), beta) *
                                (
                                    fermidist(-1mu + disp(k2), beta) -
                                    fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
                                ) *
                                (-1 + fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta))
                            ) *
                            (
                                (
                                    mu - vp + w - disp(k1) + disp(k3) -
                                    disp(-1k1 + k3 - kp + q)
                                ) *
                                (
                                    mu + vp - disp(k3) - disp(-1k1 + q) +
                                    disp(-1k1 + k3 - kp + q)
                                ) *
                                (
                                    -1v + vp + disp(k2) - disp(k3) -
                                    disp(-1k - k1 + k2 + q) + disp(-1k1 + k3 - kp + q)
                                )
                            ) ^ -1
                        )
                    )
                ) * 2 ^ -1,
            ),
        ),
    ) +
    ifelse(
        and(
            isless(abs2(2mu + w - disp(k1) - disp(-1k1 + q)), atol),
            isless(
                abs2(
                    v + vp - w - disp(k2) + disp(-1k + k1 + k2) - disp(k3) +
                    disp(-1k1 + k3 - kp + q),
                ),
                atol,
            ),
        ),
        (
            cnorm ^ 3 *
            u ^ 4 *
            (
                fermidist(-1mu + disp(k2), beta) -
                fermidist(-1mu + disp(-1k + k1 + k2), beta)
            ) *
            (
                -(
                    (v + vp) *
                    (
                        (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) *
                        fermidist(-1mu + disp(k1), beta) +
                        (mu - vp - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) * (
                            dfermidist(-1mu + disp(k1), beta) *
                            (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) +
                            fermidist(-1mu + disp(k1), beta)
                        )
                    ) *
                    (
                        fermidist(-1mu + disp(k3), beta) -
                        fermidist(-1mu + disp(k2) - disp(-1k + k1 + k2) + disp(k3), beta)
                    )
                ) +
                (-1mu + vp + disp(k1) + disp(k2) - disp(-1k + k1 + k2)) ^ 2 *
                fermidist(-1mu + disp(k3), beta) *
                (-1 + fermidist(-1mu + disp(k2) - disp(-1k + k1 + k2) + disp(k3), beta)) -
                (
                    (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) ^ 2 *
                    fermidist(-1mu + disp(k3), beta) *
                    (-1 + fermidist(-1mu + disp(k2) - disp(-1k + k1 + k2) + disp(k3), beta))
                )
            )
        ) *
        (
            2 *
            (v + vp) *
            (-1mu + vp + disp(k1) + disp(k2) - disp(-1k + k1 + k2)) ^ 2 *
            (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) ^ 2
        ) ^ -1,
        ifelse(
            isless(
                abs2(
                    v + vp - w - disp(k2) + disp(-1k + k1 + k2) - disp(k3) +
                    disp(-1k1 + k3 - kp + q),
                ),
                atol,
            ),
            (
                cnorm ^ 3 *
                u ^ 4 *
                (
                    -(
                        (
                            dfermidist(
                                -1mu + disp(k2) - disp(-1k + k1 + k2) + disp(k3),
                                beta,
                            ) *
                            (-1 + fermidist(-1mu + disp(k2), beta)) *
                            fermidist(-1mu + disp(-1k + k1 + k2), beta)
                        ) *
                        (
                            (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) *
                            (mu + vp + disp(k2) - disp(-1k + k1 + k2) - disp(-1k1 + q))
                        ) ^ -1
                    ) - (
                        (
                            (
                                fermidist(-1mu + disp(k2), beta) -
                                fermidist(-1mu + disp(-1k + k1 + k2), beta)
                            ) * (
                                (
                                    mu + vp + disp(k2) - disp(-1k + k1 + k2) -
                                    disp(-1k1 + q)
                                ) * (
                                    -(
                                        (
                                            mu + v - disp(k1) - disp(k2) +
                                            disp(-1k + k1 + k2)
                                        ) * (
                                            dfermidist(
                                                -1mu + disp(k2) - disp(-1k + k1 + k2) +
                                                disp(k3),
                                                beta,
                                            ) * fermidist(
                                                v - disp(k2) + disp(-1k + k1 + k2),
                                                beta,
                                            ) +
                                            dfermidist(
                                                v - disp(k2) + disp(-1k + k1 + k2),
                                                beta,
                                            ) * (
                                                fermidist(-1mu + disp(k3), beta) -
                                                fermidist(
                                                    -1mu + disp(k2) - disp(-1k + k1 + k2) +
                                                    disp(k3),
                                                    beta,
                                                )
                                            )
                                        )
                                    ) - (
                                        fermidist(-1mu + disp(k3), beta) * (
                                            -1 + fermidist(
                                                -1mu + disp(k2) - disp(-1k + k1 + k2) +
                                                disp(k3),
                                                beta,
                                            )
                                        )
                                    )
                                ) +
                                (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) *
                                fermidist(-1mu + disp(k3), beta) *
                                (
                                    -1 + fermidist(
                                        -1mu + disp(k2) - disp(-1k + k1 + k2) + disp(k3),
                                        beta,
                                    )
                                )
                            )
                        ) *
                        (
                            (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) ^ 2 *
                            (mu + vp + disp(k2) - disp(-1k + k1 + k2) - disp(-1k1 + q)) ^ 2
                        ) ^ -1
                    ) +
                    (
                        fermidist(-1mu + disp(k1), beta) *
                        (
                            fermidist(-1mu + disp(k2), beta) -
                            fermidist(-1mu + disp(-1k + k1 + k2), beta)
                        ) *
                        (
                            fermidist(-1mu + disp(k3), beta) - fermidist(
                                -1mu + disp(k2) - disp(-1k + k1 + k2) + disp(k3),
                                beta,
                            )
                        )
                    ) *
                    (
                        (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) ^ 2 *
                        (2mu + v + vp - disp(k1) - disp(-1k1 + q))
                    ) ^ -1 - (
                        (
                            (
                                fermidist(-1mu + disp(k2), beta) -
                                fermidist(-1mu + disp(-1k + k1 + k2), beta)
                            ) *
                            (
                                fermidist(-1mu + disp(k3), beta) - fermidist(
                                    -1mu + disp(k2) - disp(-1k + k1 + k2) + disp(k3),
                                    beta,
                                )
                            ) *
                            fermidist(mu + v + vp - disp(-1k1 + q), beta)
                        ) *
                        (
                            (2mu + v + vp - disp(k1) - disp(-1k1 + q)) *
                            (mu + vp + disp(k2) - disp(-1k + k1 + k2) - disp(-1k1 + q)) ^ 2
                        ) ^ -1
                    )
                )
            ) * 2 ^ -1,
            ifelse(
                isless(abs2(2mu + w - disp(k1) - disp(-1k1 + q)), atol),
                (-2 ^ -1) * (
                    (
                        cnorm ^ 3 *
                        u ^ 4 *
                        (
                            (
                                v + vp - disp(k2) + disp(-1k + k1 + k2) - disp(k3) +
                                disp(-1k1 + k3 - kp + q)
                            ) *
                            (
                                (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) *
                                fermidist(-1mu + disp(k1), beta) +
                                (mu - vp - disp(k1) + disp(k3) - disp(-1k1 + k3 - kp + q)) *
                                (
                                    dfermidist(-1mu + disp(k1), beta) *
                                    (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) +
                                    fermidist(-1mu + disp(k1), beta)
                                )
                            ) *
                            (
                                fermidist(-1mu + disp(k2), beta) -
                                fermidist(-1mu + disp(-1k + k1 + k2), beta)
                            ) *
                            (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                            ) +
                            (-1mu + vp + disp(k1) - disp(k3) + disp(-1k1 + k3 - kp + q)) ^
                            2 *
                            (-1 + fermidist(-1mu + disp(k2), beta)) *
                            fermidist(-1mu + disp(-1k + k1 + k2), beta) *
                            (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                            ) +
                            (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) ^ 2 *
                            (
                                fermidist(-1mu + disp(k2), beta) -
                                fermidist(-1mu + disp(-1k + k1 + k2), beta)
                            ) *
                            fermidist(-1mu + disp(k3), beta) *
                            (-1 + fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta))
                        )
                    ) *
                    (
                        (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) ^ 2 *
                        (-1mu + vp + disp(k1) - disp(k3) + disp(-1k1 + k3 - kp + q)) ^ 2 *
                        (
                            v + vp - disp(k2) + disp(-1k + k1 + k2) - disp(k3) +
                            disp(-1k1 + k3 - kp + q)
                        )
                    ) ^ -1
                ),
                (
                    cnorm ^ 3 *
                    u ^ 4 *
                    (
                        (
                            fermidist(-1mu + disp(k1), beta) *
                            (
                                fermidist(-1mu + disp(k2), beta) -
                                fermidist(-1mu + disp(-1k + k1 + k2), beta)
                            ) *
                            (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                            )
                        ) *
                        (
                            (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) *
                            (2mu + w - disp(k1) - disp(-1k1 + q)) *
                            (mu - vp + w - disp(k1) + disp(k3) - disp(-1k1 + k3 - kp + q))
                        ) ^ -1 +
                        (
                            (-1 + fermidist(-1mu + disp(k2), beta)) *
                            fermidist(-1mu + disp(-1k + k1 + k2), beta) *
                            (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                            )
                        ) *
                        (
                            (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) *
                            (mu - v + w + disp(k2) - disp(-1k + k1 + k2) - disp(-1k1 + q)) *
                            (
                                v + vp - w - disp(k2) + disp(-1k + k1 + k2) - disp(k3) +
                                disp(-1k1 + k3 - kp + q)
                            )
                        ) ^ -1 - (
                            (
                                (
                                    fermidist(-1mu + disp(k2), beta) -
                                    fermidist(-1mu + disp(-1k + k1 + k2), beta)
                                ) *
                                fermidist(mu + w - disp(-1k1 + q), beta) *
                                (
                                    fermidist(-1mu + disp(k3), beta) -
                                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                                )
                            ) *
                            (
                                (2mu + w - disp(k1) - disp(-1k1 + q)) *
                                (
                                    mu - v + w + disp(k2) - disp(-1k + k1 + k2) -
                                    disp(-1k1 + q)
                                ) *
                                (
                                    mu + vp - disp(k3) - disp(-1k1 + q) +
                                    disp(-1k1 + k3 - kp + q)
                                )
                            ) ^ -1
                        ) +
                        (
                            (
                                fermidist(-1mu + disp(k2), beta) -
                                fermidist(-1mu + disp(-1k + k1 + k2), beta)
                            ) *
                            fermidist(-1mu + disp(k3), beta) *
                            (-1 + fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta))
                        ) *
                        (
                            (mu - vp + w - disp(k1) + disp(k3) - disp(-1k1 + k3 - kp + q)) *
                            (
                                v + vp - w - disp(k2) + disp(-1k + k1 + k2) - disp(k3) +
                                disp(-1k1 + k3 - kp + q)
                            ) *
                            (mu + vp - disp(k3) - disp(-1k1 + q) + disp(-1k1 + k3 - kp + q))
                        ) ^ -1
                    )
                ) * 2 ^ -1,
            ),
        ),
    ) +
    ifelse(
        isless(abs2(2mu + w - disp(k1) - disp(-1k1 + q)), atol),
        2 * cnorm ^ 3 * u ^ 2 * dfermidist(-1mu + disp(k1), beta),
        (
            -2 *
            cnorm ^ 3 *
            u ^ 2 *
            (fermidist(-1mu + disp(k1), beta) - fermidist(mu + w - disp(-1k1 + q), beta))
        ) * (2mu + w - disp(k1) - disp(-1k1 + q)) ^ -1,
    ) +
    ifelse(
        isless(abs2(2mu + w - disp(k1) - disp(-1k1 + q)), atol),
        (
            cnorm ^ 3 *
            u ^ 3 *
            (
                fermidist(-1mu + disp(k2), beta) -
                fermidist(-1mu + disp(-1k + k1 + k2), beta)
            ) *
            (
                dfermidist(-1mu + disp(k1), beta) *
                (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) +
                fermidist(-1mu + disp(k1), beta) -
                fermidist(v - disp(k2) + disp(-1k + k1 + k2), beta)
            )
        ) * ((mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) ^ 2) ^ -1,
        (
            cnorm ^ 3 *
            u ^ 3 *
            (
                -(
                    (mu - v + w + disp(k2) - disp(-1k + k1 + k2) - disp(-1k1 + q)) *
                    fermidist(-1mu + disp(k1), beta) *
                    (
                        fermidist(-1mu + disp(k2), beta) -
                        fermidist(-1mu + disp(-1k + k1 + k2), beta)
                    )
                ) +
                (2mu + w - disp(k1) - disp(-1k1 + q)) *
                (-1 + fermidist(-1mu + disp(k2), beta)) *
                fermidist(-1mu + disp(-1k + k1 + k2), beta) - (
                    (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) *
                    (
                        fermidist(-1mu + disp(k2), beta) -
                        fermidist(-1mu + disp(-1k + k1 + k2), beta)
                    ) *
                    fermidist(mu + w - disp(-1k1 + q), beta)
                )
            )
        ) *
        (
            (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) *
            (2mu + w - disp(k1) - disp(-1k1 + q)) *
            (mu - v + w + disp(k2) - disp(-1k + k1 + k2) - disp(-1k1 + q))
        ) ^ -1,
    ) +
    ifelse(
        isless(abs2(2mu + w - disp(k1) - disp(-1k1 + q)), atol),
        -(
            (
                cnorm ^ 3 *
                u ^ 3 *
                (
                    dfermidist(-1mu + disp(k1), beta) *
                    (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) +
                    fermidist(-1mu + disp(k1), beta) -
                    fermidist(vp + disp(k3) - disp(-1k1 + k3 + kp), beta)
                ) *
                (
                    fermidist(-1mu + disp(k3), beta) -
                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                )
            ) * ((mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) ^ 2) ^ -1
        ),
        (
            cnorm ^ 3 *
            u ^ 3 *
            (
                (mu - vp + w - disp(k3) + disp(-1k1 + k3 + kp) - disp(-1k1 + q)) *
                fermidist(-1mu + disp(k1), beta) *
                (
                    fermidist(-1mu + disp(k3), beta) -
                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                ) +
                (2mu + w - disp(k1) - disp(-1k1 + q)) *
                fermidist(-1mu + disp(k3), beta) *
                (-1 + fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
                (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) *
                (
                    fermidist(-1mu + disp(k3), beta) -
                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                ) *
                fermidist(mu + w - disp(-1k1 + q), beta)
            )
        ) *
        (
            (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) *
            (2mu + w - disp(k1) - disp(-1k1 + q)) *
            (mu - vp + w - disp(k3) + disp(-1k1 + k3 + kp) - disp(-1k1 + q))
        ) ^ -1,
    ) +
    ifelse(
        isless(abs2(2mu + w - disp(k1) - disp(-1k1 + q)), atol),
        -(
            (
                cnorm ^ 3 *
                u ^ 3 *
                (
                    dfermidist(-1mu + disp(k1), beta) *
                    (mu - v - disp(k1) + disp(k2) - disp(-1k - k1 + k2 + q)) +
                    fermidist(-1mu + disp(k1), beta) -
                    fermidist(-1v + disp(k2) - disp(-1k - k1 + k2 + q), beta)
                ) *
                (
                    fermidist(-1mu + disp(k2), beta) -
                    fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
                )
            ) * ((-1mu + v + disp(k1) - disp(k2) + disp(-1k - k1 + k2 + q)) ^ 2) ^ -1
        ),
        (
            cnorm ^ 3 *
            u ^ 3 *
            (
                (mu + v - disp(k2) - disp(-1k1 + q) + disp(-1k - k1 + k2 + q)) *
                fermidist(-1mu + disp(k1), beta) *
                (
                    fermidist(-1mu + disp(k2), beta) -
                    fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
                ) +
                (mu - v + w - disp(k1) + disp(k2) - disp(-1k - k1 + k2 + q)) *
                fermidist(mu + w - disp(-1k1 + q), beta) *
                (
                    fermidist(-1mu + disp(k2), beta) -
                    fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
                ) +
                (2mu + w - disp(k1) - disp(-1k1 + q)) *
                fermidist(-1mu + disp(k2), beta) *
                (-1 + fermidist(-1mu + disp(-1k - k1 + k2 + q), beta))
            )
        ) *
        (
            (2mu + w - disp(k1) - disp(-1k1 + q)) *
            (mu - v + w - disp(k1) + disp(k2) - disp(-1k - k1 + k2 + q)) *
            (mu + v - disp(k2) - disp(-1k1 + q) + disp(-1k - k1 + k2 + q))
        ) ^ -1,
    ) +
    ifelse(
        isless(abs2(2mu + w - disp(k1) - disp(-1k1 + q)), atol),
        -(
            (
                cnorm ^ 3 *
                u ^ 3 *
                (
                    dfermidist(-1mu + disp(k1), beta) *
                    (mu - vp - disp(k1) + disp(k3) - disp(-1k1 + k3 - kp + q)) +
                    fermidist(-1mu + disp(k1), beta) -
                    fermidist(-1vp + disp(k3) - disp(-1k1 + k3 - kp + q), beta)
                ) *
                (
                    fermidist(-1mu + disp(k3), beta) -
                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                )
            ) * ((-1mu + vp + disp(k1) - disp(k3) + disp(-1k1 + k3 - kp + q)) ^ 2) ^ -1
        ),
        (
            cnorm ^ 3 *
            u ^ 3 *
            (
                (mu + vp - disp(k3) - disp(-1k1 + q) + disp(-1k1 + k3 - kp + q)) *
                fermidist(-1mu + disp(k1), beta) *
                (
                    fermidist(-1mu + disp(k3), beta) -
                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                ) +
                (mu - vp + w - disp(k1) + disp(k3) - disp(-1k1 + k3 - kp + q)) *
                fermidist(mu + w - disp(-1k1 + q), beta) *
                (
                    fermidist(-1mu + disp(k3), beta) -
                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                ) +
                (2mu + w - disp(k1) - disp(-1k1 + q)) *
                fermidist(-1mu + disp(k3), beta) *
                (-1 + fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta))
            )
        ) *
        (
            (2mu + w - disp(k1) - disp(-1k1 + q)) *
            (mu - vp + w - disp(k1) + disp(k3) - disp(-1k1 + k3 - kp + q)) *
            (mu + vp - disp(k3) - disp(-1k1 + q) + disp(-1k1 + k3 - kp + q))
        ) ^ -1,
    )
end
