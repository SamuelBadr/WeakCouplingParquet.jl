@cse function phi2_m(
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
            isless(abs2(w + disp(k1) - disp(k1 + q)), atol),
            isless(abs2(w + disp(k3) - disp(k3 + q)), atol),
        ),
        -(
            cnorm ^ 3 *
            u ^ 3 *
            dfermidist(-1mu + disp(k1), beta) *
            dfermidist(-1mu + disp(k3), beta)
        ),
        ifelse(
            isless(abs2(w + disp(k3) - disp(k3 + q)), atol),
            -(
                (
                    cnorm ^ 3 *
                    u ^ 3 *
                    dfermidist(-1mu + disp(k3), beta) *
                    (
                        fermidist(-1mu + disp(k1), beta) -
                        fermidist(-1mu + disp(k1 + q), beta)
                    )
                ) * (disp(k1) - disp(k1 + q)) ^ -1
            ),
            ifelse(
                isless(abs2(w + disp(k1) - disp(k1 + q)), atol),
                -(
                    (
                        cnorm ^ 3 *
                        u ^ 3 *
                        dfermidist(-1mu + disp(k1), beta) *
                        (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(-1mu + disp(k3 + q), beta)
                        )
                    ) * (disp(k3) - disp(k3 + q)) ^ -1
                ),
                -(
                    (
                        cnorm ^ 3 *
                        u ^ 3 *
                        (
                            fermidist(-1mu + disp(k1), beta) -
                            fermidist(-1mu - w + disp(k1 + q), beta)
                        ) *
                        (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(-1mu + disp(k3 + q), beta)
                        )
                    ) *
                    ((w + disp(k1) - disp(k1 + q)) * (w + disp(k3) - disp(k3 + q))) ^ -1
                ),
            ),
        ),
    ) +
    ifelse(
        and(
            isless(abs2(w + disp(k1) - disp(k1 + q)), atol),
            isless(abs2(w + disp(k3) - disp(k3 + q)), atol),
        ),
        (
            cnorm ^ 3 *
            u ^ 4 *
            dfermidist(-1mu + disp(k3), beta) *
            (
                fermidist(-1mu + disp(k2), beta) -
                fermidist(mu - disp(k + k1 - k2 + q), beta)
            ) *
            (
                dfermidist(-1mu + disp(k1), beta) *
                (mu + v + disp(k1) - disp(k2) - disp(k + k1 - k2 + q)) -
                fermidist(-1mu + disp(k1), beta) +
                fermidist(-2mu - v + disp(k2) + disp(k + k1 - k2 + q), beta)
            )
        ) * ((mu + v + disp(k1) - disp(k2) - disp(k + k1 - k2 + q)) ^ 2) ^ -1,
        ifelse(
            isless(abs2(w + disp(k3) - disp(k3 + q)), atol),
            (
                cnorm ^ 3 *
                u ^ 4 *
                dfermidist(-1mu + disp(k3), beta) *
                (
                    (mu + v - disp(k2) + disp(k1 + q) - disp(k + k1 - k2 + q)) *
                    fermidist(-1mu + disp(k1), beta) *
                    (
                        fermidist(-1mu + disp(k2), beta) -
                        fermidist(mu - disp(k + k1 - k2 + q), beta)
                    ) - (
                        (mu + v + disp(k1) - disp(k2) - disp(k + k1 - k2 + q)) *
                        fermidist(-1mu + disp(k1 + q), beta) *
                        (
                            fermidist(-1mu + disp(k2), beta) -
                            fermidist(mu - disp(k + k1 - k2 + q), beta)
                        )
                    ) - (
                        (disp(k1) - disp(k1 + q)) *
                        fermidist(-1mu + disp(k2), beta) *
                        (-1 + fermidist(mu - disp(k + k1 - k2 + q), beta))
                    )
                )
            ) *
            (
                (disp(k1) - disp(k1 + q)) *
                (mu + v + disp(k1) - disp(k2) - disp(k + k1 - k2 + q)) *
                (mu + v - disp(k2) + disp(k1 + q) - disp(k + k1 - k2 + q))
            ) ^ -1,
            ifelse(
                isless(abs2(w + disp(k1) - disp(k1 + q)), atol),
                (
                    cnorm ^ 3 *
                    u ^ 4 *
                    (
                        fermidist(-1mu + disp(k2), beta) -
                        fermidist(mu - disp(k + k1 - k2 + q), beta)
                    ) *
                    (
                        dfermidist(-1mu + disp(k1), beta) *
                        (mu + v + disp(k1) - disp(k2) - disp(k + k1 - k2 + q)) -
                        fermidist(-1mu + disp(k1), beta) +
                        fermidist(-2mu - v + disp(k2) + disp(k + k1 - k2 + q), beta)
                    ) *
                    (
                        fermidist(-1mu + disp(k3), beta) -
                        fermidist(-1mu + disp(k3 + q), beta)
                    )
                ) *
                (
                    (mu + v + disp(k1) - disp(k2) - disp(k + k1 - k2 + q)) ^ 2 *
                    (disp(k3) - disp(k3 + q))
                ) ^ -1,
                (
                    cnorm ^ 3 *
                    u ^ 4 *
                    (
                        (mu + v - disp(k2) + disp(k1 + q) - disp(k + k1 - k2 + q)) *
                        fermidist(-1mu + disp(k1), beta) *
                        (
                            fermidist(-1mu + disp(k2), beta) -
                            fermidist(mu - disp(k + k1 - k2 + q), beta)
                        ) - (
                            (mu + v + w + disp(k1) - disp(k2) - disp(k + k1 - k2 + q)) *
                            fermidist(-1mu - w + disp(k1 + q), beta) *
                            (
                                fermidist(-1mu + disp(k2), beta) -
                                fermidist(mu - disp(k + k1 - k2 + q), beta)
                            )
                        ) - (
                            (w + disp(k1) - disp(k1 + q)) *
                            fermidist(-1mu + disp(k2), beta) *
                            (-1 + fermidist(mu - disp(k + k1 - k2 + q), beta))
                        )
                    ) *
                    (
                        fermidist(-1mu + disp(k3), beta) -
                        fermidist(-1mu + disp(k3 + q), beta)
                    )
                ) *
                (
                    (w + disp(k1) - disp(k1 + q)) *
                    (mu + v + w + disp(k1) - disp(k2) - disp(k + k1 - k2 + q)) *
                    (mu + v - disp(k2) + disp(k1 + q) - disp(k + k1 - k2 + q)) *
                    (w + disp(k3) - disp(k3 + q))
                ) ^ -1,
            ),
        ),
    ) +
    ifelse(
        and(
            isless(abs2(w + disp(k1) - disp(k1 + q)), atol),
            isless(
                abs2(
                    v - vp - disp(k2) + disp(k3) - disp(k + k1 - k2 + q) +
                    disp(k1 - k3 + kp + q),
                ),
                atol,
            ),
        ),
        (
            cnorm ^ 3 *
            u ^ 4 *
            (
                dfermidist(mu + disp(k2) - disp(k3) - disp(k1 - k3 + kp + q), beta) *
                (mu + v + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q)) *
                fermidist(-1mu + disp(k3), beta) *
                (-1 + fermidist(mu - disp(k1 - k3 + kp + q), beta)) - (
                    (
                        -(
                            dfermidist(-1mu + disp(k1), beta) *
                            (mu + v + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q))
                        ) + 2 * fermidist(-1mu + disp(k1), beta)
                    ) *
                    (
                        fermidist(-1mu + disp(k3), beta) -
                        fermidist(mu - disp(k1 - k3 + kp + q), beta)
                    ) *
                    (
                        fermidist(-1mu + disp(k2), beta) -
                        fermidist(mu + disp(k2) - disp(k3) - disp(k1 - k3 + kp + q), beta)
                    )
                ) +
                (
                    fermidist(-1mu + disp(k3), beta) -
                    fermidist(mu - disp(k1 - k3 + kp + q), beta)
                ) * (
                    -2 *
                    fermidist(-1mu + disp(k2), beta) *
                    (
                        -1 +
                        fermidist(mu + disp(k2) - disp(k3) - disp(k1 - k3 + kp + q), beta)
                    ) +
                    (mu + v + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q)) * (
                        dfermidist(-2mu - v + disp(k3) + disp(k1 - k3 + kp + q), beta) * (
                            fermidist(-1mu + disp(k2), beta) - fermidist(
                                mu + disp(k2) - disp(k3) - disp(k1 - k3 + kp + q),
                                beta,
                            )
                        ) +
                        dfermidist(
                            mu + disp(k2) - disp(k3) - disp(k1 - k3 + kp + q),
                            beta,
                        ) * fermidist(-2mu - v + disp(k3) + disp(k1 - k3 + kp + q), beta)
                    )
                )
            )
        ) * ((mu + v + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q)) ^ 3) ^ -1,
        ifelse(
            isless(
                abs2(
                    v - vp - disp(k2) + disp(k3) - disp(k + k1 - k2 + q) +
                    disp(k1 - k3 + kp + q),
                ),
                atol,
            ),
            cnorm ^ 3 *
            u ^ 4 *
            (
                (
                    dfermidist(mu + disp(k2) - disp(k3) - disp(k1 - k3 + kp + q), beta) *
                    fermidist(-1mu + disp(k3), beta) *
                    (-1 + fermidist(mu - disp(k1 - k3 + kp + q), beta))
                ) *
                (
                    (mu + v + w + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q)) *
                    (mu + v - disp(k3) + disp(k1 + q) - disp(k1 - k3 + kp + q))
                ) ^ -1 +
                (
                    fermidist(-1mu + disp(k1), beta) *
                    (
                        fermidist(-1mu + disp(k3), beta) -
                        fermidist(mu - disp(k1 - k3 + kp + q), beta)
                    ) *
                    (
                        fermidist(-1mu + disp(k2), beta) -
                        fermidist(mu + disp(k2) - disp(k3) - disp(k1 - k3 + kp + q), beta)
                    )
                ) *
                (
                    (w + disp(k1) - disp(k1 + q)) *
                    (mu + v + w + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q)) ^ 2
                ) ^ -1 - (
                    (
                        fermidist(-1mu - w + disp(k1 + q), beta) *
                        (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(mu - disp(k1 - k3 + kp + q), beta)
                        ) *
                        (
                            fermidist(-1mu + disp(k2), beta) - fermidist(
                                mu + disp(k2) - disp(k3) - disp(k1 - k3 + kp + q),
                                beta,
                            )
                        )
                    ) *
                    (
                        (w + disp(k1) - disp(k1 + q)) *
                        (mu + v - disp(k3) + disp(k1 + q) - disp(k1 - k3 + kp + q)) ^ 2
                    ) ^ -1
                ) +
                (
                    (
                        fermidist(-1mu + disp(k3), beta) -
                        fermidist(mu - disp(k1 - k3 + kp + q), beta)
                    ) * (
                        -(
                            (mu + v + w + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q)) *
                            fermidist(-1mu + disp(k2), beta) *
                            (
                                -1 + fermidist(
                                    mu + disp(k2) - disp(k3) - disp(k1 - k3 + kp + q),
                                    beta,
                                )
                            )
                        ) +
                        (mu + v - disp(k3) + disp(k1 + q) - disp(k1 - k3 + kp + q)) * (
                            -(
                                fermidist(-1mu + disp(k2), beta) * (
                                    -1 + fermidist(
                                        mu + disp(k2) - disp(k3) - disp(k1 - k3 + kp + q),
                                        beta,
                                    )
                                )
                            ) +
                            (mu + v + w + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q)) * (
                                dfermidist(
                                    -2mu - v - w + disp(k3) + disp(k1 - k3 + kp + q),
                                    beta,
                                ) * (
                                    fermidist(-1mu + disp(k2), beta) - fermidist(
                                        mu + disp(k2) - disp(k3) - disp(k1 - k3 + kp + q),
                                        beta,
                                    )
                                ) +
                                dfermidist(
                                    mu + disp(k2) - disp(k3) - disp(k1 - k3 + kp + q),
                                    beta,
                                ) * fermidist(
                                    -2mu - v - w + disp(k3) + disp(k1 - k3 + kp + q),
                                    beta,
                                )
                            )
                        )
                    )
                ) *
                (
                    (mu + v + w + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q)) ^ 2 *
                    (mu + v - disp(k3) + disp(k1 + q) - disp(k1 - k3 + kp + q)) ^ 2
                ) ^ -1
            ),
            ifelse(
                isless(abs2(w + disp(k1) - disp(k1 + q)), atol),
                (
                    cnorm ^ 3 *
                    u ^ 4 *
                    (
                        -(
                            (
                                v - vp - disp(k2) + disp(k3) - disp(k + k1 - k2 + q) +
                                disp(k1 - k3 + kp + q)
                            ) *
                            (
                                (mu + v + disp(k1) - disp(k2) - disp(k + k1 - k2 + q)) *
                                fermidist(-1mu + disp(k1), beta) +
                                (mu + vp + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q)) *
                                (
                                    -(
                                        dfermidist(-1mu + disp(k1), beta) * (
                                            mu + v + disp(k1) - disp(k2) -
                                            disp(k + k1 - k2 + q)
                                        )
                                    ) + fermidist(-1mu + disp(k1), beta)
                                )
                            ) *
                            (
                                fermidist(-1mu + disp(k2), beta) -
                                fermidist(mu - disp(k + k1 - k2 + q), beta)
                            ) *
                            (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(mu - disp(k1 - k3 + kp + q), beta)
                            )
                        ) +
                        (mu + vp + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q)) ^ 2 *
                        fermidist(-1mu + disp(k2), beta) *
                        (-1 + fermidist(mu - disp(k + k1 - k2 + q), beta)) *
                        (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(mu - disp(k1 - k3 + kp + q), beta)
                        ) - (
                            (mu + v + disp(k1) - disp(k2) - disp(k + k1 - k2 + q)) ^ 2 *
                            fermidist(-1mu + disp(k3), beta) *
                            (
                                fermidist(-1mu + disp(k2), beta) -
                                fermidist(mu - disp(k + k1 - k2 + q), beta)
                            ) *
                            (-1 + fermidist(mu - disp(k1 - k3 + kp + q), beta))
                        )
                    )
                ) *
                (
                    (mu + v + disp(k1) - disp(k2) - disp(k + k1 - k2 + q)) ^ 2 *
                    (mu + vp + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q)) ^ 2 *
                    (
                        v - vp - disp(k2) + disp(k3) - disp(k + k1 - k2 + q) +
                        disp(k1 - k3 + kp + q)
                    )
                ) ^ -1,
                cnorm ^ 3 *
                u ^ 4 *
                (
                    (
                        fermidist(-1mu + disp(k1), beta) *
                        (
                            fermidist(-1mu + disp(k2), beta) -
                            fermidist(mu - disp(k + k1 - k2 + q), beta)
                        ) *
                        (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(mu - disp(k1 - k3 + kp + q), beta)
                        )
                    ) *
                    (
                        (w + disp(k1) - disp(k1 + q)) *
                        (mu + v + w + disp(k1) - disp(k2) - disp(k + k1 - k2 + q)) *
                        (mu + vp + w + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q))
                    ) ^ -1 - (
                        (
                            fermidist(-1mu - w + disp(k1 + q), beta) *
                            (
                                fermidist(-1mu + disp(k2), beta) -
                                fermidist(mu - disp(k + k1 - k2 + q), beta)
                            ) *
                            (
                                fermidist(-1mu + disp(k3), beta) -
                                fermidist(mu - disp(k1 - k3 + kp + q), beta)
                            )
                        ) *
                        (
                            (w + disp(k1) - disp(k1 + q)) *
                            (mu + v - disp(k2) + disp(k1 + q) - disp(k + k1 - k2 + q)) *
                            (mu + vp - disp(k3) + disp(k1 + q) - disp(k1 - k3 + kp + q))
                        ) ^ -1
                    ) +
                    (
                        fermidist(-1mu + disp(k2), beta) *
                        (-1 + fermidist(mu - disp(k + k1 - k2 + q), beta)) *
                        (
                            fermidist(-1mu + disp(k3), beta) -
                            fermidist(mu - disp(k1 - k3 + kp + q), beta)
                        )
                    ) *
                    (
                        (mu + v + w + disp(k1) - disp(k2) - disp(k + k1 - k2 + q)) *
                        (mu + v - disp(k2) + disp(k1 + q) - disp(k + k1 - k2 + q)) *
                        (
                            v - vp - disp(k2) + disp(k3) - disp(k + k1 - k2 + q) +
                            disp(k1 - k3 + kp + q)
                        )
                    ) ^ -1 +
                    (
                        fermidist(-1mu + disp(k3), beta) *
                        (
                            fermidist(-1mu + disp(k2), beta) -
                            fermidist(mu - disp(k + k1 - k2 + q), beta)
                        ) *
                        (-1 + fermidist(mu - disp(k1 - k3 + kp + q), beta))
                    ) *
                    (
                        (mu + vp + w + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q)) *
                        (mu + vp - disp(k3) + disp(k1 + q) - disp(k1 - k3 + kp + q)) *
                        (
                            -1v + vp + disp(k2) - disp(k3) + disp(k + k1 - k2 + q) -
                            disp(k1 - k3 + kp + q)
                        )
                    ) ^ -1
                ),
            ),
        ),
    ) +
    ifelse(
        isless(abs2(w + disp(k1) - disp(k1 + q)), atol),
        cnorm ^ 3 * u ^ 2 * dfermidist(-1mu + disp(k1), beta),
        (
            cnorm ^ 3 *
            u ^ 2 *
            (fermidist(-1mu + disp(k1), beta) - fermidist(-1mu - w + disp(k1 + q), beta))
        ) * (w + disp(k1) - disp(k1 + q)) ^ -1,
    ) +
    ifelse(
        isless(abs2(w + disp(k1) - disp(k1 + q)), atol),
        (
            cnorm ^ 3 *
            u ^ 3 *
            (
                fermidist(-1mu + disp(k2), beta) -
                fermidist(mu - disp(k + k1 - k2 + q), beta)
            ) *
            (
                -(
                    dfermidist(-1mu + disp(k1), beta) *
                    (mu + v + disp(k1) - disp(k2) - disp(k + k1 - k2 + q))
                ) + fermidist(-1mu + disp(k1), beta) -
                fermidist(-2mu - v + disp(k2) + disp(k + k1 - k2 + q), beta)
            )
        ) * ((mu + v + disp(k1) - disp(k2) - disp(k + k1 - k2 + q)) ^ 2) ^ -1,
        (
            cnorm ^ 3 *
            u ^ 3 *
            (
                -(
                    (mu + v - disp(k2) + disp(k1 + q) - disp(k + k1 - k2 + q)) *
                    fermidist(-1mu + disp(k1), beta) *
                    (
                        fermidist(-1mu + disp(k2), beta) -
                        fermidist(mu - disp(k + k1 - k2 + q), beta)
                    )
                ) +
                (mu + v + w + disp(k1) - disp(k2) - disp(k + k1 - k2 + q)) *
                fermidist(-1mu - w + disp(k1 + q), beta) *
                (
                    fermidist(-1mu + disp(k2), beta) -
                    fermidist(mu - disp(k + k1 - k2 + q), beta)
                ) +
                (w + disp(k1) - disp(k1 + q)) *
                fermidist(-1mu + disp(k2), beta) *
                (-1 + fermidist(mu - disp(k + k1 - k2 + q), beta))
            )
        ) *
        (
            (w + disp(k1) - disp(k1 + q)) *
            (mu + v + w + disp(k1) - disp(k2) - disp(k + k1 - k2 + q)) *
            (mu + v - disp(k2) + disp(k1 + q) - disp(k + k1 - k2 + q))
        ) ^ -1,
    ) +
    ifelse(
        isless(abs2(w + disp(k1) - disp(k1 + q)), atol),
        (
            cnorm ^ 3 *
            u ^ 3 *
            (
                fermidist(-1mu + disp(k3), beta) -
                fermidist(mu - disp(k1 - k3 + kp + q), beta)
            ) *
            (
                -(
                    dfermidist(-1mu + disp(k1), beta) *
                    (mu + vp + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q))
                ) + fermidist(-1mu + disp(k1), beta) -
                fermidist(-2mu - vp + disp(k3) + disp(k1 - k3 + kp + q), beta)
            )
        ) * ((mu + vp + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q)) ^ 2) ^ -1,
        (
            cnorm ^ 3 *
            u ^ 3 *
            (
                -(
                    (mu + vp - disp(k3) + disp(k1 + q) - disp(k1 - k3 + kp + q)) *
                    fermidist(-1mu + disp(k1), beta) *
                    (
                        fermidist(-1mu + disp(k3), beta) -
                        fermidist(mu - disp(k1 - k3 + kp + q), beta)
                    )
                ) +
                (mu + vp + w + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q)) *
                fermidist(-1mu - w + disp(k1 + q), beta) *
                (
                    fermidist(-1mu + disp(k3), beta) -
                    fermidist(mu - disp(k1 - k3 + kp + q), beta)
                ) +
                (w + disp(k1) - disp(k1 + q)) *
                fermidist(-1mu + disp(k3), beta) *
                (-1 + fermidist(mu - disp(k1 - k3 + kp + q), beta))
            )
        ) *
        (
            (w + disp(k1) - disp(k1 + q)) *
            (mu + vp + w + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q)) *
            (mu + vp - disp(k3) + disp(k1 + q) - disp(k1 - k3 + kp + q))
        ) ^ -1,
    )
end
