function phi2_0D_d(u, mu, beta, v, vp, w)

    -(
        (
            u ^ 2 *
            fermidist(-2mu, beta) *
            (fermidist(-1mu, beta) - fermidist(mu, beta)) *
            (
                u ^ 2 * fermidist(-1mu, beta) - (u ^ 2 * fermidist(mu, beta)) - (
                    (v - vp) * (
                        u - (
                            2 * ifelse(
                                iszero(-2mu - v - vp - w),
                                u ^ 2 * dfermidist(-1mu, beta),
                                0,
                            )
                        )
                    )
                )
            )
        ) *
        ((mu + vp) * (-1v + vp) * (mu + vp + w) * (-1 + 2 * fermidist(-2mu, beta))) ^ -1
    ) +
    (
        fermidist(-1mu, beta) *
        (
            -(u ^ 2 * fermidist(-1mu, beta)) +
            u ^ 2 * fermidist(mu, beta) +
            (mu + v + w) *
            (u - (2 * ifelse(iszero(-1mu - v), u ^ 2 * dfermidist(-1mu, beta), 0)))
        ) *
        (
            -(u ^ 2 * fermidist(-1mu, beta)) +
            u ^ 2 * fermidist(mu, beta) +
            (mu + vp + w) * (
                u - (2 * ifelse(iszero(mu + vp), u ^ 2 * dfermidist(-1mu, beta), 0)) +
                ifelse(iszero(w), u ^ 2 * dfermidist(-1mu, beta), 0)
            )
        )
    ) * (w * (mu + v + w) * (mu + vp + w)) ^ -1 - (
        (
            fermidist(-1mu, beta) *
            (
                -(u ^ 2 * fermidist(-1mu, beta)) +
                u ^ 2 * fermidist(mu, beta) +
                (mu + v) *
                (u - (2 * ifelse(iszero(-1mu - v - w), u ^ 2 * dfermidist(-1mu, beta), 0)))
            ) *
            (
                -(u ^ 2 * fermidist(-1mu, beta)) +
                u ^ 2 * fermidist(mu, beta) +
                (mu + vp) * (
                    u + ifelse(iszero(w), u ^ 2 * dfermidist(-1mu, beta), 0) -
                    (2 * ifelse(iszero(mu + vp + w), u ^ 2 * dfermidist(-1mu, beta), 0))
                )
            )
        ) * ((mu + v) * (mu + vp) * w) ^ -1
    ) - (
        (
            u ^ 2 *
            fermidist(-2mu, beta) *
            (fermidist(-1mu, beta) - fermidist(mu, beta)) *
            (
                u ^ 2 * fermidist(-1mu, beta) - (u ^ 2 * fermidist(mu, beta)) +
                (v - vp) * (
                    u + ifelse(iszero(w), u ^ 2 * dfermidist(-1mu, beta), 0) - (
                        2 *
                        ifelse(iszero(2mu + v + vp + w), u ^ 2 * dfermidist(-1mu, beta), 0)
                    )
                )
            )
        ) * ((mu + v) * (v - vp) * (mu + v + w) * (-1 + 2 * fermidist(-2mu, beta))) ^ -1
    )
end
