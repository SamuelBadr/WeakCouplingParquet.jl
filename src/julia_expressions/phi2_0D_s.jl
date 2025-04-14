function phi2_0D_s(u, mu, beta, v, vp, w)

    (
        (
            2u +
            ifelse(iszero(-1mu - v), -(u ^ 2 * dfermidist(-1mu, beta)), 0) +
            ifelse(iszero(mu - v + w), -(u ^ 2 * dfermidist(-1mu, beta)), 0)
        ) * (
            2 * u ^ 2 * fermidist(-1mu, beta) ^ 2 +
            fermidist(mu, beta) * (
                2 * u ^ 2 * fermidist(mu, beta) +
                (2mu + w) * (
                    2u +
                    ifelse(iszero(-1mu - vp), -(u ^ 2 * dfermidist(-1mu, beta)), 0) +
                    ifelse(iszero(-1mu + vp - w), -(u ^ 2 * dfermidist(-1mu, beta)), 0)
                )
            ) - (
                fermidist(-1mu, beta) * (
                    4 * u ^ 2 * fermidist(mu, beta) +
                    (2mu + w) * (
                        2u +
                        ifelse(iszero(mu + vp), -(u ^ 2 * dfermidist(-1mu, beta)), 0) +
                        ifelse(iszero(mu - vp + w), -(u ^ 2 * dfermidist(-1mu, beta)), 0)
                    )
                )
            )
        )
    ) * (2 * (2mu + w) ^ 2) ^ -1
end
