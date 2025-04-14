function phi2_0D_t(u, mu, beta, v, vp, w)

    (
        -(
            fermidist(mu, beta) *
            (
                ifelse(iszero(-1mu - vp), u ^ 2 * dfermidist(-1mu, beta), 0) +
                ifelse(iszero(-1mu + vp - w), -(u ^ 2 * dfermidist(-1mu, beta)), 0)
            ) *
            (
                ifelse(iszero(-1mu - v), u ^ 2 * dfermidist(-1mu, beta), 0) +
                ifelse(iszero(mu - v + w), -(u ^ 2 * dfermidist(-1mu, beta)), 0)
            )
        ) +
        fermidist(-1mu, beta) *
        (
            ifelse(iszero(-1mu - v), -(u ^ 2 * dfermidist(-1mu, beta)), 0) +
            ifelse(iszero(mu - v + w), u ^ 2 * dfermidist(-1mu, beta), 0)
        ) *
        (
            ifelse(iszero(mu + vp), -(u ^ 2 * dfermidist(-1mu, beta)), 0) +
            ifelse(iszero(mu - vp + w), u ^ 2 * dfermidist(-1mu, beta), 0)
        )
    ) * (2 * (2mu + w)) ^ -1
end
