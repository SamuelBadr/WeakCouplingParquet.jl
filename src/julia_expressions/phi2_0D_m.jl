function phi2_0D_m(u, mu, beta, v, vp, w)

    -(
        (
            u ^ 2 *
            (
                -fermidist(-1mu, beta) +
                fermidist(-2mu, beta) * (-1 + 2 * fermidist(-1mu, beta))
            ) *
            (fermidist(-1mu, beta) - fermidist(mu, beta)) *
            (
                -2 * mu ^ 2 * u - (2 * mu * u * v) - (u * v ^ 2) - (2 * mu * u * vp) -
                (u * vp ^ 2) - (2 * mu * u * w) - (u * v * w) - (u * vp * w) +
                u ^ 2 * (2mu + v + vp + w) * fermidist(-1mu, beta) -
                (u ^ 2 * (2mu + v + vp + w) * fermidist(mu, beta)) +
                mu ^ 2 * ifelse(iszero(w), u ^ 2 * dfermidist(-1mu, beta), 0) +
                2 * mu * vp * ifelse(iszero(w), u ^ 2 * dfermidist(-1mu, beta), 0) +
                vp ^ 2 * ifelse(iszero(w), u ^ 2 * dfermidist(-1mu, beta), 0) +
                mu * w * ifelse(iszero(w), u ^ 2 * dfermidist(-1mu, beta), 0) +
                vp * w * ifelse(iszero(w), u ^ 2 * dfermidist(-1mu, beta), 0)
            )
        ) *
        (
            (mu + v) *
            (mu + vp) *
            (mu + v + w) *
            (mu + vp + w) *
            (-1 + 2 * fermidist(-2mu, beta))
        ) ^ -1
    )
end
