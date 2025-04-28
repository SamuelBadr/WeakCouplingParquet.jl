function full2_0D_d(u, mu, beta, v, vp, w)

    u +
    (u ^ 2 * (-fermidist(-1mu, beta) + fermidist(mu + v + vp + w, beta))) *
    (2mu + v + vp + w) ^ -1 - (
        2 * ifelse(
            iszero(-1v + vp),
            u ^ 2 * dfermidist(-1mu, beta),
            (u ^ 2 * (-fermidist(-1mu, beta) + fermidist(-1mu + v - vp, beta))) *
            (v - vp) ^ -1,
        )
    ) + ifelse(iszero(w), u ^ 2 * dfermidist(-1mu, beta), 0)
end
