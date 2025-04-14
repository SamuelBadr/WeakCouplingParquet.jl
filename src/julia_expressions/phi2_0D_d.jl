function phi2_0D_d(u, mu, beta, v, vp, w)
    (
        u ^ 3 *
        (bosedist(-2mu, beta) + fermidist(-1mu, beta)) *
        (fermidist(-1mu, beta) - fermidist(mu, beta)) *
        (
            2 * mu ^ 2 +
            2 * mu * v +
            v ^ 2 +
            2 * mu * vp +
            vp ^ 2 +
            2 * mu * w +
            v * w +
            vp * w - (u * (2mu + v + vp + w) * fermidist(-1mu, beta)) +
            u * (2mu + v + vp + w) * fermidist(mu, beta)
        )
    ) * ((mu + v) * (mu + vp) * (mu + v + w) * (mu + vp + w)) ^ -1
end
