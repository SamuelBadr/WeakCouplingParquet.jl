function phi2_0D_s(u, mu, beta, v, vp, w)
    (
        -2 *
        u ^ 2 *
        (fermidist(-1mu, beta) - fermidist(mu, beta)) *
        (2mu + w - (u * fermidist(-1mu, beta)) + u * fermidist(mu, beta))
    ) * ((2mu + w) ^ 2) ^ -1
end
