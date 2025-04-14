function full2_0D_s(u, mu, beta, v, vp, w)
    u * (
        2 +
        u *
        ((v - vp) ^ -1 + (v + vp - w) ^ -1 - (2 * (2mu + w) ^ -1)) *
        fermidist(-1mu, beta) +
        (2 * u * fermidist(mu, beta)) * (2mu + w) ^ -1 -
        ((u * fermidist(-1mu + v - vp, beta)) * (v - vp) ^ -1) -
        ((u * fermidist(-1mu + v + vp - w, beta)) * (v + vp - w) ^ -1)
    )
end
