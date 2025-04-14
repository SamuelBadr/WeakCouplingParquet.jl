function gamma2_0D_d(u, mu, beta, v, vp, w)
    (
        u * (
            u * (4mu + v + 3vp + 2w) * fermidist(-1mu, beta) -
            (2 * u * (2mu + v + vp + w) * fermidist(-1mu + v - vp, beta)) +
            (v - vp) * (2mu + v + vp + w + u * fermidist(mu + v + vp + w, beta))
        )
    ) * ((v - vp) * (2mu + v + vp + w)) ^ -1
end
