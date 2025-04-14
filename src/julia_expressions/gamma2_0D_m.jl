function gamma2_0D_m(u, mu, beta, v, vp, w)

    (
        u * (
            -2mu - v - vp - w + u * fermidist(-1mu, beta) -
            (u * fermidist(mu + v + vp + w, beta))
        )
    ) * (2mu + v + vp + w) ^ -1
end
