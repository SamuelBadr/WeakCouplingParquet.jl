function gamma2_0D_s(u, mu, beta, v, vp, w)
    (
        u * (
            u * (2v - w) * fermidist(-1mu, beta) -
            (u * (v + vp - w) * fermidist(-1mu + v - vp, beta)) +
            (v - vp) * (2 * (v + vp - w) - (u * fermidist(-1mu + v + vp - w, beta)))
        )
    ) * ((v - vp) * (v + vp - w)) ^ -1
end
