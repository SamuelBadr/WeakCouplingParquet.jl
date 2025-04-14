function gamma2_0D_t(u, mu, beta, v, vp, w)
    -(
        (
            u ^ 2 * (
                (-2vp + w) * fermidist(-1mu, beta) +
                (v + vp - w) * fermidist(-1mu + v - vp, beta) +
                (-1v + vp) * fermidist(-1mu + v + vp - w, beta)
            )
        ) * ((v - vp) * (v + vp - w)) ^ -1
    )
end
