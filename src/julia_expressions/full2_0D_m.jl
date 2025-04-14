function full2_0D_m(u, mu, beta, v, vp, w)

    -1u + (u ^ 2 * fermidist(-1mu, beta)) * (2mu + v + vp + w) ^ -1 -
    ((u ^ 2 * fermidist(mu + v + vp + w, beta)) * (2mu + v + vp + w) ^ -1) +
    ifelse(iszero(w), u ^ 2 * dfermidist(-1mu, beta), 0)
end
