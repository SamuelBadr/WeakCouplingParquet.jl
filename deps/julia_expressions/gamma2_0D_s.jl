function gamma2_0D_s(u, mu, beta, v, vp, w)

    2u +
    ifelse(
        iszero(-1v + vp),
        -(u ^ 2 * dfermidist(-1mu, beta)),
        (u ^ 2 * (fermidist(-1mu, beta) - fermidist(-1mu + v - vp, beta))) * (v - vp) ^ -1,
    ) +
    ifelse(
        iszero(-1v - vp + w),
        -(u ^ 2 * dfermidist(-1mu, beta)),
        (u ^ 2 * (fermidist(-1mu, beta) - fermidist(-1mu + v + vp - w, beta))) *
        (v + vp - w) ^ -1,
    )
end
