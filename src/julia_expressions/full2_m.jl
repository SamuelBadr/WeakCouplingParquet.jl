function full2_m(u, mu, beta, v, vp, w, k::SVector, kp::SVector, q::SVector, k1::SVector)
    cnorm = 1 / (2pi)^length(k)
    (
        u * (
            -2mu - v - vp - w +
            disp(k1) +
            disp(k - k1 + kp + q) +
            cnorm *
            u *
            (
                fermidist(-1mu + disp(k1), beta) -
                fermidist(mu + v + vp + w - disp(k - k1 + kp + q), beta)
            )
        )
    ) * (2mu + v + vp + w - disp(k1) - disp(k - k1 + kp + q)) ^ -1 +
    cnorm * ifelse(
        and(iszero(w), iszero(q)),
        u ^ 2 * dfermidist(-1mu + disp(k1), beta),
        (
            u ^ 2 *
            (fermidist(-1mu + disp(k1), beta) - fermidist(-1mu + disp(k1 + q), beta))
        ) * (w + disp(k1) - disp(k1 + q)) ^ -1,
    )
end
