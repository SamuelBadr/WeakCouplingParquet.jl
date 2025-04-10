function full2_s(u, mu, beta, v, vp, w, k::SVector, kp::SVector, q::SVector, k1::SVector)
    2u +
    (u ^ 2 * fermidist(-1mu + disp(k1), beta)) *
    (v - vp - disp(k1) + disp(-1k + k1 + kp)) ^ -1 - (
        (2 * u ^ 2 * fermidist(-1mu + disp(k1), beta)) *
        (2mu + w - disp(k1) - disp(-1k1 + q)) ^ -1
    ) +
    (u ^ 2 * fermidist(-1mu + disp(k1), beta)) *
    (v + vp - w - disp(k1) + disp(-1k + k1 - kp + q)) ^ -1 - (
        (u ^ 2 * fermidist(-1mu + v - vp + disp(-1k + k1 + kp), beta)) *
        (v - vp - disp(k1) + disp(-1k + k1 + kp)) ^ -1
    ) +
    (2 * u ^ 2 * fermidist(mu - disp(-1k1 + q), beta)) *
    (2mu + w - disp(k1) - disp(-1k1 + q)) ^ -1 - (
        (u ^ 2 * fermidist(-1mu + v + vp - w + disp(-1k + k1 - kp + q), beta)) *
        (v + vp - w - disp(k1) + disp(-1k + k1 - kp + q)) ^ -1
    )
end
