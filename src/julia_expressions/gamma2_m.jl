function gamma2_m(u, mu, beta, v, vp, w, k::SVector, kp::SVector, q::SVector, k1::SVector)
    -1u +
    (
        (u ^ 2 * fermidist(-1mu + disp(k1), beta)) *
        (v - vp - disp(k1) + disp(-1k + k1 + kp)) ^ -1 - (
            (u ^ 2 * fermidist(-1mu + v - vp + disp(-1k + k1 + kp), beta)) *
            (v - vp - disp(k1) + disp(-1k + k1 + kp)) ^ -1
        )
    ) * 2 ^ -1 +
    (
        -(
            (u ^ 2 * fermidist(-1mu + disp(k1), beta)) *
            (v - vp - disp(k1) + disp(-1k + k1 + kp)) ^ -1
        ) +
        (u ^ 2 * fermidist(-1mu + v - vp + disp(-1k + k1 + kp), beta)) *
        (v - vp - disp(k1) + disp(-1k + k1 + kp)) ^ -1
    ) * 2 ^ -1 +
    (
        (2 * u ^ 2 * fermidist(-1mu + disp(k1), beta)) *
        (2mu + v + vp + w - disp(k1) - disp(k - k1 + kp + q)) ^ -1 - (
            (2 * u ^ 2 * fermidist(mu + v + vp + w - disp(k - k1 + kp + q), beta)) *
            (2mu + v + vp + w - disp(k1) - disp(k - k1 + kp + q)) ^ -1
        )
    ) * 2 ^ -1
end
