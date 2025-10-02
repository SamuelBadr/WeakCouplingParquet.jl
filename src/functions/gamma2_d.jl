@cse function gamma2_d(
    u,
    mu,
    beta,
    v,
    vp,
    w,
    k::SVector,
    kp::SVector,
    q::SVector,
    k1::SVector,
)
    atol = 1e-6 #eps()
    cnorm = 1 / (2pi)^length(k)
    cnorm * u - (
        2 * ifelse(
            isless(abs2(v - vp - disp(k1) + disp(-1k + k1 + kp)), atol),
            cnorm * u ^ 2 * dfermidist(-1mu + disp(k1), beta),
            (
                cnorm *
                u ^ 2 *
                (
                    fermidist(-1mu + disp(k1), beta) -
                    fermidist(-1mu + disp(-1k + k1 + kp), beta)
                )
            ) * (-1v + vp + disp(k1) - disp(-1k + k1 + kp)) ^ -1,
        )
    ) +
    ifelse(
        isless(abs2(2mu + v + vp + w - disp(k1) - disp(k - k1 + kp + q)), atol),
        2 * cnorm * u ^ 2 * dfermidist(-1mu + disp(k1), beta),
        (
            2 *
            cnorm *
            u ^ 2 *
            (
                -fermidist(-1mu + disp(k1), beta) +
                fermidist(mu - disp(k - k1 + kp + q), beta)
            )
        ) * (2mu + v + vp + w - disp(k1) - disp(k - k1 + kp + q)) ^ -1,
    ) * 2 ^ -1
end
