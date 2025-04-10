function phi2_s(
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
    k2::SVector,
    k3::SVector,
)
    (
        u ^ 2 *
        bosedist(disp(k3) - disp(-1k1 + k3 + kp), beta) *
        (
            4 * mu * u * vp * fermidist(-1mu + disp(k3), beta) -
            (2 * mu * u * w * fermidist(-1mu + disp(k3), beta)) +
            2 * u * vp * w * fermidist(-1mu + disp(k3), beta) -
            (u * w ^ 2 * fermidist(-1mu + disp(k3), beta)) -
            (2 * u * vp * disp(k3) * fermidist(-1mu + disp(k3), beta)) +
            u * w * disp(k3) * fermidist(-1mu + disp(k3), beta) -
            (2 * mu * u * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) -
            (u * w * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) +
            u * disp(k3) * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta) -
            (2 * u * vp * disp(-1k3 + q) * fermidist(-1mu + disp(k3), beta)) +
            u * w * disp(-1k3 + q) * fermidist(-1mu + disp(k3), beta) +
            u * disp(-1k1 + k3 + kp) * disp(-1k3 + q) * fermidist(-1mu + disp(k3), beta) +
            2 * mu * u * disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta) +
            u * w * disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta) -
            (u * disp(k3) * disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta)) -
            (
                u *
                disp(-1k3 + q) *
                disp(-1k1 + k3 - kp + q) *
                fermidist(-1mu + disp(k3), beta)
            ) - (4 * mu * u * vp * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
            2 * mu * u * w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
            (2 * u * vp * w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
            u * w ^ 2 * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            2 * u * vp * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
            (u * w * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
            2 *
            mu *
            u *
            disp(-1k1 + k3 + kp) *
            fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            u * w * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
            (
                u *
                disp(k3) *
                disp(-1k1 + k3 + kp) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
            ) + 2 * u * vp * disp(-1k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
            (u * w * disp(-1k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) - (
                u *
                disp(-1k1 + k3 + kp) *
                disp(-1k3 + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
            ) - (
                2 *
                mu *
                u *
                disp(-1k1 + k3 - kp + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
            ) - (
                u *
                w *
                disp(-1k1 + k3 - kp + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
            ) +
            u *
            disp(k3) *
            disp(-1k1 + k3 - kp + q) *
            fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            u *
            disp(-1k3 + q) *
            disp(-1k1 + k3 - kp + q) *
            fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
        ) *
        (
            -2 * v ^ 2 + 2 * vp ^ 2 + 2 * v * w - (2 * vp * w) + 4 * v * disp(k2) -
            (2 * w * disp(k2)) - (2 * disp(k2) ^ 2) - (2 * v * disp(-1k + k1 + k2)) -
            (2 * vp * disp(-1k + k1 + k2)) +
            2 * w * disp(-1k + k1 + k2) +
            2 * disp(k2) * disp(-1k + k1 + k2) +
            4 * vp * disp(k3) - (2 * w * disp(k3)) - (2 * disp(-1k + k1 + k2) * disp(k3)) +
            2 * disp(k3) ^ 2 - (4 * vp * disp(-1k1 + k3 + kp)) +
            2 * w * disp(-1k1 + k3 + kp) +
            2 * disp(-1k + k1 + k2) * disp(-1k1 + k3 + kp) -
            (4 * disp(k3) * disp(-1k1 + k3 + kp)) + 2 * disp(-1k1 + k3 + kp) ^ 2 -
            (2 * v * disp(-1k - k1 + k2 + q)) +
            2 * vp * disp(-1k - k1 + k2 + q) +
            2 * disp(k2) * disp(-1k - k1 + k2 + q) -
            (2 * disp(-1k + k1 + k2) * disp(-1k - k1 + k2 + q)) +
            2 * disp(k3) * disp(-1k - k1 + k2 + q) -
            (2 * disp(-1k1 + k3 + kp) * disp(-1k - k1 + k2 + q)) -
            (2 * u * v * fermidist(-1mu + disp(k2), beta)) +
            u * w * fermidist(-1mu + disp(k2), beta) +
            2 * u * disp(k2) * fermidist(-1mu + disp(k2), beta) -
            (u * disp(-1k + k1 + k2) * fermidist(-1mu + disp(k2), beta)) -
            (u * disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(k2), beta)) +
            u * v * fermidist(-1mu + disp(-1k + k1 + k2), beta) +
            u * vp * fermidist(-1mu + disp(-1k + k1 + k2), beta) -
            (u * w * fermidist(-1mu + disp(-1k + k1 + k2), beta)) -
            (u * disp(k2) * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
            u * disp(k3) * fermidist(-1mu + disp(-1k + k1 + k2), beta) -
            (u * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
            u * disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(-1k + k1 + k2), beta) +
            u * v * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta) -
            (u * vp * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) -
            (u * disp(k2) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) +
            u * disp(-1k + k1 + k2) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta) -
            (u * disp(k3) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) +
            u * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
        )
    ) *
    (
        2 *
        (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) *
        (-1v + vp + disp(k2) - disp(-1k + k1 + k2) + disp(k3) - disp(-1k1 + k3 + kp)) *
        (mu - vp + w - disp(k3) + disp(-1k1 + k3 + kp) - disp(-1k1 + q)) *
        (
            v + vp - w - disp(k2) + disp(k3) - disp(-1k1 + k3 + kp) +
            disp(-1k - k1 + k2 + q)
        ) *
        (2mu + w - disp(k3) - disp(-1k3 + q)) *
        (2vp - w - disp(-1k1 + k3 + kp) + disp(-1k1 + k3 - kp + q))
    ) ^ -1 - (
        (
            u ^ 2 *
            fermidist(-1mu + disp(k1), beta) *
            (
                -2 * mu ^ 2 + 2 * v ^ 2 - (2 * mu * w) - (2 * v * w) +
                4 * mu * disp(k1) +
                2 * w * disp(k1) - (2 * disp(k1) ^ 2) - (4 * v * disp(k2)) +
                2 * w * disp(k2) +
                2 * disp(k2) ^ 2 - (2 * mu * disp(-1k + k1 + k2)) +
                2 * v * disp(-1k + k1 + k2) - (2 * w * disp(-1k + k1 + k2)) +
                2 * disp(k1) * disp(-1k + k1 + k2) - (2 * disp(k2) * disp(-1k + k1 + k2)) +
                2 * mu * disp(-1k - k1 + k2 + q) +
                2 * v * disp(-1k - k1 + k2 + q) - (2 * disp(k1) * disp(-1k - k1 + k2 + q)) -
                (2 * disp(k2) * disp(-1k - k1 + k2 + q)) +
                2 * disp(-1k + k1 + k2) * disp(-1k - k1 + k2 + q) +
                2 * u * v * fermidist(-1mu + disp(k2), beta) -
                (u * w * fermidist(-1mu + disp(k2), beta)) -
                (2 * u * disp(k2) * fermidist(-1mu + disp(k2), beta)) +
                u * disp(-1k + k1 + k2) * fermidist(-1mu + disp(k2), beta) +
                u * disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(k2), beta) +
                mu * u * fermidist(-1mu + disp(-1k + k1 + k2), beta) -
                (u * v * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
                u * w * fermidist(-1mu + disp(-1k + k1 + k2), beta) -
                (u * disp(k1) * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
                u * disp(k2) * fermidist(-1mu + disp(-1k + k1 + k2), beta) - (
                    u *
                    disp(-1k - k1 + k2 + q) *
                    fermidist(-1mu + disp(-1k + k1 + k2), beta)
                ) - (mu * u * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) -
                (u * v * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) +
                u * disp(k1) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta) +
                u * disp(k2) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta) -
                (u * disp(-1k + k1 + k2) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta))
            ) *
            (
                -4 * mu ^ 3 + 4 * mu * vp ^ 2 - (6 * mu ^ 2 * w) - (4 * mu * vp * w) +
                2 * vp ^ 2 * w - (2 * mu * w ^ 2) - (2 * vp * w ^ 2) +
                8 * mu ^ 2 * disp(k1) +
                8 * mu * w * disp(k1) +
                2 * w ^ 2 * disp(k1) - (4 * mu * disp(k1) ^ 2) - (2 * w * disp(k1) ^ 2) -
                (6 * mu ^ 2 * disp(k3)) - (2 * vp ^ 2 * disp(k3)) -
                (6 * mu * w * disp(k3)) + 2 * vp * w * disp(k3) - (2 * w ^ 2 * disp(k3)) +
                4 * mu * disp(k1) * disp(k3) +
                2 * w * disp(k1) * disp(k3) +
                2 * disp(k1) ^ 2 * disp(k3) - (4 * disp(k1) * disp(k3) ^ 2) +
                2 * disp(k3) ^ 3 +
                4 * mu ^ 2 * disp(-1k1 + k3 + kp) - (4 * mu * vp * disp(-1k1 + k3 + kp)) +
                6 * mu * w * disp(-1k1 + k3 + kp) - (2 * vp * w * disp(-1k1 + k3 + kp)) +
                2 * w ^ 2 * disp(-1k1 + k3 + kp) -
                (4 * mu * disp(k1) * disp(-1k1 + k3 + kp)) -
                (2 * w * disp(k1) * disp(-1k1 + k3 + kp)) +
                2 * mu * disp(k3) * disp(-1k1 + k3 + kp) +
                2 * vp * disp(k3) * disp(-1k1 + k3 + kp) +
                2 * disp(k1) * disp(k3) * disp(-1k1 + k3 + kp) -
                (2 * disp(k3) ^ 2 * disp(-1k1 + k3 + kp)) + 2 * mu ^ 2 * disp(-1k3 + q) -
                (2 * vp ^ 2 * disp(-1k3 + q)) +
                2 * mu * w * disp(-1k3 + q) +
                2 * vp * w * disp(-1k3 + q) - (4 * mu * disp(k1) * disp(-1k3 + q)) -
                (2 * w * disp(k1) * disp(-1k3 + q)) +
                2 * disp(k1) ^ 2 * disp(-1k3 + q) +
                4 * mu * disp(k3) * disp(-1k3 + q) +
                2 * w * disp(k3) * disp(-1k3 + q) -
                (4 * disp(k1) * disp(k3) * disp(-1k3 + q)) +
                2 * disp(k3) ^ 2 * disp(-1k3 + q) -
                (2 * mu * disp(-1k1 + k3 + kp) * disp(-1k3 + q)) +
                2 * vp * disp(-1k1 + k3 + kp) * disp(-1k3 + q) -
                (2 * w * disp(-1k1 + k3 + kp) * disp(-1k3 + q)) +
                2 * disp(k1) * disp(-1k1 + k3 + kp) * disp(-1k3 + q) -
                (2 * disp(k3) * disp(-1k1 + k3 + kp) * disp(-1k3 + q)) +
                4 * mu ^ 2 * disp(-1k1 + k3 - kp + q) +
                4 * mu * vp * disp(-1k1 + k3 - kp + q) +
                2 * mu * w * disp(-1k1 + k3 - kp + q) +
                2 * vp * w * disp(-1k1 + k3 - kp + q) -
                (4 * mu * disp(k1) * disp(-1k1 + k3 - kp + q)) -
                (2 * w * disp(k1) * disp(-1k1 + k3 - kp + q)) +
                2 * mu * disp(k3) * disp(-1k1 + k3 - kp + q) -
                (2 * vp * disp(k3) * disp(-1k1 + k3 - kp + q)) +
                2 * w * disp(k3) * disp(-1k1 + k3 - kp + q) +
                2 * disp(k1) * disp(k3) * disp(-1k1 + k3 - kp + q) -
                (2 * disp(k3) ^ 2 * disp(-1k1 + k3 - kp + q)) -
                (4 * mu * disp(-1k1 + k3 + kp) * disp(-1k1 + k3 - kp + q)) -
                (2 * w * disp(-1k1 + k3 + kp) * disp(-1k1 + k3 - kp + q)) +
                2 * disp(k3) * disp(-1k1 + k3 + kp) * disp(-1k1 + k3 - kp + q) -
                (2 * mu * disp(-1k3 + q) * disp(-1k1 + k3 - kp + q)) -
                (2 * vp * disp(-1k3 + q) * disp(-1k1 + k3 - kp + q)) +
                2 * disp(k1) * disp(-1k3 + q) * disp(-1k1 + k3 - kp + q) -
                (2 * disp(k3) * disp(-1k3 + q) * disp(-1k1 + k3 - kp + q)) +
                2 * disp(-1k1 + k3 + kp) * disp(-1k3 + q) * disp(-1k1 + k3 - kp + q) +
                6 * mu ^ 2 * u * fermidist(-1mu + disp(k3), beta) -
                (2 * u * vp ^ 2 * fermidist(-1mu + disp(k3), beta)) +
                6 * mu * u * w * fermidist(-1mu + disp(k3), beta) +
                2 * u * vp * w * fermidist(-1mu + disp(k3), beta) +
                u * w ^ 2 * fermidist(-1mu + disp(k3), beta) -
                (8 * mu * u * disp(k1) * fermidist(-1mu + disp(k3), beta)) -
                (4 * u * w * disp(k1) * fermidist(-1mu + disp(k3), beta)) +
                2 * u * disp(k1) ^ 2 * fermidist(-1mu + disp(k3), beta) +
                6 * mu * u * disp(k3) * fermidist(-1mu + disp(k3), beta) +
                3 * u * w * disp(k3) * fermidist(-1mu + disp(k3), beta) -
                (2 * u * disp(k1) * disp(k3) * fermidist(-1mu + disp(k3), beta)) -
                (4 * mu * u * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) +
                2 * u * vp * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta) -
                (3 * u * w * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) +
                2 * u * disp(k1) * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta) -
                (u * disp(k3) * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) -
                (2 * mu * u * disp(-1k3 + q) * fermidist(-1mu + disp(k3), beta)) -
                (u * w * disp(-1k3 + q) * fermidist(-1mu + disp(k3), beta)) +
                2 * u * disp(k1) * disp(-1k3 + q) * fermidist(-1mu + disp(k3), beta) -
                (2 * u * disp(k3) * disp(-1k3 + q) * fermidist(-1mu + disp(k3), beta)) +
                u *
                disp(-1k1 + k3 + kp) *
                disp(-1k3 + q) *
                fermidist(-1mu + disp(k3), beta) -
                (4 * mu * u * disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta)) -
                (2 * u * vp * disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta)) -
                (u * w * disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta)) +
                2 *
                u *
                disp(k1) *
                disp(-1k1 + k3 - kp + q) *
                fermidist(-1mu + disp(k3), beta) - (
                    u *
                    disp(k3) *
                    disp(-1k1 + k3 - kp + q) *
                    fermidist(-1mu + disp(k3), beta)
                ) +
                2 *
                u *
                disp(-1k1 + k3 + kp) *
                disp(-1k1 + k3 - kp + q) *
                fermidist(-1mu + disp(k3), beta) +
                u *
                disp(-1k3 + q) *
                disp(-1k1 + k3 - kp + q) *
                fermidist(-1mu + disp(k3), beta) -
                (2 * mu ^ 2 * u * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
                2 * mu * u * vp * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
                (3 * mu * u * w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
                u * vp * w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
                (u * w ^ 2 * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
                2 * mu * u * disp(k1) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
                u * w * disp(k1) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
                (mu * u * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
                (u * vp * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
                (u * disp(k1) * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
                u * disp(k3) ^ 2 * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
                mu * u * disp(-1k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
                (u * vp * disp(-1k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
                u * w * disp(-1k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) - (
                    u *
                    disp(k1) *
                    disp(-1k3 + q) *
                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                ) +
                u *
                disp(k3) *
                disp(-1k3 + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
                2 *
                mu *
                u *
                disp(-1k1 + k3 - kp + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
                u *
                w *
                disp(-1k1 + k3 - kp + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta) - (
                    u *
                    disp(k3) *
                    disp(-1k1 + k3 - kp + q) *
                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                ) - (
                    u *
                    disp(-1k3 + q) *
                    disp(-1k1 + k3 - kp + q) *
                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                ) - (2 * mu ^ 2 * u * fermidist(mu - disp(-1k3 + q), beta)) +
                2 * u * vp ^ 2 * fermidist(mu - disp(-1k3 + q), beta) -
                (2 * mu * u * w * fermidist(mu - disp(-1k3 + q), beta)) -
                (2 * u * vp * w * fermidist(mu - disp(-1k3 + q), beta)) +
                4 * mu * u * disp(k1) * fermidist(mu - disp(-1k3 + q), beta) +
                2 * u * w * disp(k1) * fermidist(mu - disp(-1k3 + q), beta) -
                (2 * u * disp(k1) ^ 2 * fermidist(mu - disp(-1k3 + q), beta)) -
                (4 * mu * u * disp(k3) * fermidist(mu - disp(-1k3 + q), beta)) -
                (2 * u * w * disp(k3) * fermidist(mu - disp(-1k3 + q), beta)) +
                4 * u * disp(k1) * disp(k3) * fermidist(mu - disp(-1k3 + q), beta) -
                (2 * u * disp(k3) ^ 2 * fermidist(mu - disp(-1k3 + q), beta)) +
                2 * mu * u * disp(-1k1 + k3 + kp) * fermidist(mu - disp(-1k3 + q), beta) -
                (2 * u * vp * disp(-1k1 + k3 + kp) * fermidist(mu - disp(-1k3 + q), beta)) +
                2 * u * w * disp(-1k1 + k3 + kp) * fermidist(mu - disp(-1k3 + q), beta) -
                (
                    2 *
                    u *
                    disp(k1) *
                    disp(-1k1 + k3 + kp) *
                    fermidist(mu - disp(-1k3 + q), beta)
                ) +
                2 *
                u *
                disp(k3) *
                disp(-1k1 + k3 + kp) *
                fermidist(mu - disp(-1k3 + q), beta) +
                2 *
                mu *
                u *
                disp(-1k1 + k3 - kp + q) *
                fermidist(mu - disp(-1k3 + q), beta) +
                2 * u * vp * disp(-1k1 + k3 - kp + q) * fermidist(mu - disp(-1k3 + q), beta) -
                (
                    2 *
                    u *
                    disp(k1) *
                    disp(-1k1 + k3 - kp + q) *
                    fermidist(mu - disp(-1k3 + q), beta)
                ) +
                2 *
                u *
                disp(k3) *
                disp(-1k1 + k3 - kp + q) *
                fermidist(mu - disp(-1k3 + q), beta) - (
                    2 *
                    u *
                    disp(-1k1 + k3 + kp) *
                    disp(-1k1 + k3 - kp + q) *
                    fermidist(mu - disp(-1k3 + q), beta)
                ) - (2 * mu ^ 2 * u * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) -
                (2 * mu * u * vp * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) -
                (mu * u * w * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) -
                (u * vp * w * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
                2 * mu * u * disp(k1) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) +
                u * w * disp(k1) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
                (mu * u * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
                u * vp * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
                (u * w * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) - (
                    u *
                    disp(k1) *
                    disp(k3) *
                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                ) +
                u * disp(k3) ^ 2 * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) +
                2 *
                mu *
                u *
                disp(-1k1 + k3 + kp) *
                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) +
                u *
                w *
                disp(-1k1 + k3 + kp) *
                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) - (
                    u *
                    disp(k3) *
                    disp(-1k1 + k3 + kp) *
                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                ) +
                mu * u * disp(-1k3 + q) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) +
                u * vp * disp(-1k3 + q) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
                (
                    u *
                    disp(k1) *
                    disp(-1k3 + q) *
                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                ) +
                u *
                disp(k3) *
                disp(-1k3 + q) *
                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) - (
                    u *
                    disp(-1k1 + k3 + kp) *
                    disp(-1k3 + q) *
                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                )
            )
        ) *
        (
            2 *
            (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) *
            (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) *
            (2mu + w - disp(k1) - disp(-1k1 + q)) *
            (mu - v + w - disp(k1) + disp(k2) - disp(-1k - k1 + k2 + q)) *
            (2mu + w - disp(k3) - disp(-1k3 + q)) *
            (mu - vp + w - disp(k1) + disp(k3) - disp(-1k1 + k3 - kp + q))
        ) ^ -1
    ) - (
        (
            u ^ 2 *
            bosedist(-disp(k2) + disp(-1k + k1 + k2), beta) *
            (
                2 * u * v * fermidist(-1mu + disp(k2), beta) -
                (u * w * fermidist(-1mu + disp(k2), beta)) -
                (2 * u * disp(k2) * fermidist(-1mu + disp(k2), beta)) +
                u * disp(-1k + k1 + k2) * fermidist(-1mu + disp(k2), beta) +
                u * disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(k2), beta) -
                (2 * u * v * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
                u * w * fermidist(-1mu + disp(-1k + k1 + k2), beta) +
                2 * u * disp(k2) * fermidist(-1mu + disp(-1k + k1 + k2), beta) -
                (u * disp(-1k + k1 + k2) * fermidist(-1mu + disp(-1k + k1 + k2), beta)) -
                (u * disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(-1k + k1 + k2), beta))
            ) *
            (
                4 * mu * v ^ 2 - (4 * mu * vp ^ 2) - (4 * mu * v * w) +
                2 * v ^ 2 * w +
                4 * mu * vp * w - (2 * vp ^ 2 * w) - (2 * v * w ^ 2) + 2 * vp * w ^ 2 -
                (8 * mu * v * disp(k2)) + 4 * mu * w * disp(k2) - (4 * v * w * disp(k2)) +
                2 * w ^ 2 * disp(k2) +
                4 * mu * disp(k2) ^ 2 +
                2 * w * disp(k2) ^ 2 +
                8 * mu * v * disp(-1k + k1 + k2) - (4 * mu * w * disp(-1k + k1 + k2)) +
                4 * v * w * disp(-1k + k1 + k2) - (2 * w ^ 2 * disp(-1k + k1 + k2)) -
                (8 * mu * disp(k2) * disp(-1k + k1 + k2)) -
                (4 * w * disp(k2) * disp(-1k + k1 + k2)) +
                4 * mu * disp(-1k + k1 + k2) ^ 2 +
                2 * w * disp(-1k + k1 + k2) ^ 2 - (8 * mu * v * disp(k3)) -
                (2 * v ^ 2 * disp(k3)) +
                2 * vp ^ 2 * disp(k3) +
                4 * mu * w * disp(k3) - (2 * v * w * disp(k3)) - (2 * vp * w * disp(k3)) +
                2 * w ^ 2 * disp(k3) +
                8 * mu * disp(k2) * disp(k3) +
                4 * v * disp(k2) * disp(k3) +
                2 * w * disp(k2) * disp(k3) - (2 * disp(k2) ^ 2 * disp(k3)) -
                (8 * mu * disp(-1k + k1 + k2) * disp(k3)) -
                (4 * v * disp(-1k + k1 + k2) * disp(k3)) -
                (2 * w * disp(-1k + k1 + k2) * disp(k3)) +
                4 * disp(k2) * disp(-1k + k1 + k2) * disp(k3) -
                (2 * disp(-1k + k1 + k2) ^ 2 * disp(k3)) +
                4 * mu * disp(k3) ^ 2 +
                4 * v * disp(k3) ^ 2 - (4 * disp(k2) * disp(k3) ^ 2) +
                4 * disp(-1k + k1 + k2) * disp(k3) ^ 2 - (2 * disp(k3) ^ 3) +
                4 * mu * v * disp(-1k1 + k3 + kp) +
                4 * mu * vp * disp(-1k1 + k3 + kp) - (4 * mu * w * disp(-1k1 + k3 + kp)) +
                2 * v * w * disp(-1k1 + k3 + kp) +
                2 * vp * w * disp(-1k1 + k3 + kp) - (2 * w ^ 2 * disp(-1k1 + k3 + kp)) -
                (4 * mu * disp(k2) * disp(-1k1 + k3 + kp)) -
                (2 * w * disp(k2) * disp(-1k1 + k3 + kp)) +
                4 * mu * disp(-1k + k1 + k2) * disp(-1k1 + k3 + kp) +
                2 * w * disp(-1k + k1 + k2) * disp(-1k1 + k3 + kp) -
                (4 * mu * disp(k3) * disp(-1k1 + k3 + kp)) -
                (2 * v * disp(k3) * disp(-1k1 + k3 + kp)) -
                (2 * vp * disp(k3) * disp(-1k1 + k3 + kp)) +
                2 * disp(k2) * disp(k3) * disp(-1k1 + k3 + kp) -
                (2 * disp(-1k + k1 + k2) * disp(k3) * disp(-1k1 + k3 + kp)) +
                2 * disp(k3) ^ 2 * disp(-1k1 + k3 + kp) - (2 * v ^ 2 * disp(-1k3 + q)) +
                2 * vp ^ 2 * disp(-1k3 + q) +
                2 * v * w * disp(-1k3 + q) - (2 * vp * w * disp(-1k3 + q)) +
                4 * v * disp(k2) * disp(-1k3 + q) - (2 * w * disp(k2) * disp(-1k3 + q)) -
                (2 * disp(k2) ^ 2 * disp(-1k3 + q)) -
                (4 * v * disp(-1k + k1 + k2) * disp(-1k3 + q)) +
                2 * w * disp(-1k + k1 + k2) * disp(-1k3 + q) +
                4 * disp(k2) * disp(-1k + k1 + k2) * disp(-1k3 + q) -
                (2 * disp(-1k + k1 + k2) ^ 2 * disp(-1k3 + q)) +
                4 * v * disp(k3) * disp(-1k3 + q) - (2 * w * disp(k3) * disp(-1k3 + q)) -
                (4 * disp(k2) * disp(k3) * disp(-1k3 + q)) +
                4 * disp(-1k + k1 + k2) * disp(k3) * disp(-1k3 + q) -
                (2 * disp(k3) ^ 2 * disp(-1k3 + q)) -
                (2 * v * disp(-1k1 + k3 + kp) * disp(-1k3 + q)) -
                (2 * vp * disp(-1k1 + k3 + kp) * disp(-1k3 + q)) +
                2 * w * disp(-1k1 + k3 + kp) * disp(-1k3 + q) +
                2 * disp(k2) * disp(-1k1 + k3 + kp) * disp(-1k3 + q) -
                (2 * disp(-1k + k1 + k2) * disp(-1k1 + k3 + kp) * disp(-1k3 + q)) +
                2 * disp(k3) * disp(-1k1 + k3 + kp) * disp(-1k3 + q) +
                4 * mu * v * disp(-1k1 + k3 - kp + q) -
                (4 * mu * vp * disp(-1k1 + k3 - kp + q)) +
                2 * v * w * disp(-1k1 + k3 - kp + q) -
                (2 * vp * w * disp(-1k1 + k3 - kp + q)) -
                (4 * mu * disp(k2) * disp(-1k1 + k3 - kp + q)) -
                (2 * w * disp(k2) * disp(-1k1 + k3 - kp + q)) +
                4 * mu * disp(-1k + k1 + k2) * disp(-1k1 + k3 - kp + q) +
                2 * w * disp(-1k + k1 + k2) * disp(-1k1 + k3 - kp + q) -
                (4 * mu * disp(k3) * disp(-1k1 + k3 - kp + q)) -
                (2 * v * disp(k3) * disp(-1k1 + k3 - kp + q)) +
                2 * vp * disp(k3) * disp(-1k1 + k3 - kp + q) -
                (2 * w * disp(k3) * disp(-1k1 + k3 - kp + q)) +
                2 * disp(k2) * disp(k3) * disp(-1k1 + k3 - kp + q) -
                (2 * disp(-1k + k1 + k2) * disp(k3) * disp(-1k1 + k3 - kp + q)) +
                2 * disp(k3) ^ 2 * disp(-1k1 + k3 - kp + q) +
                4 * mu * disp(-1k1 + k3 + kp) * disp(-1k1 + k3 - kp + q) +
                2 * w * disp(-1k1 + k3 + kp) * disp(-1k1 + k3 - kp + q) -
                (2 * disp(k3) * disp(-1k1 + k3 + kp) * disp(-1k1 + k3 - kp + q)) -
                (2 * v * disp(-1k3 + q) * disp(-1k1 + k3 - kp + q)) +
                2 * vp * disp(-1k3 + q) * disp(-1k1 + k3 - kp + q) +
                2 * disp(k2) * disp(-1k3 + q) * disp(-1k1 + k3 - kp + q) -
                (2 * disp(-1k + k1 + k2) * disp(-1k3 + q) * disp(-1k1 + k3 - kp + q)) +
                2 * disp(k3) * disp(-1k3 + q) * disp(-1k1 + k3 - kp + q) -
                (2 * disp(-1k1 + k3 + kp) * disp(-1k3 + q) * disp(-1k1 + k3 - kp + q)) +
                4 * mu * u * v * fermidist(-1mu + disp(k3), beta) -
                (2 * u * v ^ 2 * fermidist(-1mu + disp(k3), beta)) +
                2 * u * vp ^ 2 * fermidist(-1mu + disp(k3), beta) -
                (2 * mu * u * w * fermidist(-1mu + disp(k3), beta)) +
                4 * u * v * w * fermidist(-1mu + disp(k3), beta) -
                (2 * u * vp * w * fermidist(-1mu + disp(k3), beta)) -
                (u * w ^ 2 * fermidist(-1mu + disp(k3), beta)) -
                (4 * mu * u * disp(k2) * fermidist(-1mu + disp(k3), beta)) +
                4 * u * v * disp(k2) * fermidist(-1mu + disp(k3), beta) -
                (4 * u * w * disp(k2) * fermidist(-1mu + disp(k3), beta)) -
                (2 * u * disp(k2) ^ 2 * fermidist(-1mu + disp(k3), beta)) +
                4 * mu * u * disp(-1k + k1 + k2) * fermidist(-1mu + disp(k3), beta) -
                (4 * u * v * disp(-1k + k1 + k2) * fermidist(-1mu + disp(k3), beta)) +
                4 * u * w * disp(-1k + k1 + k2) * fermidist(-1mu + disp(k3), beta) +
                4 * u * disp(k2) * disp(-1k + k1 + k2) * fermidist(-1mu + disp(k3), beta) -
                (2 * u * disp(-1k + k1 + k2) ^ 2 * fermidist(-1mu + disp(k3), beta)) -
                (4 * mu * u * disp(k3) * fermidist(-1mu + disp(k3), beta)) +
                2 * u * v * disp(k3) * fermidist(-1mu + disp(k3), beta) -
                (3 * u * w * disp(k3) * fermidist(-1mu + disp(k3), beta)) -
                (2 * u * disp(k2) * disp(k3) * fermidist(-1mu + disp(k3), beta)) +
                2 * u * disp(-1k + k1 + k2) * disp(k3) * fermidist(-1mu + disp(k3), beta) +
                2 * mu * u * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta) -
                (2 * u * v * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) -
                (2 * u * vp * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) +
                3 * u * w * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta) +
                2 * u * disp(k2) * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta) -
                (
                    2 *
                    u *
                    disp(-1k + k1 + k2) *
                    disp(-1k1 + k3 + kp) *
                    fermidist(-1mu + disp(k3), beta)
                ) + u * disp(k3) * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta) -
                (2 * u * v * disp(-1k3 + q) * fermidist(-1mu + disp(k3), beta)) +
                u * w * disp(-1k3 + q) * fermidist(-1mu + disp(k3), beta) +
                2 * u * disp(k2) * disp(-1k3 + q) * fermidist(-1mu + disp(k3), beta) - (
                    2 *
                    u *
                    disp(-1k + k1 + k2) *
                    disp(-1k3 + q) *
                    fermidist(-1mu + disp(k3), beta)
                ) + 2 * u * disp(k3) * disp(-1k3 + q) * fermidist(-1mu + disp(k3), beta) -
                (
                    u *
                    disp(-1k1 + k3 + kp) *
                    disp(-1k3 + q) *
                    fermidist(-1mu + disp(k3), beta)
                ) +
                2 * mu * u * disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta) -
                (2 * u * v * disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta)) +
                2 * u * vp * disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta) +
                u * w * disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta) +
                2 *
                u *
                disp(k2) *
                disp(-1k1 + k3 - kp + q) *
                fermidist(-1mu + disp(k3), beta) - (
                    2 *
                    u *
                    disp(-1k + k1 + k2) *
                    disp(-1k1 + k3 - kp + q) *
                    fermidist(-1mu + disp(k3), beta)
                ) +
                u * disp(k3) * disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta) -
                (
                    2 *
                    u *
                    disp(-1k1 + k3 + kp) *
                    disp(-1k1 + k3 - kp + q) *
                    fermidist(-1mu + disp(k3), beta)
                ) - (
                    u *
                    disp(-1k3 + q) *
                    disp(-1k1 + k3 - kp + q) *
                    fermidist(-1mu + disp(k3), beta)
                ) - (2 * mu * u * v * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
                (2 * mu * u * vp * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
                2 * mu * u * w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
                (u * v * w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
                (u * vp * w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
                u * w ^ 2 * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
                2 * mu * u * disp(k2) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
                u * w * disp(k2) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) - (
                    2 *
                    mu *
                    u *
                    disp(-1k + k1 + k2) *
                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                ) - (
                    u *
                    w *
                    disp(-1k + k1 + k2) *
                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                ) +
                2 * mu * u * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
                u * v * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
                u * vp * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
                (u * disp(k2) * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
                u *
                disp(-1k + k1 + k2) *
                disp(k3) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
                (u * disp(k3) ^ 2 * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
                u * v * disp(-1k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
                u * vp * disp(-1k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
                (u * w * disp(-1k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
                (
                    u *
                    disp(k2) *
                    disp(-1k3 + q) *
                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                ) +
                u *
                disp(-1k + k1 + k2) *
                disp(-1k3 + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta) - (
                    u *
                    disp(k3) *
                    disp(-1k3 + q) *
                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                ) - (
                    2 *
                    mu *
                    u *
                    disp(-1k1 + k3 - kp + q) *
                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                ) - (
                    u *
                    w *
                    disp(-1k1 + k3 - kp + q) *
                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                ) +
                u *
                disp(k3) *
                disp(-1k1 + k3 - kp + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
                u *
                disp(-1k3 + q) *
                disp(-1k1 + k3 - kp + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
                2 * u * v ^ 2 * fermidist(mu - disp(-1k3 + q), beta) -
                (2 * u * vp ^ 2 * fermidist(mu - disp(-1k3 + q), beta)) -
                (2 * u * v * w * fermidist(mu - disp(-1k3 + q), beta)) +
                2 * u * vp * w * fermidist(mu - disp(-1k3 + q), beta) -
                (4 * u * v * disp(k2) * fermidist(mu - disp(-1k3 + q), beta)) +
                2 * u * w * disp(k2) * fermidist(mu - disp(-1k3 + q), beta) +
                2 * u * disp(k2) ^ 2 * fermidist(mu - disp(-1k3 + q), beta) +
                4 * u * v * disp(-1k + k1 + k2) * fermidist(mu - disp(-1k3 + q), beta) -
                (2 * u * w * disp(-1k + k1 + k2) * fermidist(mu - disp(-1k3 + q), beta)) -
                (
                    4 *
                    u *
                    disp(k2) *
                    disp(-1k + k1 + k2) *
                    fermidist(mu - disp(-1k3 + q), beta)
                ) + 2 * u * disp(-1k + k1 + k2) ^ 2 * fermidist(mu - disp(-1k3 + q), beta) -
                (4 * u * v * disp(k3) * fermidist(mu - disp(-1k3 + q), beta)) +
                2 * u * w * disp(k3) * fermidist(mu - disp(-1k3 + q), beta) +
                4 * u * disp(k2) * disp(k3) * fermidist(mu - disp(-1k3 + q), beta) - (
                    4 *
                    u *
                    disp(-1k + k1 + k2) *
                    disp(k3) *
                    fermidist(mu - disp(-1k3 + q), beta)
                ) +
                2 * u * disp(k3) ^ 2 * fermidist(mu - disp(-1k3 + q), beta) +
                2 * u * v * disp(-1k1 + k3 + kp) * fermidist(mu - disp(-1k3 + q), beta) +
                2 * u * vp * disp(-1k1 + k3 + kp) * fermidist(mu - disp(-1k3 + q), beta) -
                (2 * u * w * disp(-1k1 + k3 + kp) * fermidist(mu - disp(-1k3 + q), beta)) -
                (
                    2 *
                    u *
                    disp(k2) *
                    disp(-1k1 + k3 + kp) *
                    fermidist(mu - disp(-1k3 + q), beta)
                ) +
                2 *
                u *
                disp(-1k + k1 + k2) *
                disp(-1k1 + k3 + kp) *
                fermidist(mu - disp(-1k3 + q), beta) - (
                    2 *
                    u *
                    disp(k3) *
                    disp(-1k1 + k3 + kp) *
                    fermidist(mu - disp(-1k3 + q), beta)
                ) +
                2 *
                u *
                v *
                disp(-1k1 + k3 - kp + q) *
                fermidist(mu - disp(-1k3 + q), beta) - (
                    2 *
                    u *
                    vp *
                    disp(-1k1 + k3 - kp + q) *
                    fermidist(mu - disp(-1k3 + q), beta)
                ) - (
                    2 *
                    u *
                    disp(k2) *
                    disp(-1k1 + k3 - kp + q) *
                    fermidist(mu - disp(-1k3 + q), beta)
                ) +
                2 *
                u *
                disp(-1k + k1 + k2) *
                disp(-1k1 + k3 - kp + q) *
                fermidist(mu - disp(-1k3 + q), beta) - (
                    2 *
                    u *
                    disp(k3) *
                    disp(-1k1 + k3 - kp + q) *
                    fermidist(mu - disp(-1k3 + q), beta)
                ) +
                2 *
                u *
                disp(-1k1 + k3 + kp) *
                disp(-1k1 + k3 - kp + q) *
                fermidist(mu - disp(-1k3 + q), beta) -
                (2 * mu * u * v * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
                2 * mu * u * vp * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
                (u * v * w * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
                u * vp * w * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) +
                2 * mu * u * disp(k2) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) +
                u * w * disp(k2) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) - (
                    2 *
                    mu *
                    u *
                    disp(-1k + k1 + k2) *
                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                ) - (
                    u *
                    w *
                    disp(-1k + k1 + k2) *
                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                ) +
                2 * mu * u * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) +
                u * v * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
                (u * vp * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
                u * w * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) - (
                    u *
                    disp(k2) *
                    disp(k3) *
                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                ) +
                u *
                disp(-1k + k1 + k2) *
                disp(k3) *
                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
                (u * disp(k3) ^ 2 * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) - (
                    2 *
                    mu *
                    u *
                    disp(-1k1 + k3 + kp) *
                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                ) - (
                    u *
                    w *
                    disp(-1k1 + k3 + kp) *
                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                ) +
                u *
                disp(k3) *
                disp(-1k1 + k3 + kp) *
                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) +
                u * v * disp(-1k3 + q) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
                (
                    u *
                    vp *
                    disp(-1k3 + q) *
                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                ) - (
                    u *
                    disp(k2) *
                    disp(-1k3 + q) *
                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                ) +
                u *
                disp(-1k + k1 + k2) *
                disp(-1k3 + q) *
                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) - (
                    u *
                    disp(k3) *
                    disp(-1k3 + q) *
                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                ) +
                u *
                disp(-1k1 + k3 + kp) *
                disp(-1k3 + q) *
                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
            )
        ) *
        (
            2 *
            (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) *
            (v - vp - disp(k2) + disp(-1k + k1 + k2) - disp(k3) + disp(-1k1 + k3 + kp)) *
            (mu - v + w + disp(k2) - disp(-1k + k1 + k2) - disp(-1k1 + q)) *
            (2v - w - (2 * disp(k2)) + disp(-1k + k1 + k2) + disp(-1k - k1 + k2 + q)) *
            (2mu + w - disp(k3) - disp(-1k3 + q)) *
            (
                v + vp - w - disp(k2) + disp(-1k + k1 + k2) - disp(k3) +
                disp(-1k1 + k3 - kp + q)
            )
        ) ^ -1
    ) +
    (
        u ^ 2 *
        fermidist(mu - disp(-1k1 + q), beta) *
        (
            2 * mu ^ 2 - (2 * v ^ 2) + 2 * mu * w + 2 * v * w + 4 * v * disp(k2) -
            (2 * w * disp(k2)) - (2 * disp(k2) ^ 2) - (2 * mu * disp(-1k + k1 + k2)) -
            (2 * v * disp(-1k + k1 + k2)) + 2 * disp(k2) * disp(-1k + k1 + k2) -
            (4 * mu * disp(-1k1 + q)) - (2 * w * disp(-1k1 + q)) +
            2 * disp(-1k + k1 + k2) * disp(-1k1 + q) +
            2 * disp(-1k1 + q) ^ 2 +
            2 * mu * disp(-1k - k1 + k2 + q) - (2 * v * disp(-1k - k1 + k2 + q)) +
            2 * w * disp(-1k - k1 + k2 + q) +
            2 * disp(k2) * disp(-1k - k1 + k2 + q) -
            (2 * disp(-1k + k1 + k2) * disp(-1k - k1 + k2 + q)) -
            (2 * disp(-1k1 + q) * disp(-1k - k1 + k2 + q)) -
            (2 * u * v * fermidist(-1mu + disp(k2), beta)) +
            u * w * fermidist(-1mu + disp(k2), beta) +
            2 * u * disp(k2) * fermidist(-1mu + disp(k2), beta) -
            (u * disp(-1k + k1 + k2) * fermidist(-1mu + disp(k2), beta)) -
            (u * disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(k2), beta)) +
            mu * u * fermidist(-1mu + disp(-1k + k1 + k2), beta) +
            u * v * fermidist(-1mu + disp(-1k + k1 + k2), beta) -
            (u * disp(k2) * fermidist(-1mu + disp(-1k + k1 + k2), beta)) -
            (u * disp(-1k1 + q) * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
            u * disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(-1k + k1 + k2), beta) -
            (mu * u * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) +
            u * v * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta) -
            (u * w * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) -
            (u * disp(k2) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) +
            u * disp(-1k + k1 + k2) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta) +
            u * disp(-1k1 + q) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
        ) *
        (
            4 * mu ^ 3 - (4 * mu * vp ^ 2) + 6 * mu ^ 2 * w + 4 * mu * vp * w -
            (2 * vp ^ 2 * w) +
            2 * mu * w ^ 2 +
            2 * vp * w ^ 2 - (10 * mu ^ 2 * disp(k3)) + 2 * vp ^ 2 * disp(k3) -
            (10 * mu * w * disp(k3)) - (2 * vp * w * disp(k3)) - (2 * w ^ 2 * disp(k3)) +
            8 * mu * disp(k3) ^ 2 +
            4 * w * disp(k3) ^ 2 - (2 * disp(k3) ^ 3) +
            4 * mu ^ 2 * disp(-1k1 + k3 + kp) +
            4 * mu * vp * disp(-1k1 + k3 + kp) +
            2 * mu * w * disp(-1k1 + k3 + kp) +
            2 * vp * w * disp(-1k1 + k3 + kp) - (6 * mu * disp(k3) * disp(-1k1 + k3 + kp)) -
            (2 * vp * disp(k3) * disp(-1k1 + k3 + kp)) -
            (2 * w * disp(k3) * disp(-1k1 + k3 + kp)) +
            2 * disp(k3) ^ 2 * disp(-1k1 + k3 + kp) - (8 * mu ^ 2 * disp(-1k1 + q)) -
            (8 * mu * w * disp(-1k1 + q)) - (2 * w ^ 2 * disp(-1k1 + q)) +
            12 * mu * disp(k3) * disp(-1k1 + q) +
            6 * w * disp(k3) * disp(-1k1 + q) - (4 * disp(k3) ^ 2 * disp(-1k1 + q)) -
            (4 * mu * disp(-1k1 + k3 + kp) * disp(-1k1 + q)) -
            (2 * w * disp(-1k1 + k3 + kp) * disp(-1k1 + q)) +
            2 * disp(k3) * disp(-1k1 + k3 + kp) * disp(-1k1 + q) +
            4 * mu * disp(-1k1 + q) ^ 2 +
            2 * w * disp(-1k1 + q) ^ 2 - (2 * disp(k3) * disp(-1k1 + q) ^ 2) -
            (2 * mu ^ 2 * disp(-1k3 + q)) + 2 * vp ^ 2 * disp(-1k3 + q) -
            (2 * mu * w * disp(-1k3 + q)) - (2 * vp * w * disp(-1k3 + q)) +
            4 * mu * disp(k3) * disp(-1k3 + q) +
            2 * w * disp(k3) * disp(-1k3 + q) - (2 * disp(k3) ^ 2 * disp(-1k3 + q)) -
            (2 * mu * disp(-1k1 + k3 + kp) * disp(-1k3 + q)) -
            (2 * vp * disp(-1k1 + k3 + kp) * disp(-1k3 + q)) +
            2 * disp(k3) * disp(-1k1 + k3 + kp) * disp(-1k3 + q) +
            4 * mu * disp(-1k1 + q) * disp(-1k3 + q) +
            2 * w * disp(-1k1 + q) * disp(-1k3 + q) -
            (4 * disp(k3) * disp(-1k1 + q) * disp(-1k3 + q)) +
            2 * disp(-1k1 + k3 + kp) * disp(-1k1 + q) * disp(-1k3 + q) -
            (2 * disp(-1k1 + q) ^ 2 * disp(-1k3 + q)) +
            4 * mu ^ 2 * disp(-1k1 + k3 - kp + q) -
            (4 * mu * vp * disp(-1k1 + k3 - kp + q)) +
            6 * mu * w * disp(-1k1 + k3 - kp + q) -
            (2 * vp * w * disp(-1k1 + k3 - kp + q)) + 2 * w ^ 2 * disp(-1k1 + k3 - kp + q) -
            (6 * mu * disp(k3) * disp(-1k1 + k3 - kp + q)) +
            2 * vp * disp(k3) * disp(-1k1 + k3 - kp + q) -
            (4 * w * disp(k3) * disp(-1k1 + k3 - kp + q)) +
            2 * disp(k3) ^ 2 * disp(-1k1 + k3 - kp + q) +
            4 * mu * disp(-1k1 + k3 + kp) * disp(-1k1 + k3 - kp + q) +
            2 * w * disp(-1k1 + k3 + kp) * disp(-1k1 + k3 - kp + q) -
            (2 * disp(k3) * disp(-1k1 + k3 + kp) * disp(-1k1 + k3 - kp + q)) -
            (4 * mu * disp(-1k1 + q) * disp(-1k1 + k3 - kp + q)) -
            (2 * w * disp(-1k1 + q) * disp(-1k1 + k3 - kp + q)) +
            2 * disp(k3) * disp(-1k1 + q) * disp(-1k1 + k3 - kp + q) -
            (2 * mu * disp(-1k3 + q) * disp(-1k1 + k3 - kp + q)) +
            2 * vp * disp(-1k3 + q) * disp(-1k1 + k3 - kp + q) -
            (2 * w * disp(-1k3 + q) * disp(-1k1 + k3 - kp + q)) +
            2 * disp(k3) * disp(-1k3 + q) * disp(-1k1 + k3 - kp + q) -
            (2 * disp(-1k1 + k3 + kp) * disp(-1k3 + q) * disp(-1k1 + k3 - kp + q)) +
            2 * disp(-1k1 + q) * disp(-1k3 + q) * disp(-1k1 + k3 - kp + q) +
            2 * mu ^ 2 * u * fermidist(-1mu + disp(k3), beta) +
            2 * u * vp ^ 2 * fermidist(-1mu + disp(k3), beta) +
            2 * mu * u * w * fermidist(-1mu + disp(k3), beta) -
            (2 * u * vp * w * fermidist(-1mu + disp(k3), beta)) +
            u * w ^ 2 * fermidist(-1mu + disp(k3), beta) -
            (2 * mu * u * disp(k3) * fermidist(-1mu + disp(k3), beta)) -
            (u * w * disp(k3) * fermidist(-1mu + disp(k3), beta)) -
            (2 * u * vp * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) +
            u * w * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta) +
            u * disp(k3) * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta) -
            (2 * u * disp(k3) * disp(-1k1 + q) * fermidist(-1mu + disp(k3), beta)) +
            2 *
            u *
            disp(-1k1 + k3 + kp) *
            disp(-1k1 + q) *
            fermidist(-1mu + disp(k3), beta) -
            (2 * u * disp(-1k1 + q) ^ 2 * fermidist(-1mu + disp(k3), beta)) -
            (2 * mu * u * disp(-1k3 + q) * fermidist(-1mu + disp(k3), beta)) -
            (u * w * disp(-1k3 + q) * fermidist(-1mu + disp(k3), beta)) +
            2 * u * disp(k3) * disp(-1k3 + q) * fermidist(-1mu + disp(k3), beta) -
            (u * disp(-1k1 + k3 + kp) * disp(-1k3 + q) * fermidist(-1mu + disp(k3), beta)) +
            2 * u * disp(-1k1 + q) * disp(-1k3 + q) * fermidist(-1mu + disp(k3), beta) +
            2 * u * vp * disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta) -
            (u * w * disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta)) +
            u * disp(k3) * disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta) - (
                2 *
                u *
                disp(-1k1 + k3 + kp) *
                disp(-1k1 + k3 - kp + q) *
                fermidist(-1mu + disp(k3), beta)
            ) +
            2 *
            u *
            disp(-1k1 + q) *
            disp(-1k1 + k3 - kp + q) *
            fermidist(-1mu + disp(k3), beta) - (
                u *
                disp(-1k3 + q) *
                disp(-1k1 + k3 - kp + q) *
                fermidist(-1mu + disp(k3), beta)
            ) - (2 * mu ^ 2 * u * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
            (2 * mu * u * vp * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
            (mu * u * w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
            (u * vp * w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
            3 * mu * u * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            u * vp * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            u * w * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
            (u * disp(k3) ^ 2 * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
            2 * mu * u * disp(-1k1 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            u * w * disp(-1k1 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
            (u * disp(k3) * disp(-1k1 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
            mu * u * disp(-1k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            u * vp * disp(-1k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
            (u * disp(k3) * disp(-1k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
            (
                u *
                disp(-1k1 + q) *
                disp(-1k3 + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
            ) - (
                2 *
                mu *
                u *
                disp(-1k1 + k3 - kp + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
            ) - (
                u *
                w *
                disp(-1k1 + k3 - kp + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
            ) +
            u *
            disp(k3) *
            disp(-1k1 + k3 - kp + q) *
            fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            u *
            disp(-1k3 + q) *
            disp(-1k1 + k3 - kp + q) *
            fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            2 * mu ^ 2 * u * fermidist(mu - disp(-1k3 + q), beta) -
            (2 * u * vp ^ 2 * fermidist(mu - disp(-1k3 + q), beta)) +
            2 * mu * u * w * fermidist(mu - disp(-1k3 + q), beta) +
            2 * u * vp * w * fermidist(mu - disp(-1k3 + q), beta) -
            (4 * mu * u * disp(k3) * fermidist(mu - disp(-1k3 + q), beta)) -
            (2 * u * w * disp(k3) * fermidist(mu - disp(-1k3 + q), beta)) +
            2 * u * disp(k3) ^ 2 * fermidist(mu - disp(-1k3 + q), beta) +
            2 * mu * u * disp(-1k1 + k3 + kp) * fermidist(mu - disp(-1k3 + q), beta) +
            2 * u * vp * disp(-1k1 + k3 + kp) * fermidist(mu - disp(-1k3 + q), beta) - (
                2 *
                u *
                disp(k3) *
                disp(-1k1 + k3 + kp) *
                fermidist(mu - disp(-1k3 + q), beta)
            ) - (4 * mu * u * disp(-1k1 + q) * fermidist(mu - disp(-1k3 + q), beta)) -
            (2 * u * w * disp(-1k1 + q) * fermidist(mu - disp(-1k3 + q), beta)) +
            4 * u * disp(k3) * disp(-1k1 + q) * fermidist(mu - disp(-1k3 + q), beta) - (
                2 *
                u *
                disp(-1k1 + k3 + kp) *
                disp(-1k1 + q) *
                fermidist(mu - disp(-1k3 + q), beta)
            ) +
            2 * u * disp(-1k1 + q) ^ 2 * fermidist(mu - disp(-1k3 + q), beta) +
            2 * mu * u * disp(-1k1 + k3 - kp + q) * fermidist(mu - disp(-1k3 + q), beta) -
            (2 * u * vp * disp(-1k1 + k3 - kp + q) * fermidist(mu - disp(-1k3 + q), beta)) +
            2 * u * w * disp(-1k1 + k3 - kp + q) * fermidist(mu - disp(-1k3 + q), beta) -
            (
                2 *
                u *
                disp(k3) *
                disp(-1k1 + k3 - kp + q) *
                fermidist(mu - disp(-1k3 + q), beta)
            ) +
            2 *
            u *
            disp(-1k1 + k3 + kp) *
            disp(-1k1 + k3 - kp + q) *
            fermidist(mu - disp(-1k3 + q), beta) - (
                2 *
                u *
                disp(-1k1 + q) *
                disp(-1k1 + k3 - kp + q) *
                fermidist(mu - disp(-1k3 + q), beta)
            ) - (2 * mu ^ 2 * u * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
            2 * mu * u * vp * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
            (3 * mu * u * w * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
            u * vp * w * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
            (u * w ^ 2 * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
            3 * mu * u * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
            (u * vp * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
            2 * u * w * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
            (u * disp(k3) ^ 2 * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) - (
                2 *
                mu *
                u *
                disp(-1k1 + k3 + kp) *
                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
            ) - (
                u *
                w *
                disp(-1k1 + k3 + kp) *
                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
            ) +
            u *
            disp(k3) *
            disp(-1k1 + k3 + kp) *
            fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) +
            2 * mu * u * disp(-1k1 + q) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) +
            u * w * disp(-1k1 + q) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) - (
                u *
                disp(k3) *
                disp(-1k1 + q) *
                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
            ) + mu * u * disp(-1k3 + q) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
            (u * vp * disp(-1k3 + q) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
            u * w * disp(-1k3 + q) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) - (
                u *
                disp(k3) *
                disp(-1k3 + q) *
                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
            ) +
            u *
            disp(-1k1 + k3 + kp) *
            disp(-1k3 + q) *
            fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) - (
                u *
                disp(-1k1 + q) *
                disp(-1k3 + q) *
                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
            )
        )
    ) *
    (
        2 *
        (2mu + w - disp(k1) - disp(-1k1 + q)) *
        (mu - v + w + disp(k2) - disp(-1k + k1 + k2) - disp(-1k1 + q)) *
        (mu - vp + w - disp(k3) + disp(-1k1 + k3 + kp) - disp(-1k1 + q)) *
        (mu + v - disp(k2) - disp(-1k1 + q) + disp(-1k - k1 + k2 + q)) *
        (2mu + w - disp(k3) - disp(-1k3 + q)) *
        (mu + vp - disp(k3) - disp(-1k1 + q) + disp(-1k1 + k3 - kp + q))
    ) ^ -1 - (
        (
            u ^ 2 *
            bosedist(disp(k2) - disp(-1k - k1 + k2 + q), beta) *
            (
                2 * u * v * fermidist(-1mu + disp(k2), beta) -
                (u * w * fermidist(-1mu + disp(k2), beta)) -
                (2 * u * disp(k2) * fermidist(-1mu + disp(k2), beta)) +
                u * disp(-1k + k1 + k2) * fermidist(-1mu + disp(k2), beta) +
                u * disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(k2), beta) -
                (2 * u * v * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) +
                u * w * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta) +
                2 * u * disp(k2) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta) - (
                    u *
                    disp(-1k + k1 + k2) *
                    fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
                ) - (
                    u *
                    disp(-1k - k1 + k2 + q) *
                    fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
                )
            ) *
            (
                -4 * mu * v ^ 2 + 4 * mu * vp ^ 2 + 4 * mu * v * w - (2 * v ^ 2 * w) -
                (4 * mu * vp * w) +
                2 * vp ^ 2 * w +
                2 * v * w ^ 2 - (2 * vp * w ^ 2) + 8 * mu * v * disp(k2) -
                (4 * mu * w * disp(k2)) + 4 * v * w * disp(k2) - (2 * w ^ 2 * disp(k2)) -
                (4 * mu * disp(k2) ^ 2) - (2 * w * disp(k2) ^ 2) - (8 * mu * v * disp(k3)) +
                2 * v ^ 2 * disp(k3) - (2 * vp ^ 2 * disp(k3)) + 4 * mu * w * disp(k3) -
                (6 * v * w * disp(k3)) +
                2 * vp * w * disp(k3) +
                2 * w ^ 2 * disp(k3) +
                8 * mu * disp(k2) * disp(k3) - (4 * v * disp(k2) * disp(k3)) +
                6 * w * disp(k2) * disp(k3) +
                2 * disp(k2) ^ 2 * disp(k3) - (4 * mu * disp(k3) ^ 2) +
                4 * v * disp(k3) ^ 2 - (4 * w * disp(k3) ^ 2) -
                (4 * disp(k2) * disp(k3) ^ 2) +
                2 * disp(k3) ^ 3 +
                4 * mu * v * disp(-1k1 + k3 + kp) - (4 * mu * vp * disp(-1k1 + k3 + kp)) +
                2 * v * w * disp(-1k1 + k3 + kp) - (2 * vp * w * disp(-1k1 + k3 + kp)) -
                (4 * mu * disp(k2) * disp(-1k1 + k3 + kp)) -
                (2 * w * disp(k2) * disp(-1k1 + k3 + kp)) +
                4 * mu * disp(k3) * disp(-1k1 + k3 + kp) -
                (2 * v * disp(k3) * disp(-1k1 + k3 + kp)) +
                2 * vp * disp(k3) * disp(-1k1 + k3 + kp) +
                2 * w * disp(k3) * disp(-1k1 + k3 + kp) +
                2 * disp(k2) * disp(k3) * disp(-1k1 + k3 + kp) -
                (2 * disp(k3) ^ 2 * disp(-1k1 + k3 + kp)) -
                (8 * mu * v * disp(-1k - k1 + k2 + q)) +
                4 * mu * w * disp(-1k - k1 + k2 + q) -
                (4 * v * w * disp(-1k - k1 + k2 + q)) +
                2 * w ^ 2 * disp(-1k - k1 + k2 + q) +
                8 * mu * disp(k2) * disp(-1k - k1 + k2 + q) +
                4 * w * disp(k2) * disp(-1k - k1 + k2 + q) -
                (8 * mu * disp(k3) * disp(-1k - k1 + k2 + q)) +
                4 * v * disp(k3) * disp(-1k - k1 + k2 + q) -
                (6 * w * disp(k3) * disp(-1k - k1 + k2 + q)) -
                (4 * disp(k2) * disp(k3) * disp(-1k - k1 + k2 + q)) +
                4 * disp(k3) ^ 2 * disp(-1k - k1 + k2 + q) +
                4 * mu * disp(-1k1 + k3 + kp) * disp(-1k - k1 + k2 + q) +
                2 * w * disp(-1k1 + k3 + kp) * disp(-1k - k1 + k2 + q) -
                (2 * disp(k3) * disp(-1k1 + k3 + kp) * disp(-1k - k1 + k2 + q)) -
                (4 * mu * disp(-1k - k1 + k2 + q) ^ 2) -
                (2 * w * disp(-1k - k1 + k2 + q) ^ 2) +
                2 * disp(k3) * disp(-1k - k1 + k2 + q) ^ 2 +
                2 * v ^ 2 * disp(-1k3 + q) - (2 * vp ^ 2 * disp(-1k3 + q)) -
                (2 * v * w * disp(-1k3 + q)) + 2 * vp * w * disp(-1k3 + q) -
                (4 * v * disp(k2) * disp(-1k3 + q)) +
                2 * w * disp(k2) * disp(-1k3 + q) +
                2 * disp(k2) ^ 2 * disp(-1k3 + q) +
                4 * v * disp(k3) * disp(-1k3 + q) - (2 * w * disp(k3) * disp(-1k3 + q)) -
                (4 * disp(k2) * disp(k3) * disp(-1k3 + q)) +
                2 * disp(k3) ^ 2 * disp(-1k3 + q) -
                (2 * v * disp(-1k1 + k3 + kp) * disp(-1k3 + q)) +
                2 * vp * disp(-1k1 + k3 + kp) * disp(-1k3 + q) +
                2 * disp(k2) * disp(-1k1 + k3 + kp) * disp(-1k3 + q) -
                (2 * disp(k3) * disp(-1k1 + k3 + kp) * disp(-1k3 + q)) +
                4 * v * disp(-1k - k1 + k2 + q) * disp(-1k3 + q) -
                (2 * w * disp(-1k - k1 + k2 + q) * disp(-1k3 + q)) -
                (4 * disp(k2) * disp(-1k - k1 + k2 + q) * disp(-1k3 + q)) +
                4 * disp(k3) * disp(-1k - k1 + k2 + q) * disp(-1k3 + q) -
                (2 * disp(-1k1 + k3 + kp) * disp(-1k - k1 + k2 + q) * disp(-1k3 + q)) +
                2 * disp(-1k - k1 + k2 + q) ^ 2 * disp(-1k3 + q) +
                4 * mu * v * disp(-1k1 + k3 - kp + q) +
                4 * mu * vp * disp(-1k1 + k3 - kp + q) -
                (4 * mu * w * disp(-1k1 + k3 - kp + q)) +
                2 * v * w * disp(-1k1 + k3 - kp + q) +
                2 * vp * w * disp(-1k1 + k3 - kp + q) -
                (2 * w ^ 2 * disp(-1k1 + k3 - kp + q)) -
                (4 * mu * disp(k2) * disp(-1k1 + k3 - kp + q)) -
                (2 * w * disp(k2) * disp(-1k1 + k3 - kp + q)) +
                4 * mu * disp(k3) * disp(-1k1 + k3 - kp + q) -
                (2 * v * disp(k3) * disp(-1k1 + k3 - kp + q)) -
                (2 * vp * disp(k3) * disp(-1k1 + k3 - kp + q)) +
                4 * w * disp(k3) * disp(-1k1 + k3 - kp + q) +
                2 * disp(k2) * disp(k3) * disp(-1k1 + k3 - kp + q) -
                (2 * disp(k3) ^ 2 * disp(-1k1 + k3 - kp + q)) -
                (4 * mu * disp(-1k1 + k3 + kp) * disp(-1k1 + k3 - kp + q)) -
                (2 * w * disp(-1k1 + k3 + kp) * disp(-1k1 + k3 - kp + q)) +
                2 * disp(k3) * disp(-1k1 + k3 + kp) * disp(-1k1 + k3 - kp + q) +
                4 * mu * disp(-1k - k1 + k2 + q) * disp(-1k1 + k3 - kp + q) +
                2 * w * disp(-1k - k1 + k2 + q) * disp(-1k1 + k3 - kp + q) -
                (2 * disp(k3) * disp(-1k - k1 + k2 + q) * disp(-1k1 + k3 - kp + q)) -
                (2 * v * disp(-1k3 + q) * disp(-1k1 + k3 - kp + q)) -
                (2 * vp * disp(-1k3 + q) * disp(-1k1 + k3 - kp + q)) +
                2 * w * disp(-1k3 + q) * disp(-1k1 + k3 - kp + q) +
                2 * disp(k2) * disp(-1k3 + q) * disp(-1k1 + k3 - kp + q) -
                (2 * disp(k3) * disp(-1k3 + q) * disp(-1k1 + k3 - kp + q)) +
                2 * disp(-1k1 + k3 + kp) * disp(-1k3 + q) * disp(-1k1 + k3 - kp + q) -
                (2 * disp(-1k - k1 + k2 + q) * disp(-1k3 + q) * disp(-1k1 + k3 - kp + q)) +
                4 * mu * u * v * fermidist(-1mu + disp(k3), beta) +
                2 * u * v ^ 2 * fermidist(-1mu + disp(k3), beta) -
                (2 * u * vp ^ 2 * fermidist(-1mu + disp(k3), beta)) -
                (2 * mu * u * w * fermidist(-1mu + disp(k3), beta)) +
                2 * u * vp * w * fermidist(-1mu + disp(k3), beta) -
                (u * w ^ 2 * fermidist(-1mu + disp(k3), beta)) -
                (4 * mu * u * disp(k2) * fermidist(-1mu + disp(k3), beta)) -
                (4 * u * v * disp(k2) * fermidist(-1mu + disp(k3), beta)) +
                2 * u * disp(k2) ^ 2 * fermidist(-1mu + disp(k3), beta) +
                4 * mu * u * disp(k3) * fermidist(-1mu + disp(k3), beta) +
                2 * u * v * disp(k3) * fermidist(-1mu + disp(k3), beta) +
                u * w * disp(k3) * fermidist(-1mu + disp(k3), beta) -
                (2 * u * disp(k2) * disp(k3) * fermidist(-1mu + disp(k3), beta)) -
                (2 * mu * u * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) -
                (2 * u * v * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) +
                2 * u * vp * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta) -
                (u * w * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) +
                2 * u * disp(k2) * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta) -
                (u * disp(k3) * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) +
                4 * mu * u * disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(k3), beta) +
                4 * u * v * disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(k3), beta) - (
                    4 *
                    u *
                    disp(k2) *
                    disp(-1k - k1 + k2 + q) *
                    fermidist(-1mu + disp(k3), beta)
                ) +
                2 *
                u *
                disp(k3) *
                disp(-1k - k1 + k2 + q) *
                fermidist(-1mu + disp(k3), beta) - (
                    2 *
                    u *
                    disp(-1k1 + k3 + kp) *
                    disp(-1k - k1 + k2 + q) *
                    fermidist(-1mu + disp(k3), beta)
                ) + 2 * u * disp(-1k - k1 + k2 + q) ^ 2 * fermidist(-1mu + disp(k3), beta) -
                (2 * u * v * disp(-1k3 + q) * fermidist(-1mu + disp(k3), beta)) +
                u * w * disp(-1k3 + q) * fermidist(-1mu + disp(k3), beta) +
                2 * u * disp(k2) * disp(-1k3 + q) * fermidist(-1mu + disp(k3), beta) -
                (2 * u * disp(k3) * disp(-1k3 + q) * fermidist(-1mu + disp(k3), beta)) +
                u *
                disp(-1k1 + k3 + kp) *
                disp(-1k3 + q) *
                fermidist(-1mu + disp(k3), beta) - (
                    2 *
                    u *
                    disp(-1k - k1 + k2 + q) *
                    disp(-1k3 + q) *
                    fermidist(-1mu + disp(k3), beta)
                ) -
                (2 * mu * u * disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta)) -
                (2 * u * v * disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta)) -
                (2 * u * vp * disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta)) +
                u * w * disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta) +
                2 *
                u *
                disp(k2) *
                disp(-1k1 + k3 - kp + q) *
                fermidist(-1mu + disp(k3), beta) - (
                    u *
                    disp(k3) *
                    disp(-1k1 + k3 - kp + q) *
                    fermidist(-1mu + disp(k3), beta)
                ) +
                2 *
                u *
                disp(-1k1 + k3 + kp) *
                disp(-1k1 + k3 - kp + q) *
                fermidist(-1mu + disp(k3), beta) - (
                    2 *
                    u *
                    disp(-1k - k1 + k2 + q) *
                    disp(-1k1 + k3 - kp + q) *
                    fermidist(-1mu + disp(k3), beta)
                ) +
                u *
                disp(-1k3 + q) *
                disp(-1k1 + k3 - kp + q) *
                fermidist(-1mu + disp(k3), beta) -
                (2 * mu * u * v * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
                2 * mu * u * vp * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
                (u * v * w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
                u * vp * w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
                2 * mu * u * disp(k2) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
                u * w * disp(k2) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
                (2 * mu * u * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
                u * v * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
                (u * vp * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
                (u * w * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
                (u * disp(k2) * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
                u * disp(k3) ^ 2 * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) - (
                    2 *
                    mu *
                    u *
                    disp(-1k - k1 + k2 + q) *
                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                ) - (
                    u *
                    w *
                    disp(-1k - k1 + k2 + q) *
                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                ) +
                u *
                disp(k3) *
                disp(-1k - k1 + k2 + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
                u * v * disp(-1k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
                (u * vp * disp(-1k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
                (
                    u *
                    disp(k2) *
                    disp(-1k3 + q) *
                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                ) +
                u *
                disp(k3) *
                disp(-1k3 + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
                u *
                disp(-1k - k1 + k2 + q) *
                disp(-1k3 + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
                2 *
                mu *
                u *
                disp(-1k1 + k3 - kp + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
                u *
                w *
                disp(-1k1 + k3 - kp + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta) - (
                    u *
                    disp(k3) *
                    disp(-1k1 + k3 - kp + q) *
                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                ) - (
                    u *
                    disp(-1k3 + q) *
                    disp(-1k1 + k3 - kp + q) *
                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                ) - (2 * u * v ^ 2 * fermidist(mu - disp(-1k3 + q), beta)) +
                2 * u * vp ^ 2 * fermidist(mu - disp(-1k3 + q), beta) +
                2 * u * v * w * fermidist(mu - disp(-1k3 + q), beta) -
                (2 * u * vp * w * fermidist(mu - disp(-1k3 + q), beta)) +
                4 * u * v * disp(k2) * fermidist(mu - disp(-1k3 + q), beta) -
                (2 * u * w * disp(k2) * fermidist(mu - disp(-1k3 + q), beta)) -
                (2 * u * disp(k2) ^ 2 * fermidist(mu - disp(-1k3 + q), beta)) -
                (4 * u * v * disp(k3) * fermidist(mu - disp(-1k3 + q), beta)) +
                2 * u * w * disp(k3) * fermidist(mu - disp(-1k3 + q), beta) +
                4 * u * disp(k2) * disp(k3) * fermidist(mu - disp(-1k3 + q), beta) -
                (2 * u * disp(k3) ^ 2 * fermidist(mu - disp(-1k3 + q), beta)) +
                2 * u * v * disp(-1k1 + k3 + kp) * fermidist(mu - disp(-1k3 + q), beta) -
                (2 * u * vp * disp(-1k1 + k3 + kp) * fermidist(mu - disp(-1k3 + q), beta)) -
                (
                    2 *
                    u *
                    disp(k2) *
                    disp(-1k1 + k3 + kp) *
                    fermidist(mu - disp(-1k3 + q), beta)
                ) +
                2 *
                u *
                disp(k3) *
                disp(-1k1 + k3 + kp) *
                fermidist(mu - disp(-1k3 + q), beta) - (
                    4 *
                    u *
                    v *
                    disp(-1k - k1 + k2 + q) *
                    fermidist(mu - disp(-1k3 + q), beta)
                ) +
                2 * u * w * disp(-1k - k1 + k2 + q) * fermidist(mu - disp(-1k3 + q), beta) +
                4 *
                u *
                disp(k2) *
                disp(-1k - k1 + k2 + q) *
                fermidist(mu - disp(-1k3 + q), beta) - (
                    4 *
                    u *
                    disp(k3) *
                    disp(-1k - k1 + k2 + q) *
                    fermidist(mu - disp(-1k3 + q), beta)
                ) +
                2 *
                u *
                disp(-1k1 + k3 + kp) *
                disp(-1k - k1 + k2 + q) *
                fermidist(mu - disp(-1k3 + q), beta) - (
                    2 *
                    u *
                    disp(-1k - k1 + k2 + q) ^ 2 *
                    fermidist(mu - disp(-1k3 + q), beta)
                ) +
                2 *
                u *
                v *
                disp(-1k1 + k3 - kp + q) *
                fermidist(mu - disp(-1k3 + q), beta) +
                2 * u * vp * disp(-1k1 + k3 - kp + q) * fermidist(mu - disp(-1k3 + q), beta) -
                (
                    2 *
                    u *
                    w *
                    disp(-1k1 + k3 - kp + q) *
                    fermidist(mu - disp(-1k3 + q), beta)
                ) - (
                    2 *
                    u *
                    disp(k2) *
                    disp(-1k1 + k3 - kp + q) *
                    fermidist(mu - disp(-1k3 + q), beta)
                ) +
                2 *
                u *
                disp(k3) *
                disp(-1k1 + k3 - kp + q) *
                fermidist(mu - disp(-1k3 + q), beta) - (
                    2 *
                    u *
                    disp(-1k1 + k3 + kp) *
                    disp(-1k1 + k3 - kp + q) *
                    fermidist(mu - disp(-1k3 + q), beta)
                ) +
                2 *
                u *
                disp(-1k - k1 + k2 + q) *
                disp(-1k1 + k3 - kp + q) *
                fermidist(mu - disp(-1k3 + q), beta) -
                (2 * mu * u * v * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) -
                (2 * mu * u * vp * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
                2 * mu * u * w * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
                (u * v * w * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) -
                (u * vp * w * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
                u * w ^ 2 * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) +
                2 * mu * u * disp(k2) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) +
                u * w * disp(k2) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
                (2 * mu * u * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
                u * v * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) +
                u * vp * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
                (2 * u * w * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) -
                (
                    u *
                    disp(k2) *
                    disp(k3) *
                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                ) +
                u * disp(k3) ^ 2 * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) +
                2 *
                mu *
                u *
                disp(-1k1 + k3 + kp) *
                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) +
                u *
                w *
                disp(-1k1 + k3 + kp) *
                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) - (
                    u *
                    disp(k3) *
                    disp(-1k1 + k3 + kp) *
                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                ) - (
                    2 *
                    mu *
                    u *
                    disp(-1k - k1 + k2 + q) *
                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                ) - (
                    u *
                    w *
                    disp(-1k - k1 + k2 + q) *
                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                ) +
                u *
                disp(k3) *
                disp(-1k - k1 + k2 + q) *
                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) +
                u * v * disp(-1k3 + q) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) +
                u * vp * disp(-1k3 + q) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
                (
                    u *
                    w *
                    disp(-1k3 + q) *
                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                ) - (
                    u *
                    disp(k2) *
                    disp(-1k3 + q) *
                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                ) +
                u *
                disp(k3) *
                disp(-1k3 + q) *
                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) - (
                    u *
                    disp(-1k1 + k3 + kp) *
                    disp(-1k3 + q) *
                    fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
                ) +
                u *
                disp(-1k - k1 + k2 + q) *
                disp(-1k3 + q) *
                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
            )
        ) *
        (
            2 *
            (mu - v + w - disp(k1) + disp(k2) - disp(-1k - k1 + k2 + q)) *
            (2v - w - (2 * disp(k2)) + disp(-1k + k1 + k2) + disp(-1k - k1 + k2 + q)) *
            (
                v + vp - w - disp(k2) + disp(k3) - disp(-1k1 + k3 + kp) +
                disp(-1k - k1 + k2 + q)
            ) *
            (mu + v - disp(k2) - disp(-1k1 + q) + disp(-1k - k1 + k2 + q)) *
            (2mu + w - disp(k3) - disp(-1k3 + q)) *
            (
                v - vp - disp(k2) + disp(k3) + disp(-1k - k1 + k2 + q) -
                disp(-1k1 + k3 - kp + q)
            )
        ) ^ -1
    ) +
    (
        u ^ 2 *
        bosedist(disp(k3) - disp(-1k1 + k3 - kp + q), beta) *
        (
            -2 * v ^ 2 + 2 * vp ^ 2 + 2 * v * w - (2 * vp * w) + 4 * v * disp(k2) -
            (2 * w * disp(k2)) - (2 * disp(k2) ^ 2) - (2 * v * disp(-1k + k1 + k2)) +
            2 * vp * disp(-1k + k1 + k2) +
            2 * disp(k2) * disp(-1k + k1 + k2) - (4 * vp * disp(k3)) + 2 * w * disp(k3) -
            (2 * disp(-1k + k1 + k2) * disp(k3)) + 2 * disp(k3) ^ 2 -
            (2 * v * disp(-1k - k1 + k2 + q)) - (2 * vp * disp(-1k - k1 + k2 + q)) +
            2 * w * disp(-1k - k1 + k2 + q) +
            2 * disp(k2) * disp(-1k - k1 + k2 + q) -
            (2 * disp(-1k + k1 + k2) * disp(-1k - k1 + k2 + q)) +
            2 * disp(k3) * disp(-1k - k1 + k2 + q) +
            4 * vp * disp(-1k1 + k3 - kp + q) - (2 * w * disp(-1k1 + k3 - kp + q)) +
            2 * disp(-1k + k1 + k2) * disp(-1k1 + k3 - kp + q) -
            (4 * disp(k3) * disp(-1k1 + k3 - kp + q)) -
            (2 * disp(-1k - k1 + k2 + q) * disp(-1k1 + k3 - kp + q)) +
            2 * disp(-1k1 + k3 - kp + q) ^ 2 -
            (2 * u * v * fermidist(-1mu + disp(k2), beta)) +
            u * w * fermidist(-1mu + disp(k2), beta) +
            2 * u * disp(k2) * fermidist(-1mu + disp(k2), beta) -
            (u * disp(-1k + k1 + k2) * fermidist(-1mu + disp(k2), beta)) -
            (u * disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(k2), beta)) +
            u * v * fermidist(-1mu + disp(-1k + k1 + k2), beta) -
            (u * vp * fermidist(-1mu + disp(-1k + k1 + k2), beta)) -
            (u * disp(k2) * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
            u * disp(k3) * fermidist(-1mu + disp(-1k + k1 + k2), beta) +
            u * disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(-1k + k1 + k2), beta) -
            (u * disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
            u * v * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta) +
            u * vp * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta) -
            (u * w * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) -
            (u * disp(k2) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) +
            u * disp(-1k + k1 + k2) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta) -
            (u * disp(k3) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) +
            u * disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
        ) *
        (
            -4 * mu * u * vp * fermidist(-1mu + disp(k3), beta) +
            2 * mu * u * w * fermidist(-1mu + disp(k3), beta) -
            (2 * u * vp * w * fermidist(-1mu + disp(k3), beta)) +
            u * w ^ 2 * fermidist(-1mu + disp(k3), beta) +
            2 * u * vp * disp(k3) * fermidist(-1mu + disp(k3), beta) -
            (u * w * disp(k3) * fermidist(-1mu + disp(k3), beta)) +
            2 * mu * u * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta) +
            u * w * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta) -
            (u * disp(k3) * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) +
            2 * u * vp * disp(-1k3 + q) * fermidist(-1mu + disp(k3), beta) -
            (u * w * disp(-1k3 + q) * fermidist(-1mu + disp(k3), beta)) -
            (u * disp(-1k1 + k3 + kp) * disp(-1k3 + q) * fermidist(-1mu + disp(k3), beta)) -
            (2 * mu * u * disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta)) -
            (u * w * disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta)) +
            u * disp(k3) * disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta) +
            u *
            disp(-1k3 + q) *
            disp(-1k1 + k3 - kp + q) *
            fermidist(-1mu + disp(k3), beta) +
            4 * mu * u * vp * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
            (2 * mu * u * w * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
            2 * u * vp * w * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
            (u * w ^ 2 * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) -
            (2 * u * vp * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
            u * w * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) - (
                2 *
                mu *
                u *
                disp(-1k1 + k3 + kp) *
                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
            ) - (
                u *
                w *
                disp(-1k1 + k3 + kp) *
                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
            ) +
            u *
            disp(k3) *
            disp(-1k1 + k3 + kp) *
            fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) - (
                2 *
                u *
                vp *
                disp(-1k3 + q) *
                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
            ) +
            u * w * disp(-1k3 + q) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) +
            u *
            disp(-1k1 + k3 + kp) *
            disp(-1k3 + q) *
            fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) +
            2 *
            mu *
            u *
            disp(-1k1 + k3 - kp + q) *
            fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) +
            u *
            w *
            disp(-1k1 + k3 - kp + q) *
            fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) - (
                u *
                disp(k3) *
                disp(-1k1 + k3 - kp + q) *
                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
            ) - (
                u *
                disp(-1k3 + q) *
                disp(-1k1 + k3 - kp + q) *
                fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
            )
        )
    ) *
    (
        2 *
        (2mu + w - disp(k3) - disp(-1k3 + q)) *
        (mu - vp + w - disp(k1) + disp(k3) - disp(-1k1 + k3 - kp + q)) *
        (
            -1v - vp + w + disp(k2) - disp(-1k + k1 + k2) + disp(k3) -
            disp(-1k1 + k3 - kp + q)
        ) *
        (-2vp + w + disp(-1k1 + k3 + kp) - disp(-1k1 + k3 - kp + q)) *
        (
            v - vp - disp(k2) + disp(k3) + disp(-1k - k1 + k2 + q) -
            disp(-1k1 + k3 - kp + q)
        ) *
        (mu + vp - disp(k3) - disp(-1k1 + q) + disp(-1k1 + k3 - kp + q))
    ) ^ -1
end
