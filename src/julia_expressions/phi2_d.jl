function phi2_d(
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
            4 * mu * u * w * fermidist(-1mu + disp(k3), beta) +
            4 * u * vp * w * fermidist(-1mu + disp(k3), beta) +
            2 * u * w ^ 2 * fermidist(-1mu + disp(k3), beta) +
            4 * mu * u * disp(k3) * fermidist(-1mu + disp(k3), beta) +
            4 * u * vp * disp(k3) * fermidist(-1mu + disp(k3), beta) +
            2 * u * w * disp(k3) * fermidist(-1mu + disp(k3), beta) -
            (2 * u * w * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) -
            (2 * u * disp(k3) * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) -
            (4 * mu * u * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) -
            (4 * u * vp * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) -
            (2 * u * w * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) +
            2 * u * disp(-1k1 + k3 + kp) * disp(k3 + q) * fermidist(-1mu + disp(k3), beta) -
            (2 * u * w * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta)) -
            (2 * u * disp(k3) * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta)) +
            2 *
            u *
            disp(k3 + q) *
            disp(k1 - k3 + kp + q) *
            fermidist(-1mu + disp(k3), beta) -
            (4 * mu * u * w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
            (4 * u * vp * w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
            (2 * u * w ^ 2 * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
            (4 * mu * u * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
            (4 * u * vp * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
            (2 * u * w * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
            2 *
            u *
            w *
            disp(-1k1 + k3 + kp) *
            fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            2 *
            u *
            disp(k3) *
            disp(-1k1 + k3 + kp) *
            fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            4 * mu * u * disp(k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            4 * u * vp * disp(k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            2 * u * w * disp(k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) - (
                2 *
                u *
                disp(-1k1 + k3 + kp) *
                disp(k3 + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
            ) +
            2 *
            u *
            w *
            disp(k1 - k3 + kp + q) *
            fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            2 *
            u *
            disp(k3) *
            disp(k1 - k3 + kp + q) *
            fermidist(-1mu + disp(-1k1 + k3 + kp), beta) - (
                2 *
                u *
                disp(k3 + q) *
                disp(k1 - k3 + kp + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
            )
        ) *
        (
            2 * mu * v + v ^ 2 - (2 * mu * vp) - vp ^ 2 + v * w - (vp * w) -
            (2 * mu * disp(k2)) - (2 * v * disp(k2)) - (w * disp(k2)) +
            disp(k2) ^ 2 +
            2 * mu * disp(-1k + k1 + k2) +
            v * disp(-1k + k1 + k2) +
            vp * disp(-1k + k1 + k2) +
            w * disp(-1k + k1 + k2) - (disp(k2) * disp(-1k + k1 + k2)) -
            (2 * mu * disp(k3)) - (2 * vp * disp(k3)) - (w * disp(k3)) +
            disp(-1k + k1 + k2) * disp(k3) - disp(k3) ^ 2 +
            2 * mu * disp(-1k1 + k3 + kp) +
            2 * vp * disp(-1k1 + k3 + kp) +
            w * disp(-1k1 + k3 + kp) - (disp(-1k + k1 + k2) * disp(-1k1 + k3 + kp)) +
            2 * disp(k3) * disp(-1k1 + k3 + kp) - disp(-1k1 + k3 + kp) ^ 2 -
            (v * disp(k + k1 - k2 + q)) +
            vp * disp(k + k1 - k2 + q) +
            disp(k2) * disp(k + k1 - k2 + q) -
            (disp(-1k + k1 + k2) * disp(k + k1 - k2 + q)) +
            disp(k3) * disp(k + k1 - k2 + q) -
            (disp(-1k1 + k3 + kp) * disp(k + k1 - k2 + q)) +
            4 * mu * u * fermidist(-1mu + disp(k2), beta) +
            u * v * fermidist(-1mu + disp(k2), beta) +
            3 * u * vp * fermidist(-1mu + disp(k2), beta) +
            2 * u * w * fermidist(-1mu + disp(k2), beta) -
            (u * disp(k2) * fermidist(-1mu + disp(k2), beta)) -
            (u * disp(-1k + k1 + k2) * fermidist(-1mu + disp(k2), beta)) +
            3 * u * disp(k3) * fermidist(-1mu + disp(k2), beta) -
            (3 * u * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k2), beta)) -
            (2 * u * disp(k + k1 - k2 + q) * fermidist(-1mu + disp(k2), beta)) -
            (4 * mu * u * fermidist(-1mu + disp(-1k + k1 + k2), beta)) -
            (2 * u * v * fermidist(-1mu + disp(-1k + k1 + k2), beta)) -
            (2 * u * vp * fermidist(-1mu + disp(-1k + k1 + k2), beta)) -
            (2 * u * w * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
            2 * u * disp(k2) * fermidist(-1mu + disp(-1k + k1 + k2), beta) -
            (2 * u * disp(k3) * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
            2 * u * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(-1k + k1 + k2), beta) +
            2 * u * disp(k + k1 - k2 + q) * fermidist(-1mu + disp(-1k + k1 + k2), beta) +
            u * v * fermidist(mu - disp(k + k1 - k2 + q), beta) -
            (u * vp * fermidist(mu - disp(k + k1 - k2 + q), beta)) -
            (u * disp(k2) * fermidist(mu - disp(k + k1 - k2 + q), beta)) +
            u * disp(-1k + k1 + k2) * fermidist(mu - disp(k + k1 - k2 + q), beta) -
            (u * disp(k3) * fermidist(mu - disp(k + k1 - k2 + q), beta)) +
            u * disp(-1k1 + k3 + kp) * fermidist(mu - disp(k + k1 - k2 + q), beta)
        )
    ) *
    (
        (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) *
        (-1v + vp + disp(k2) - disp(-1k + k1 + k2) + disp(k3) - disp(-1k1 + k3 + kp)) *
        (mu + vp + w + disp(k3) - disp(-1k1 + k3 + kp) - disp(k1 + q)) *
        (
            2mu + v + vp + w - disp(k2) + disp(k3) - disp(-1k1 + k3 + kp) -
            disp(k + k1 - k2 + q)
        ) *
        (w + disp(k3) - disp(k3 + q)) *
        (2mu + 2vp + w - disp(-1k1 + k3 + kp) - disp(k1 - k3 + kp + q))
    ) ^ -1 +
    (
        u ^ 2 *
        bosedist(-disp(k2) + disp(-1k + k1 + k2), beta) *
        (
            4 * mu * u * fermidist(-1mu + disp(k2), beta) +
            4 * u * v * fermidist(-1mu + disp(k2), beta) +
            2 * u * w * fermidist(-1mu + disp(k2), beta) -
            (4 * u * disp(k2) * fermidist(-1mu + disp(k2), beta)) +
            2 * u * disp(-1k + k1 + k2) * fermidist(-1mu + disp(k2), beta) -
            (2 * u * disp(k + k1 - k2 + q) * fermidist(-1mu + disp(k2), beta)) -
            (4 * mu * u * fermidist(-1mu + disp(-1k + k1 + k2), beta)) -
            (4 * u * v * fermidist(-1mu + disp(-1k + k1 + k2), beta)) -
            (2 * u * w * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
            4 * u * disp(k2) * fermidist(-1mu + disp(-1k + k1 + k2), beta) -
            (2 * u * disp(-1k + k1 + k2) * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
            2 * u * disp(k + k1 - k2 + q) * fermidist(-1mu + disp(-1k + k1 + k2), beta)
        ) *
        (
            2 * mu * v * w + v ^ 2 * w - (2 * mu * vp * w) - (vp ^ 2 * w) + v * w ^ 2 -
            (vp * w ^ 2) - (2 * mu * w * disp(k2)) - (2 * v * w * disp(k2)) -
            (w ^ 2 * disp(k2)) +
            w * disp(k2) ^ 2 +
            2 * mu * w * disp(-1k + k1 + k2) +
            2 * v * w * disp(-1k + k1 + k2) +
            w ^ 2 * disp(-1k + k1 + k2) - (2 * w * disp(k2) * disp(-1k + k1 + k2)) +
            w * disp(-1k + k1 + k2) ^ 2 +
            2 * mu * v * disp(k3) +
            v ^ 2 * disp(k3) - (2 * mu * vp * disp(k3)) - (vp ^ 2 * disp(k3)) -
            (2 * mu * w * disp(k3)) - (v * w * disp(k3)) - (vp * w * disp(k3)) -
            (w ^ 2 * disp(k3)) - (2 * mu * disp(k2) * disp(k3)) -
            (2 * v * disp(k2) * disp(k3)) +
            w * disp(k2) * disp(k3) +
            disp(k2) ^ 2 * disp(k3) +
            2 * mu * disp(-1k + k1 + k2) * disp(k3) +
            2 * v * disp(-1k + k1 + k2) * disp(k3) - (w * disp(-1k + k1 + k2) * disp(k3)) -
            (2 * disp(k2) * disp(-1k + k1 + k2) * disp(k3)) +
            disp(-1k + k1 + k2) ^ 2 * disp(k3) - (2 * mu * disp(k3) ^ 2) -
            (2 * v * disp(k3) ^ 2) + 2 * disp(k2) * disp(k3) ^ 2 -
            (2 * disp(-1k + k1 + k2) * disp(k3) ^ 2) +
            disp(k3) ^ 3 +
            2 * mu * w * disp(-1k1 + k3 + kp) +
            v * w * disp(-1k1 + k3 + kp) +
            vp * w * disp(-1k1 + k3 + kp) +
            w ^ 2 * disp(-1k1 + k3 + kp) - (w * disp(k2) * disp(-1k1 + k3 + kp)) +
            w * disp(-1k + k1 + k2) * disp(-1k1 + k3 + kp) +
            2 * mu * disp(k3) * disp(-1k1 + k3 + kp) +
            v * disp(k3) * disp(-1k1 + k3 + kp) +
            vp * disp(k3) * disp(-1k1 + k3 + kp) -
            (disp(k2) * disp(k3) * disp(-1k1 + k3 + kp)) +
            disp(-1k + k1 + k2) * disp(k3) * disp(-1k1 + k3 + kp) -
            (disp(k3) ^ 2 * disp(-1k1 + k3 + kp)) - (2 * mu * v * disp(k3 + q)) -
            (v ^ 2 * disp(k3 + q)) +
            2 * mu * vp * disp(k3 + q) +
            vp ^ 2 * disp(k3 + q) - (v * w * disp(k3 + q)) +
            vp * w * disp(k3 + q) +
            2 * mu * disp(k2) * disp(k3 + q) +
            2 * v * disp(k2) * disp(k3 + q) +
            w * disp(k2) * disp(k3 + q) - (disp(k2) ^ 2 * disp(k3 + q)) -
            (2 * mu * disp(-1k + k1 + k2) * disp(k3 + q)) -
            (2 * v * disp(-1k + k1 + k2) * disp(k3 + q)) -
            (w * disp(-1k + k1 + k2) * disp(k3 + q)) +
            2 * disp(k2) * disp(-1k + k1 + k2) * disp(k3 + q) -
            (disp(-1k + k1 + k2) ^ 2 * disp(k3 + q)) +
            2 * mu * disp(k3) * disp(k3 + q) +
            2 * v * disp(k3) * disp(k3 + q) +
            w * disp(k3) * disp(k3 + q) - (2 * disp(k2) * disp(k3) * disp(k3 + q)) +
            2 * disp(-1k + k1 + k2) * disp(k3) * disp(k3 + q) -
            (disp(k3) ^ 2 * disp(k3 + q)) - (2 * mu * disp(-1k1 + k3 + kp) * disp(k3 + q)) -
            (v * disp(-1k1 + k3 + kp) * disp(k3 + q)) -
            (vp * disp(-1k1 + k3 + kp) * disp(k3 + q)) -
            (w * disp(-1k1 + k3 + kp) * disp(k3 + q)) +
            disp(k2) * disp(-1k1 + k3 + kp) * disp(k3 + q) -
            (disp(-1k + k1 + k2) * disp(-1k1 + k3 + kp) * disp(k3 + q)) +
            disp(k3) * disp(-1k1 + k3 + kp) * disp(k3 + q) -
            (v * w * disp(k1 - k3 + kp + q)) +
            vp * w * disp(k1 - k3 + kp + q) +
            w * disp(k2) * disp(k1 - k3 + kp + q) -
            (w * disp(-1k + k1 + k2) * disp(k1 - k3 + kp + q)) -
            (v * disp(k3) * disp(k1 - k3 + kp + q)) +
            vp * disp(k3) * disp(k1 - k3 + kp + q) +
            w * disp(k3) * disp(k1 - k3 + kp + q) +
            disp(k2) * disp(k3) * disp(k1 - k3 + kp + q) -
            (disp(-1k + k1 + k2) * disp(k3) * disp(k1 - k3 + kp + q)) +
            disp(k3) ^ 2 * disp(k1 - k3 + kp + q) -
            (w * disp(-1k1 + k3 + kp) * disp(k1 - k3 + kp + q)) -
            (disp(k3) * disp(-1k1 + k3 + kp) * disp(k1 - k3 + kp + q)) +
            v * disp(k3 + q) * disp(k1 - k3 + kp + q) -
            (vp * disp(k3 + q) * disp(k1 - k3 + kp + q)) -
            (disp(k2) * disp(k3 + q) * disp(k1 - k3 + kp + q)) +
            disp(-1k + k1 + k2) * disp(k3 + q) * disp(k1 - k3 + kp + q) -
            (disp(k3) * disp(k3 + q) * disp(k1 - k3 + kp + q)) +
            disp(-1k1 + k3 + kp) * disp(k3 + q) * disp(k1 - k3 + kp + q) +
            2 * mu * u * v * fermidist(-1mu + disp(k3), beta) +
            u * v ^ 2 * fermidist(-1mu + disp(k3), beta) -
            (2 * mu * u * vp * fermidist(-1mu + disp(k3), beta)) -
            (u * vp ^ 2 * fermidist(-1mu + disp(k3), beta)) +
            4 * mu * u * w * fermidist(-1mu + disp(k3), beta) +
            2 * u * v * w * fermidist(-1mu + disp(k3), beta) +
            2 * u * vp * w * fermidist(-1mu + disp(k3), beta) +
            2 * u * w ^ 2 * fermidist(-1mu + disp(k3), beta) -
            (2 * mu * u * disp(k2) * fermidist(-1mu + disp(k3), beta)) -
            (2 * u * v * disp(k2) * fermidist(-1mu + disp(k3), beta)) -
            (2 * u * w * disp(k2) * fermidist(-1mu + disp(k3), beta)) +
            u * disp(k2) ^ 2 * fermidist(-1mu + disp(k3), beta) +
            2 * mu * u * disp(-1k + k1 + k2) * fermidist(-1mu + disp(k3), beta) +
            2 * u * v * disp(-1k + k1 + k2) * fermidist(-1mu + disp(k3), beta) +
            2 * u * w * disp(-1k + k1 + k2) * fermidist(-1mu + disp(k3), beta) -
            (2 * u * disp(k2) * disp(-1k + k1 + k2) * fermidist(-1mu + disp(k3), beta)) +
            u * disp(-1k + k1 + k2) ^ 2 * fermidist(-1mu + disp(k3), beta) +
            2 * mu * u * disp(k3) * fermidist(-1mu + disp(k3), beta) -
            (u * v * disp(k3) * fermidist(-1mu + disp(k3), beta)) +
            3 * u * vp * disp(k3) * fermidist(-1mu + disp(k3), beta) +
            u * disp(k2) * disp(k3) * fermidist(-1mu + disp(k3), beta) -
            (u * disp(-1k + k1 + k2) * disp(k3) * fermidist(-1mu + disp(k3), beta)) +
            2 * mu * u * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta) +
            u * v * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta) +
            u * vp * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta) -
            (u * disp(k2) * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) +
            u *
            disp(-1k + k1 + k2) *
            disp(-1k1 + k3 + kp) *
            fermidist(-1mu + disp(k3), beta) -
            (2 * u * disp(k3) * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) -
            (4 * mu * u * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) -
            (u * v * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) -
            (3 * u * vp * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) -
            (2 * u * w * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) +
            u * disp(k2) * disp(k3 + q) * fermidist(-1mu + disp(k3), beta) -
            (u * disp(-1k + k1 + k2) * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) +
            u * disp(k3) * disp(k3 + q) * fermidist(-1mu + disp(k3), beta) +
            u * disp(-1k1 + k3 + kp) * disp(k3 + q) * fermidist(-1mu + disp(k3), beta) -
            (u * v * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta)) +
            u * vp * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta) -
            (2 * u * w * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta)) +
            u * disp(k2) * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta) - (
                u *
                disp(-1k + k1 + k2) *
                disp(k1 - k3 + kp + q) *
                fermidist(-1mu + disp(k3), beta)
            ) - (u * disp(k3) * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta)) -
            (
                u *
                disp(-1k1 + k3 + kp) *
                disp(k1 - k3 + kp + q) *
                fermidist(-1mu + disp(k3), beta)
            ) +
            2 *
            u *
            disp(k3 + q) *
            disp(k1 - k3 + kp + q) *
            fermidist(-1mu + disp(k3), beta) -
            (4 * mu * u * w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
            (2 * u * v * w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
            (2 * u * vp * w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
            (2 * u * w ^ 2 * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
            2 * u * w * disp(k2) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) - (
                2 *
                u *
                w *
                disp(-1k + k1 + k2) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
            ) - (4 * mu * u * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
            (2 * u * v * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
            (2 * u * vp * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
            2 * u * disp(k2) * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) - (
                2 *
                u *
                disp(-1k + k1 + k2) *
                disp(k3) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
            ) +
            2 * u * disp(k3) ^ 2 * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            4 * mu * u * disp(k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            2 * u * v * disp(k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            2 * u * vp * disp(k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            2 * u * w * disp(k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) - (
                2 *
                u *
                disp(k2) *
                disp(k3 + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
            ) +
            2 *
            u *
            disp(-1k + k1 + k2) *
            disp(k3 + q) *
            fermidist(-1mu + disp(-1k1 + k3 + kp), beta) - (
                2 *
                u *
                disp(k3) *
                disp(k3 + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
            ) +
            2 *
            u *
            w *
            disp(k1 - k3 + kp + q) *
            fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            2 *
            u *
            disp(k3) *
            disp(k1 - k3 + kp + q) *
            fermidist(-1mu + disp(-1k1 + k3 + kp), beta) - (
                2 *
                u *
                disp(k3 + q) *
                disp(k1 - k3 + kp + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
            ) - (2 * mu * u * v * fermidist(-1mu + disp(k3 + q), beta)) -
            (u * v ^ 2 * fermidist(-1mu + disp(k3 + q), beta)) +
            2 * mu * u * vp * fermidist(-1mu + disp(k3 + q), beta) +
            u * vp ^ 2 * fermidist(-1mu + disp(k3 + q), beta) -
            (u * v * w * fermidist(-1mu + disp(k3 + q), beta)) +
            u * vp * w * fermidist(-1mu + disp(k3 + q), beta) +
            2 * mu * u * disp(k2) * fermidist(-1mu + disp(k3 + q), beta) +
            2 * u * v * disp(k2) * fermidist(-1mu + disp(k3 + q), beta) +
            u * w * disp(k2) * fermidist(-1mu + disp(k3 + q), beta) -
            (u * disp(k2) ^ 2 * fermidist(-1mu + disp(k3 + q), beta)) -
            (2 * mu * u * disp(-1k + k1 + k2) * fermidist(-1mu + disp(k3 + q), beta)) -
            (2 * u * v * disp(-1k + k1 + k2) * fermidist(-1mu + disp(k3 + q), beta)) -
            (u * w * disp(-1k + k1 + k2) * fermidist(-1mu + disp(k3 + q), beta)) +
            2 * u * disp(k2) * disp(-1k + k1 + k2) * fermidist(-1mu + disp(k3 + q), beta) -
            (u * disp(-1k + k1 + k2) ^ 2 * fermidist(-1mu + disp(k3 + q), beta)) +
            2 * mu * u * disp(k3) * fermidist(-1mu + disp(k3 + q), beta) +
            2 * u * v * disp(k3) * fermidist(-1mu + disp(k3 + q), beta) +
            u * w * disp(k3) * fermidist(-1mu + disp(k3 + q), beta) -
            (2 * u * disp(k2) * disp(k3) * fermidist(-1mu + disp(k3 + q), beta)) +
            2 * u * disp(-1k + k1 + k2) * disp(k3) * fermidist(-1mu + disp(k3 + q), beta) -
            (u * disp(k3) ^ 2 * fermidist(-1mu + disp(k3 + q), beta)) -
            (2 * mu * u * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3 + q), beta)) -
            (u * v * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3 + q), beta)) -
            (u * vp * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3 + q), beta)) -
            (u * w * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3 + q), beta)) +
            u * disp(k2) * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3 + q), beta) - (
                u *
                disp(-1k + k1 + k2) *
                disp(-1k1 + k3 + kp) *
                fermidist(-1mu + disp(k3 + q), beta)
            ) +
            u * disp(k3) * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3 + q), beta) +
            u * v * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3 + q), beta) -
            (u * vp * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3 + q), beta)) -
            (u * disp(k2) * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3 + q), beta)) +
            u *
            disp(-1k + k1 + k2) *
            disp(k1 - k3 + kp + q) *
            fermidist(-1mu + disp(k3 + q), beta) -
            (u * disp(k3) * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3 + q), beta)) +
            u *
            disp(-1k1 + k3 + kp) *
            disp(k1 - k3 + kp + q) *
            fermidist(-1mu + disp(k3 + q), beta) +
            u * v * w * fermidist(mu - disp(k1 - k3 + kp + q), beta) -
            (u * vp * w * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
            (u * w * disp(k2) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) +
            u * w * disp(-1k + k1 + k2) * fermidist(mu - disp(k1 - k3 + kp + q), beta) +
            u * v * disp(k3) * fermidist(mu - disp(k1 - k3 + kp + q), beta) -
            (u * vp * disp(k3) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
            (u * w * disp(k3) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
            (u * disp(k2) * disp(k3) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) +
            u *
            disp(-1k + k1 + k2) *
            disp(k3) *
            fermidist(mu - disp(k1 - k3 + kp + q), beta) -
            (u * disp(k3) ^ 2 * fermidist(mu - disp(k1 - k3 + kp + q), beta)) +
            u * w * disp(-1k1 + k3 + kp) * fermidist(mu - disp(k1 - k3 + kp + q), beta) +
            u *
            disp(k3) *
            disp(-1k1 + k3 + kp) *
            fermidist(mu - disp(k1 - k3 + kp + q), beta) -
            (u * v * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) +
            u * vp * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta) +
            u * disp(k2) * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta) - (
                u *
                disp(-1k + k1 + k2) *
                disp(k3 + q) *
                fermidist(mu - disp(k1 - k3 + kp + q), beta)
            ) + u * disp(k3) * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta) -
            (
                u *
                disp(-1k1 + k3 + kp) *
                disp(k3 + q) *
                fermidist(mu - disp(k1 - k3 + kp + q), beta)
            )
        )
    ) *
    (
        (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) *
        (v - vp - disp(k2) + disp(-1k + k1 + k2) - disp(k3) + disp(-1k1 + k3 + kp)) *
        (mu + v + w - disp(k2) + disp(-1k + k1 + k2) - disp(k1 + q)) *
        (2mu + 2v + w - (2 * disp(k2)) + disp(-1k + k1 + k2) - disp(k + k1 - k2 + q)) *
        (w + disp(k3) - disp(k3 + q)) *
        (
            2mu + v + vp + w - disp(k2) + disp(-1k + k1 + k2) - disp(k3) -
            disp(k1 - k3 + kp + q)
        )
    ) ^ -1 +
    (
        u ^ 2 *
        fermidist(-1mu + disp(k1), beta) *
        (
            mu ^ 2 + 2 * mu * v + v ^ 2 + mu * w + v * w - (w * disp(k1)) - disp(k1) ^ 2 -
            (2 * mu * disp(k2)) - (2 * v * disp(k2)) - (w * disp(k2)) +
            disp(k2) ^ 2 +
            mu * disp(-1k + k1 + k2) +
            v * disp(-1k + k1 + k2) +
            w * disp(-1k + k1 + k2) +
            disp(k1) * disp(-1k + k1 + k2) - (disp(k2) * disp(-1k + k1 + k2)) -
            (mu * disp(k + k1 - k2 + q)) - (v * disp(k + k1 - k2 + q)) +
            disp(k1) * disp(k + k1 - k2 + q) +
            disp(k2) * disp(k + k1 - k2 + q) -
            (disp(-1k + k1 + k2) * disp(k + k1 - k2 + q)) +
            mu * u * fermidist(-1mu + disp(k2), beta) +
            u * v * fermidist(-1mu + disp(k2), beta) +
            2 * u * w * fermidist(-1mu + disp(k2), beta) +
            3 * u * disp(k1) * fermidist(-1mu + disp(k2), beta) -
            (u * disp(k2) * fermidist(-1mu + disp(k2), beta)) -
            (u * disp(-1k + k1 + k2) * fermidist(-1mu + disp(k2), beta)) -
            (2 * u * disp(k + k1 - k2 + q) * fermidist(-1mu + disp(k2), beta)) -
            (2 * mu * u * fermidist(-1mu + disp(-1k + k1 + k2), beta)) -
            (2 * u * v * fermidist(-1mu + disp(-1k + k1 + k2), beta)) -
            (2 * u * w * fermidist(-1mu + disp(-1k + k1 + k2), beta)) -
            (2 * u * disp(k1) * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
            2 * u * disp(k2) * fermidist(-1mu + disp(-1k + k1 + k2), beta) +
            2 * u * disp(k + k1 - k2 + q) * fermidist(-1mu + disp(-1k + k1 + k2), beta) +
            mu * u * fermidist(mu - disp(k + k1 - k2 + q), beta) +
            u * v * fermidist(mu - disp(k + k1 - k2 + q), beta) -
            (u * disp(k1) * fermidist(mu - disp(k + k1 - k2 + q), beta)) -
            (u * disp(k2) * fermidist(mu - disp(k + k1 - k2 + q), beta)) +
            u * disp(-1k + k1 + k2) * fermidist(mu - disp(k + k1 - k2 + q), beta)
        ) *
        (
            mu ^ 2 * w + 2 * mu * vp * w + vp ^ 2 * w + mu * w ^ 2 + vp * w ^ 2 -
            (w ^ 2 * disp(k1)) - (w * disp(k1) ^ 2) +
            mu ^ 2 * disp(k3) +
            2 * mu * vp * disp(k3) +
            vp ^ 2 * disp(k3) +
            mu * w * disp(k3) +
            vp * w * disp(k3) +
            w ^ 2 * disp(k3) +
            w * disp(k1) * disp(k3) - (disp(k1) ^ 2 * disp(k3)) +
            2 * disp(k1) * disp(k3) ^ 2 - disp(k3) ^ 3 - (mu * w * disp(-1k1 + k3 + kp)) -
            (vp * w * disp(-1k1 + k3 + kp)) - (w ^ 2 * disp(-1k1 + k3 + kp)) -
            (w * disp(k1) * disp(-1k1 + k3 + kp)) - (mu * disp(k3) * disp(-1k1 + k3 + kp)) -
            (vp * disp(k3) * disp(-1k1 + k3 + kp)) -
            (disp(k1) * disp(k3) * disp(-1k1 + k3 + kp)) +
            disp(k3) ^ 2 * disp(-1k1 + k3 + kp) - (mu ^ 2 * disp(k3 + q)) -
            (2 * mu * vp * disp(k3 + q)) - (vp ^ 2 * disp(k3 + q)) -
            (mu * w * disp(k3 + q)) - (vp * w * disp(k3 + q)) +
            w * disp(k1) * disp(k3 + q) +
            disp(k1) ^ 2 * disp(k3 + q) - (w * disp(k3) * disp(k3 + q)) -
            (2 * disp(k1) * disp(k3) * disp(k3 + q)) +
            disp(k3) ^ 2 * disp(k3 + q) +
            mu * disp(-1k1 + k3 + kp) * disp(k3 + q) +
            vp * disp(-1k1 + k3 + kp) * disp(k3 + q) +
            w * disp(-1k1 + k3 + kp) * disp(k3 + q) +
            disp(k1) * disp(-1k1 + k3 + kp) * disp(k3 + q) -
            (disp(k3) * disp(-1k1 + k3 + kp) * disp(k3 + q)) -
            (mu * w * disp(k1 - k3 + kp + q)) - (vp * w * disp(k1 - k3 + kp + q)) +
            w * disp(k1) * disp(k1 - k3 + kp + q) -
            (mu * disp(k3) * disp(k1 - k3 + kp + q)) -
            (vp * disp(k3) * disp(k1 - k3 + kp + q)) -
            (w * disp(k3) * disp(k1 - k3 + kp + q)) +
            disp(k1) * disp(k3) * disp(k1 - k3 + kp + q) -
            (disp(k3) ^ 2 * disp(k1 - k3 + kp + q)) +
            w * disp(-1k1 + k3 + kp) * disp(k1 - k3 + kp + q) +
            disp(k3) * disp(-1k1 + k3 + kp) * disp(k1 - k3 + kp + q) +
            mu * disp(k3 + q) * disp(k1 - k3 + kp + q) +
            vp * disp(k3 + q) * disp(k1 - k3 + kp + q) -
            (disp(k1) * disp(k3 + q) * disp(k1 - k3 + kp + q)) +
            disp(k3) * disp(k3 + q) * disp(k1 - k3 + kp + q) -
            (disp(-1k1 + k3 + kp) * disp(k3 + q) * disp(k1 - k3 + kp + q)) +
            mu ^ 2 * u * fermidist(-1mu + disp(k3), beta) +
            2 * mu * u * vp * fermidist(-1mu + disp(k3), beta) +
            u * vp ^ 2 * fermidist(-1mu + disp(k3), beta) -
            (2 * mu * u * w * fermidist(-1mu + disp(k3), beta)) -
            (2 * u * vp * w * fermidist(-1mu + disp(k3), beta)) -
            (2 * u * w ^ 2 * fermidist(-1mu + disp(k3), beta)) -
            (2 * u * w * disp(k1) * fermidist(-1mu + disp(k3), beta)) -
            (u * disp(k1) ^ 2 * fermidist(-1mu + disp(k3), beta)) -
            (3 * mu * u * disp(k3) * fermidist(-1mu + disp(k3), beta)) -
            (3 * u * vp * disp(k3) * fermidist(-1mu + disp(k3), beta)) +
            u * disp(k1) * disp(k3) * fermidist(-1mu + disp(k3), beta) -
            (mu * u * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) -
            (u * vp * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) -
            (u * disp(k1) * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) +
            2 * u * disp(k3) * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta) +
            3 * mu * u * disp(k3 + q) * fermidist(-1mu + disp(k3), beta) +
            3 * u * vp * disp(k3 + q) * fermidist(-1mu + disp(k3), beta) +
            2 * u * w * disp(k3 + q) * fermidist(-1mu + disp(k3), beta) +
            u * disp(k1) * disp(k3 + q) * fermidist(-1mu + disp(k3), beta) -
            (u * disp(k3) * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) -
            (u * disp(-1k1 + k3 + kp) * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) -
            (mu * u * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta)) -
            (u * vp * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta)) +
            2 * u * w * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta) +
            u * disp(k1) * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta) +
            u * disp(k3) * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta) +
            u *
            disp(-1k1 + k3 + kp) *
            disp(k1 - k3 + kp + q) *
            fermidist(-1mu + disp(k3), beta) - (
                2 *
                u *
                disp(k3 + q) *
                disp(k1 - k3 + kp + q) *
                fermidist(-1mu + disp(k3), beta)
            ) +
            2 * mu * u * w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            2 * u * vp * w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            2 * u * w ^ 2 * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            2 * u * w * disp(k1) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            2 * mu * u * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            2 * u * vp * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            2 * u * disp(k1) * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
            (2 * u * disp(k3) ^ 2 * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
            (2 * mu * u * disp(k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
            (2 * u * vp * disp(k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
            (2 * u * w * disp(k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) - (
                2 *
                u *
                disp(k1) *
                disp(k3 + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
            ) +
            2 * u * disp(k3) * disp(k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
            (
                2 *
                u *
                w *
                disp(k1 - k3 + kp + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
            ) - (
                2 *
                u *
                disp(k3) *
                disp(k1 - k3 + kp + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
            ) +
            2 *
            u *
            disp(k3 + q) *
            disp(k1 - k3 + kp + q) *
            fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
            (mu ^ 2 * u * fermidist(-1mu + disp(k3 + q), beta)) -
            (2 * mu * u * vp * fermidist(-1mu + disp(k3 + q), beta)) -
            (u * vp ^ 2 * fermidist(-1mu + disp(k3 + q), beta)) -
            (mu * u * w * fermidist(-1mu + disp(k3 + q), beta)) -
            (u * vp * w * fermidist(-1mu + disp(k3 + q), beta)) +
            u * w * disp(k1) * fermidist(-1mu + disp(k3 + q), beta) +
            u * disp(k1) ^ 2 * fermidist(-1mu + disp(k3 + q), beta) -
            (u * w * disp(k3) * fermidist(-1mu + disp(k3 + q), beta)) -
            (2 * u * disp(k1) * disp(k3) * fermidist(-1mu + disp(k3 + q), beta)) +
            u * disp(k3) ^ 2 * fermidist(-1mu + disp(k3 + q), beta) +
            mu * u * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3 + q), beta) +
            u * vp * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3 + q), beta) +
            u * w * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3 + q), beta) +
            u * disp(k1) * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3 + q), beta) -
            (u * disp(k3) * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3 + q), beta)) +
            mu * u * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3 + q), beta) +
            u * vp * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3 + q), beta) -
            (u * disp(k1) * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3 + q), beta)) +
            u * disp(k3) * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3 + q), beta) -
            (
                u *
                disp(-1k1 + k3 + kp) *
                disp(k1 - k3 + kp + q) *
                fermidist(-1mu + disp(k3 + q), beta)
            ) +
            mu * u * w * fermidist(mu - disp(k1 - k3 + kp + q), beta) +
            u * vp * w * fermidist(mu - disp(k1 - k3 + kp + q), beta) -
            (u * w * disp(k1) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) +
            mu * u * disp(k3) * fermidist(mu - disp(k1 - k3 + kp + q), beta) +
            u * vp * disp(k3) * fermidist(mu - disp(k1 - k3 + kp + q), beta) +
            u * w * disp(k3) * fermidist(mu - disp(k1 - k3 + kp + q), beta) -
            (u * disp(k1) * disp(k3) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) +
            u * disp(k3) ^ 2 * fermidist(mu - disp(k1 - k3 + kp + q), beta) -
            (u * w * disp(-1k1 + k3 + kp) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
            (
                u *
                disp(k3) *
                disp(-1k1 + k3 + kp) *
                fermidist(mu - disp(k1 - k3 + kp + q), beta)
            ) - (mu * u * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
            (u * vp * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) +
            u * disp(k1) * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta) -
            (u * disp(k3) * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) +
            u *
            disp(-1k1 + k3 + kp) *
            disp(k3 + q) *
            fermidist(mu - disp(k1 - k3 + kp + q), beta)
        )
    ) *
    (
        (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) *
        (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) *
        (w + disp(k1) - disp(k1 + q)) *
        (mu + v + w + disp(k1) - disp(k2) - disp(k + k1 - k2 + q)) *
        (w + disp(k3) - disp(k3 + q)) *
        (mu + vp + w + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q))
    ) ^ -1 - (
        (
            u ^ 2 *
            fermidist(-1mu + disp(k1 + q), beta) *
            (
                -mu ^ 2 - (2 * mu * v) - v ^ 2 - (mu * w) - (v * w) +
                2 * mu * disp(k2) +
                2 * v * disp(k2) +
                w * disp(k2) - disp(k2) ^ 2 - (mu * disp(-1k + k1 + k2)) -
                (v * disp(-1k + k1 + k2)) + disp(k2) * disp(-1k + k1 + k2) -
                (w * disp(k1 + q)) - (disp(-1k + k1 + k2) * disp(k1 + q)) +
                disp(k1 + q) ^ 2 +
                mu * disp(k + k1 - k2 + q) +
                v * disp(k + k1 - k2 + q) +
                w * disp(k + k1 - k2 + q) - (disp(k2) * disp(k + k1 - k2 + q)) +
                disp(-1k + k1 + k2) * disp(k + k1 - k2 + q) -
                (disp(k1 + q) * disp(k + k1 - k2 + q)) -
                (mu * u * fermidist(-1mu + disp(k2), beta)) -
                (u * v * fermidist(-1mu + disp(k2), beta)) +
                u * w * fermidist(-1mu + disp(k2), beta) +
                u * disp(k2) * fermidist(-1mu + disp(k2), beta) +
                u * disp(-1k + k1 + k2) * fermidist(-1mu + disp(k2), beta) -
                (3 * u * disp(k1 + q) * fermidist(-1mu + disp(k2), beta)) +
                2 * u * disp(k + k1 - k2 + q) * fermidist(-1mu + disp(k2), beta) +
                2 * mu * u * fermidist(-1mu + disp(-1k + k1 + k2), beta) +
                2 * u * v * fermidist(-1mu + disp(-1k + k1 + k2), beta) -
                (2 * u * disp(k2) * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
                2 * u * disp(k1 + q) * fermidist(-1mu + disp(-1k + k1 + k2), beta) - (
                    2 *
                    u *
                    disp(k + k1 - k2 + q) *
                    fermidist(-1mu + disp(-1k + k1 + k2), beta)
                ) - (mu * u * fermidist(mu - disp(k + k1 - k2 + q), beta)) -
                (u * v * fermidist(mu - disp(k + k1 - k2 + q), beta)) -
                (u * w * fermidist(mu - disp(k + k1 - k2 + q), beta)) +
                u * disp(k2) * fermidist(mu - disp(k + k1 - k2 + q), beta) -
                (u * disp(-1k + k1 + k2) * fermidist(mu - disp(k + k1 - k2 + q), beta)) +
                u * disp(k1 + q) * fermidist(mu - disp(k + k1 - k2 + q), beta)
            ) *
            (
                -(mu ^ 2 * w) - (2 * mu * vp * w) - (vp ^ 2 * w) - (mu * w ^ 2) -
                (vp * w ^ 2) - (mu ^ 2 * disp(k3)) - (2 * mu * vp * disp(k3)) -
                (vp ^ 2 * disp(k3)) - (mu * w * disp(k3)) - (vp * w * disp(k3)) +
                w ^ 2 * disp(k3) +
                2 * w * disp(k3) ^ 2 +
                disp(k3) ^ 3 +
                mu * w * disp(-1k1 + k3 + kp) +
                vp * w * disp(-1k1 + k3 + kp) +
                mu * disp(k3) * disp(-1k1 + k3 + kp) +
                vp * disp(k3) * disp(-1k1 + k3 + kp) -
                (w * disp(k3) * disp(-1k1 + k3 + kp)) -
                (disp(k3) ^ 2 * disp(-1k1 + k3 + kp)) - (w ^ 2 * disp(k1 + q)) -
                (3 * w * disp(k3) * disp(k1 + q)) - (2 * disp(k3) ^ 2 * disp(k1 + q)) +
                w * disp(-1k1 + k3 + kp) * disp(k1 + q) +
                disp(k3) * disp(-1k1 + k3 + kp) * disp(k1 + q) +
                w * disp(k1 + q) ^ 2 +
                disp(k3) * disp(k1 + q) ^ 2 +
                mu ^ 2 * disp(k3 + q) +
                2 * mu * vp * disp(k3 + q) +
                vp ^ 2 * disp(k3 + q) +
                mu * w * disp(k3 + q) +
                vp * w * disp(k3 + q) - (w * disp(k3) * disp(k3 + q)) -
                (disp(k3) ^ 2 * disp(k3 + q)) - (mu * disp(-1k1 + k3 + kp) * disp(k3 + q)) -
                (vp * disp(-1k1 + k3 + kp) * disp(k3 + q)) +
                disp(k3) * disp(-1k1 + k3 + kp) * disp(k3 + q) +
                w * disp(k1 + q) * disp(k3 + q) +
                2 * disp(k3) * disp(k1 + q) * disp(k3 + q) -
                (disp(-1k1 + k3 + kp) * disp(k1 + q) * disp(k3 + q)) -
                (disp(k1 + q) ^ 2 * disp(k3 + q)) +
                mu * w * disp(k1 - k3 + kp + q) +
                vp * w * disp(k1 - k3 + kp + q) +
                w ^ 2 * disp(k1 - k3 + kp + q) +
                mu * disp(k3) * disp(k1 - k3 + kp + q) +
                vp * disp(k3) * disp(k1 - k3 + kp + q) +
                2 * w * disp(k3) * disp(k1 - k3 + kp + q) +
                disp(k3) ^ 2 * disp(k1 - k3 + kp + q) -
                (w * disp(-1k1 + k3 + kp) * disp(k1 - k3 + kp + q)) -
                (disp(k3) * disp(-1k1 + k3 + kp) * disp(k1 - k3 + kp + q)) -
                (w * disp(k1 + q) * disp(k1 - k3 + kp + q)) -
                (disp(k3) * disp(k1 + q) * disp(k1 - k3 + kp + q)) -
                (mu * disp(k3 + q) * disp(k1 - k3 + kp + q)) -
                (vp * disp(k3 + q) * disp(k1 - k3 + kp + q)) -
                (w * disp(k3 + q) * disp(k1 - k3 + kp + q)) -
                (disp(k3) * disp(k3 + q) * disp(k1 - k3 + kp + q)) +
                disp(-1k1 + k3 + kp) * disp(k3 + q) * disp(k1 - k3 + kp + q) +
                disp(k1 + q) * disp(k3 + q) * disp(k1 - k3 + kp + q) -
                (mu ^ 2 * u * fermidist(-1mu + disp(k3), beta)) -
                (2 * mu * u * vp * fermidist(-1mu + disp(k3), beta)) -
                (u * vp ^ 2 * fermidist(-1mu + disp(k3), beta)) +
                2 * mu * u * w * fermidist(-1mu + disp(k3), beta) +
                2 * u * vp * w * fermidist(-1mu + disp(k3), beta) +
                u * w ^ 2 * fermidist(-1mu + disp(k3), beta) +
                3 * mu * u * disp(k3) * fermidist(-1mu + disp(k3), beta) +
                3 * u * vp * disp(k3) * fermidist(-1mu + disp(k3), beta) +
                u * w * disp(k3) * fermidist(-1mu + disp(k3), beta) +
                mu * u * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta) +
                u * vp * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta) -
                (u * w * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) - (
                    2 *
                    u *
                    disp(k3) *
                    disp(-1k1 + k3 + kp) *
                    fermidist(-1mu + disp(k3), beta)
                ) - (u * disp(k3) * disp(k1 + q) * fermidist(-1mu + disp(k3), beta)) +
                u * disp(-1k1 + k3 + kp) * disp(k1 + q) * fermidist(-1mu + disp(k3), beta) +
                u * disp(k1 + q) ^ 2 * fermidist(-1mu + disp(k3), beta) -
                (3 * mu * u * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) -
                (3 * u * vp * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) -
                (u * w * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) +
                u * disp(k3) * disp(k3 + q) * fermidist(-1mu + disp(k3), beta) +
                u * disp(-1k1 + k3 + kp) * disp(k3 + q) * fermidist(-1mu + disp(k3), beta) -
                (u * disp(k1 + q) * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) +
                mu * u * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta) +
                u * vp * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta) -
                (u * w * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta)) -
                (u * disp(k3) * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta)) -
                (
                    u *
                    disp(-1k1 + k3 + kp) *
                    disp(k1 - k3 + kp + q) *
                    fermidist(-1mu + disp(k3), beta)
                ) - (
                    u *
                    disp(k1 + q) *
                    disp(k1 - k3 + kp + q) *
                    fermidist(-1mu + disp(k3), beta)
                ) +
                2 *
                u *
                disp(k3 + q) *
                disp(k1 - k3 + kp + q) *
                fermidist(-1mu + disp(k3), beta) -
                (2 * mu * u * w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
                (2 * u * vp * w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
                (2 * mu * u * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
                (2 * u * vp * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
                2 * u * w * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
                2 * u * disp(k3) ^ 2 * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
                (2 * u * w * disp(k1 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
                (
                    2 *
                    u *
                    disp(k3) *
                    disp(k1 + q) *
                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                ) +
                2 * mu * u * disp(k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
                2 * u * vp * disp(k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
                (
                    2 *
                    u *
                    disp(k3) *
                    disp(k3 + q) *
                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                ) +
                2 *
                u *
                disp(k1 + q) *
                disp(k3 + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
                2 *
                u *
                w *
                disp(k1 - k3 + kp + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
                2 *
                u *
                disp(k3) *
                disp(k1 - k3 + kp + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta) - (
                    2 *
                    u *
                    disp(k3 + q) *
                    disp(k1 - k3 + kp + q) *
                    fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
                ) +
                mu ^ 2 * u * fermidist(-1mu + disp(k3 + q), beta) +
                2 * mu * u * vp * fermidist(-1mu + disp(k3 + q), beta) +
                u * vp ^ 2 * fermidist(-1mu + disp(k3 + q), beta) +
                mu * u * w * fermidist(-1mu + disp(k3 + q), beta) +
                u * vp * w * fermidist(-1mu + disp(k3 + q), beta) -
                (u * w * disp(k3) * fermidist(-1mu + disp(k3 + q), beta)) -
                (u * disp(k3) ^ 2 * fermidist(-1mu + disp(k3 + q), beta)) -
                (mu * u * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3 + q), beta)) -
                (u * vp * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3 + q), beta)) +
                u * disp(k3) * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3 + q), beta) +
                u * w * disp(k1 + q) * fermidist(-1mu + disp(k3 + q), beta) +
                2 * u * disp(k3) * disp(k1 + q) * fermidist(-1mu + disp(k3 + q), beta) - (
                    u *
                    disp(-1k1 + k3 + kp) *
                    disp(k1 + q) *
                    fermidist(-1mu + disp(k3 + q), beta)
                ) - (u * disp(k1 + q) ^ 2 * fermidist(-1mu + disp(k3 + q), beta)) -
                (mu * u * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3 + q), beta)) -
                (u * vp * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3 + q), beta)) -
                (u * w * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3 + q), beta)) -
                (
                    u *
                    disp(k3) *
                    disp(k1 - k3 + kp + q) *
                    fermidist(-1mu + disp(k3 + q), beta)
                ) +
                u *
                disp(-1k1 + k3 + kp) *
                disp(k1 - k3 + kp + q) *
                fermidist(-1mu + disp(k3 + q), beta) +
                u *
                disp(k1 + q) *
                disp(k1 - k3 + kp + q) *
                fermidist(-1mu + disp(k3 + q), beta) -
                (mu * u * w * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
                (u * vp * w * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
                (u * w ^ 2 * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
                (mu * u * disp(k3) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
                (u * vp * disp(k3) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
                (2 * u * w * disp(k3) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
                (u * disp(k3) ^ 2 * fermidist(mu - disp(k1 - k3 + kp + q), beta)) +
                u *
                w *
                disp(-1k1 + k3 + kp) *
                fermidist(mu - disp(k1 - k3 + kp + q), beta) +
                u *
                disp(k3) *
                disp(-1k1 + k3 + kp) *
                fermidist(mu - disp(k1 - k3 + kp + q), beta) +
                u * w * disp(k1 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta) +
                u * disp(k3) * disp(k1 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta) +
                mu * u * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta) +
                u * vp * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta) +
                u * w * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta) +
                u * disp(k3) * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta) -
                (
                    u *
                    disp(-1k1 + k3 + kp) *
                    disp(k3 + q) *
                    fermidist(mu - disp(k1 - k3 + kp + q), beta)
                ) - (
                    u *
                    disp(k1 + q) *
                    disp(k3 + q) *
                    fermidist(mu - disp(k1 - k3 + kp + q), beta)
                )
            )
        ) *
        (
            (w + disp(k1) - disp(k1 + q)) *
            (mu + v + w - disp(k2) + disp(-1k + k1 + k2) - disp(k1 + q)) *
            (mu + vp + w + disp(k3) - disp(-1k1 + k3 + kp) - disp(k1 + q)) *
            (mu + v - disp(k2) + disp(k1 + q) - disp(k + k1 - k2 + q)) *
            (w + disp(k3) - disp(k3 + q)) *
            (mu + vp - disp(k3) + disp(k1 + q) - disp(k1 - k3 + kp + q))
        ) ^ -1
    ) +
    (
        u ^ 2 *
        bosedist(-2mu + disp(k2) + disp(k + k1 - k2 + q), beta) *
        (
            2 * mu * u * fermidist(-1mu + disp(k2), beta) +
            2 * u * v * fermidist(-1mu + disp(k2), beta) +
            u * w * fermidist(-1mu + disp(k2), beta) -
            (2 * u * disp(k2) * fermidist(-1mu + disp(k2), beta)) +
            u * disp(-1k + k1 + k2) * fermidist(-1mu + disp(k2), beta) -
            (u * disp(k + k1 - k2 + q) * fermidist(-1mu + disp(k2), beta)) -
            (2 * mu * u * fermidist(mu - disp(k + k1 - k2 + q), beta)) -
            (2 * u * v * fermidist(mu - disp(k + k1 - k2 + q), beta)) -
            (u * w * fermidist(mu - disp(k + k1 - k2 + q), beta)) +
            2 * u * disp(k2) * fermidist(mu - disp(k + k1 - k2 + q), beta) -
            (u * disp(-1k + k1 + k2) * fermidist(mu - disp(k + k1 - k2 + q), beta)) +
            u * disp(k + k1 - k2 + q) * fermidist(mu - disp(k + k1 - k2 + q), beta)
        ) *
        (
            2 * mu * v * w + v ^ 2 * w - (2 * mu * vp * w) - (vp ^ 2 * w) + v * w ^ 2 -
            (vp * w ^ 2) - (2 * mu * w * disp(k2)) - (2 * v * w * disp(k2)) -
            (w ^ 2 * disp(k2)) +
            w * disp(k2) ^ 2 +
            2 * mu * v * disp(k3) +
            v ^ 2 * disp(k3) - (2 * mu * vp * disp(k3)) - (vp ^ 2 * disp(k3)) +
            2 * mu * w * disp(k3) +
            3 * v * w * disp(k3) - (vp * w * disp(k3)) + w ^ 2 * disp(k3) -
            (2 * mu * disp(k2) * disp(k3)) - (2 * v * disp(k2) * disp(k3)) -
            (3 * w * disp(k2) * disp(k3)) +
            disp(k2) ^ 2 * disp(k3) +
            2 * mu * disp(k3) ^ 2 +
            2 * v * disp(k3) ^ 2 +
            2 * w * disp(k3) ^ 2 - (2 * disp(k2) * disp(k3) ^ 2) + disp(k3) ^ 3 -
            (v * w * disp(-1k1 + k3 + kp)) +
            vp * w * disp(-1k1 + k3 + kp) +
            w * disp(k2) * disp(-1k1 + k3 + kp) - (v * disp(k3) * disp(-1k1 + k3 + kp)) +
            vp * disp(k3) * disp(-1k1 + k3 + kp) - (w * disp(k3) * disp(-1k1 + k3 + kp)) +
            disp(k2) * disp(k3) * disp(-1k1 + k3 + kp) -
            (disp(k3) ^ 2 * disp(-1k1 + k3 + kp)) - (2 * mu * w * disp(k + k1 - k2 + q)) -
            (2 * v * w * disp(k + k1 - k2 + q)) - (w ^ 2 * disp(k + k1 - k2 + q)) +
            2 * w * disp(k2) * disp(k + k1 - k2 + q) -
            (2 * mu * disp(k3) * disp(k + k1 - k2 + q)) -
            (2 * v * disp(k3) * disp(k + k1 - k2 + q)) -
            (3 * w * disp(k3) * disp(k + k1 - k2 + q)) +
            2 * disp(k2) * disp(k3) * disp(k + k1 - k2 + q) -
            (2 * disp(k3) ^ 2 * disp(k + k1 - k2 + q)) +
            w * disp(-1k1 + k3 + kp) * disp(k + k1 - k2 + q) +
            disp(k3) * disp(-1k1 + k3 + kp) * disp(k + k1 - k2 + q) +
            w * disp(k + k1 - k2 + q) ^ 2 +
            disp(k3) * disp(k + k1 - k2 + q) ^ 2 - (2 * mu * v * disp(k3 + q)) -
            (v ^ 2 * disp(k3 + q)) +
            2 * mu * vp * disp(k3 + q) +
            vp ^ 2 * disp(k3 + q) - (v * w * disp(k3 + q)) +
            vp * w * disp(k3 + q) +
            2 * mu * disp(k2) * disp(k3 + q) +
            2 * v * disp(k2) * disp(k3 + q) +
            w * disp(k2) * disp(k3 + q) - (disp(k2) ^ 2 * disp(k3 + q)) -
            (2 * mu * disp(k3) * disp(k3 + q)) - (2 * v * disp(k3) * disp(k3 + q)) -
            (w * disp(k3) * disp(k3 + q)) + 2 * disp(k2) * disp(k3) * disp(k3 + q) -
            (disp(k3) ^ 2 * disp(k3 + q)) + v * disp(-1k1 + k3 + kp) * disp(k3 + q) -
            (vp * disp(-1k1 + k3 + kp) * disp(k3 + q)) -
            (disp(k2) * disp(-1k1 + k3 + kp) * disp(k3 + q)) +
            disp(k3) * disp(-1k1 + k3 + kp) * disp(k3 + q) +
            2 * mu * disp(k + k1 - k2 + q) * disp(k3 + q) +
            2 * v * disp(k + k1 - k2 + q) * disp(k3 + q) +
            w * disp(k + k1 - k2 + q) * disp(k3 + q) -
            (2 * disp(k2) * disp(k + k1 - k2 + q) * disp(k3 + q)) +
            2 * disp(k3) * disp(k + k1 - k2 + q) * disp(k3 + q) -
            (disp(-1k1 + k3 + kp) * disp(k + k1 - k2 + q) * disp(k3 + q)) -
            (disp(k + k1 - k2 + q) ^ 2 * disp(k3 + q)) +
            2 * mu * w * disp(k1 - k3 + kp + q) +
            v * w * disp(k1 - k3 + kp + q) +
            vp * w * disp(k1 - k3 + kp + q) +
            w ^ 2 * disp(k1 - k3 + kp + q) - (w * disp(k2) * disp(k1 - k3 + kp + q)) +
            2 * mu * disp(k3) * disp(k1 - k3 + kp + q) +
            v * disp(k3) * disp(k1 - k3 + kp + q) +
            vp * disp(k3) * disp(k1 - k3 + kp + q) +
            2 * w * disp(k3) * disp(k1 - k3 + kp + q) -
            (disp(k2) * disp(k3) * disp(k1 - k3 + kp + q)) +
            disp(k3) ^ 2 * disp(k1 - k3 + kp + q) -
            (w * disp(-1k1 + k3 + kp) * disp(k1 - k3 + kp + q)) -
            (disp(k3) * disp(-1k1 + k3 + kp) * disp(k1 - k3 + kp + q)) -
            (w * disp(k + k1 - k2 + q) * disp(k1 - k3 + kp + q)) -
            (disp(k3) * disp(k + k1 - k2 + q) * disp(k1 - k3 + kp + q)) -
            (2 * mu * disp(k3 + q) * disp(k1 - k3 + kp + q)) -
            (v * disp(k3 + q) * disp(k1 - k3 + kp + q)) -
            (vp * disp(k3 + q) * disp(k1 - k3 + kp + q)) -
            (w * disp(k3 + q) * disp(k1 - k3 + kp + q)) +
            disp(k2) * disp(k3 + q) * disp(k1 - k3 + kp + q) -
            (disp(k3) * disp(k3 + q) * disp(k1 - k3 + kp + q)) +
            disp(-1k1 + k3 + kp) * disp(k3 + q) * disp(k1 - k3 + kp + q) +
            disp(k + k1 - k2 + q) * disp(k3 + q) * disp(k1 - k3 + kp + q) +
            2 * mu * u * v * fermidist(-1mu + disp(k3), beta) +
            u * v ^ 2 * fermidist(-1mu + disp(k3), beta) -
            (2 * mu * u * vp * fermidist(-1mu + disp(k3), beta)) -
            (u * vp ^ 2 * fermidist(-1mu + disp(k3), beta)) +
            2 * mu * u * w * fermidist(-1mu + disp(k3), beta) +
            2 * u * vp * w * fermidist(-1mu + disp(k3), beta) +
            u * w ^ 2 * fermidist(-1mu + disp(k3), beta) -
            (2 * mu * u * disp(k2) * fermidist(-1mu + disp(k3), beta)) -
            (2 * u * v * disp(k2) * fermidist(-1mu + disp(k3), beta)) +
            u * disp(k2) ^ 2 * fermidist(-1mu + disp(k3), beta) +
            4 * mu * u * disp(k3) * fermidist(-1mu + disp(k3), beta) +
            u * v * disp(k3) * fermidist(-1mu + disp(k3), beta) +
            3 * u * vp * disp(k3) * fermidist(-1mu + disp(k3), beta) +
            u * w * disp(k3) * fermidist(-1mu + disp(k3), beta) -
            (u * disp(k2) * disp(k3) * fermidist(-1mu + disp(k3), beta)) -
            (u * v * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) +
            u * vp * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta) -
            (u * w * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) +
            u * disp(k2) * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta) -
            (2 * u * disp(k3) * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) -
            (2 * mu * u * disp(k + k1 - k2 + q) * fermidist(-1mu + disp(k3), beta)) -
            (2 * u * v * disp(k + k1 - k2 + q) * fermidist(-1mu + disp(k3), beta)) +
            2 * u * disp(k2) * disp(k + k1 - k2 + q) * fermidist(-1mu + disp(k3), beta) -
            (u * disp(k3) * disp(k + k1 - k2 + q) * fermidist(-1mu + disp(k3), beta)) +
            u *
            disp(-1k1 + k3 + kp) *
            disp(k + k1 - k2 + q) *
            fermidist(-1mu + disp(k3), beta) +
            u * disp(k + k1 - k2 + q) ^ 2 * fermidist(-1mu + disp(k3), beta) -
            (2 * mu * u * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) +
            u * v * disp(k3 + q) * fermidist(-1mu + disp(k3), beta) -
            (3 * u * vp * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) -
            (u * w * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) -
            (u * disp(k2) * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) +
            u * disp(k3) * disp(k3 + q) * fermidist(-1mu + disp(k3), beta) +
            u * disp(-1k1 + k3 + kp) * disp(k3 + q) * fermidist(-1mu + disp(k3), beta) -
            (u * disp(k + k1 - k2 + q) * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) +
            2 * mu * u * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta) +
            u * v * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta) +
            u * vp * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta) -
            (u * w * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta)) -
            (u * disp(k2) * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta)) -
            (u * disp(k3) * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta)) - (
                u *
                disp(-1k1 + k3 + kp) *
                disp(k1 - k3 + kp + q) *
                fermidist(-1mu + disp(k3), beta)
            ) - (
                u *
                disp(k + k1 - k2 + q) *
                disp(k1 - k3 + kp + q) *
                fermidist(-1mu + disp(k3), beta)
            ) +
            2 *
            u *
            disp(k3 + q) *
            disp(k1 - k3 + kp + q) *
            fermidist(-1mu + disp(k3), beta) +
            2 * u * v * w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
            (2 * u * vp * w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
            (2 * u * w * disp(k2) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
            2 * u * v * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
            (2 * u * vp * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
            2 * u * w * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
            (2 * u * disp(k2) * disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
            2 * u * disp(k3) ^ 2 * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) - (
                2 *
                u *
                w *
                disp(k + k1 - k2 + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
            ) - (
                2 *
                u *
                disp(k3) *
                disp(k + k1 - k2 + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
            ) - (2 * u * v * disp(k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
            2 * u * vp * disp(k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            2 * u * disp(k2) * disp(k3 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
            (
                2 *
                u *
                disp(k3) *
                disp(k3 + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
            ) +
            2 *
            u *
            disp(k + k1 - k2 + q) *
            disp(k3 + q) *
            fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            2 *
            u *
            w *
            disp(k1 - k3 + kp + q) *
            fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            2 *
            u *
            disp(k3) *
            disp(k1 - k3 + kp + q) *
            fermidist(-1mu + disp(-1k1 + k3 + kp), beta) - (
                2 *
                u *
                disp(k3 + q) *
                disp(k1 - k3 + kp + q) *
                fermidist(-1mu + disp(-1k1 + k3 + kp), beta)
            ) - (2 * mu * u * v * fermidist(-1mu + disp(k3 + q), beta)) -
            (u * v ^ 2 * fermidist(-1mu + disp(k3 + q), beta)) +
            2 * mu * u * vp * fermidist(-1mu + disp(k3 + q), beta) +
            u * vp ^ 2 * fermidist(-1mu + disp(k3 + q), beta) -
            (u * v * w * fermidist(-1mu + disp(k3 + q), beta)) +
            u * vp * w * fermidist(-1mu + disp(k3 + q), beta) +
            2 * mu * u * disp(k2) * fermidist(-1mu + disp(k3 + q), beta) +
            2 * u * v * disp(k2) * fermidist(-1mu + disp(k3 + q), beta) +
            u * w * disp(k2) * fermidist(-1mu + disp(k3 + q), beta) -
            (u * disp(k2) ^ 2 * fermidist(-1mu + disp(k3 + q), beta)) -
            (2 * mu * u * disp(k3) * fermidist(-1mu + disp(k3 + q), beta)) -
            (2 * u * v * disp(k3) * fermidist(-1mu + disp(k3 + q), beta)) -
            (u * w * disp(k3) * fermidist(-1mu + disp(k3 + q), beta)) +
            2 * u * disp(k2) * disp(k3) * fermidist(-1mu + disp(k3 + q), beta) -
            (u * disp(k3) ^ 2 * fermidist(-1mu + disp(k3 + q), beta)) +
            u * v * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3 + q), beta) -
            (u * vp * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3 + q), beta)) -
            (u * disp(k2) * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3 + q), beta)) +
            u * disp(k3) * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3 + q), beta) +
            2 * mu * u * disp(k + k1 - k2 + q) * fermidist(-1mu + disp(k3 + q), beta) +
            2 * u * v * disp(k + k1 - k2 + q) * fermidist(-1mu + disp(k3 + q), beta) +
            u * w * disp(k + k1 - k2 + q) * fermidist(-1mu + disp(k3 + q), beta) - (
                2 *
                u *
                disp(k2) *
                disp(k + k1 - k2 + q) *
                fermidist(-1mu + disp(k3 + q), beta)
            ) +
            2 *
            u *
            disp(k3) *
            disp(k + k1 - k2 + q) *
            fermidist(-1mu + disp(k3 + q), beta) - (
                u *
                disp(-1k1 + k3 + kp) *
                disp(k + k1 - k2 + q) *
                fermidist(-1mu + disp(k3 + q), beta)
            ) - (u * disp(k + k1 - k2 + q) ^ 2 * fermidist(-1mu + disp(k3 + q), beta)) -
            (2 * mu * u * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3 + q), beta)) -
            (u * v * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3 + q), beta)) -
            (u * vp * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3 + q), beta)) -
            (u * w * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3 + q), beta)) +
            u * disp(k2) * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3 + q), beta) -
            (u * disp(k3) * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3 + q), beta)) +
            u *
            disp(-1k1 + k3 + kp) *
            disp(k1 - k3 + kp + q) *
            fermidist(-1mu + disp(k3 + q), beta) +
            u *
            disp(k + k1 - k2 + q) *
            disp(k1 - k3 + kp + q) *
            fermidist(-1mu + disp(k3 + q), beta) -
            (2 * mu * u * w * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
            (u * v * w * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
            (u * vp * w * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
            (u * w ^ 2 * fermidist(mu - disp(k1 - k3 + kp + q), beta)) +
            u * w * disp(k2) * fermidist(mu - disp(k1 - k3 + kp + q), beta) -
            (2 * mu * u * disp(k3) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
            (u * v * disp(k3) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
            (u * vp * disp(k3) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
            (2 * u * w * disp(k3) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) +
            u * disp(k2) * disp(k3) * fermidist(mu - disp(k1 - k3 + kp + q), beta) -
            (u * disp(k3) ^ 2 * fermidist(mu - disp(k1 - k3 + kp + q), beta)) +
            u * w * disp(-1k1 + k3 + kp) * fermidist(mu - disp(k1 - k3 + kp + q), beta) +
            u *
            disp(k3) *
            disp(-1k1 + k3 + kp) *
            fermidist(mu - disp(k1 - k3 + kp + q), beta) +
            u * w * disp(k + k1 - k2 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta) +
            u *
            disp(k3) *
            disp(k + k1 - k2 + q) *
            fermidist(mu - disp(k1 - k3 + kp + q), beta) +
            2 * mu * u * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta) +
            u * v * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta) +
            u * vp * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta) +
            u * w * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta) -
            (u * disp(k2) * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) +
            u * disp(k3) * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta) - (
                u *
                disp(-1k1 + k3 + kp) *
                disp(k3 + q) *
                fermidist(mu - disp(k1 - k3 + kp + q), beta)
            ) - (
                u *
                disp(k + k1 - k2 + q) *
                disp(k3 + q) *
                fermidist(mu - disp(k1 - k3 + kp + q), beta)
            )
        )
    ) *
    (
        (mu + v + w + disp(k1) - disp(k2) - disp(k + k1 - k2 + q)) *
        (2mu + 2v + w - (2 * disp(k2)) + disp(-1k + k1 + k2) - disp(k + k1 - k2 + q)) *
        (
            2mu + v + vp + w - disp(k2) + disp(k3) - disp(-1k1 + k3 + kp) -
            disp(k + k1 - k2 + q)
        ) *
        (mu + v - disp(k2) + disp(k1 + q) - disp(k + k1 - k2 + q)) *
        (w + disp(k3) - disp(k3 + q)) *
        (v - vp - disp(k2) + disp(k3) - disp(k + k1 - k2 + q) + disp(k1 - k3 + kp + q))
    ) ^ -1 +
    (
        u ^ 2 *
        bosedist(-2mu + disp(k3) + disp(k1 - k3 + kp + q), beta) *
        (
            -2 * mu * v - v ^ 2 + 2 * mu * vp + vp ^ 2 - (v * w) +
            vp * w +
            2 * mu * disp(k2) +
            2 * v * disp(k2) +
            w * disp(k2) - disp(k2) ^ 2 - (v * disp(-1k + k1 + k2)) +
            vp * disp(-1k + k1 + k2) +
            disp(k2) * disp(-1k + k1 + k2) - (2 * mu * disp(k3)) - (2 * vp * disp(k3)) -
            (w * disp(k3)) - (disp(-1k + k1 + k2) * disp(k3)) +
            disp(k3) ^ 2 +
            2 * mu * disp(k + k1 - k2 + q) +
            v * disp(k + k1 - k2 + q) +
            vp * disp(k + k1 - k2 + q) +
            w * disp(k + k1 - k2 + q) - (disp(k2) * disp(k + k1 - k2 + q)) +
            disp(-1k + k1 + k2) * disp(k + k1 - k2 + q) -
            (disp(k3) * disp(k + k1 - k2 + q)) - (2 * mu * disp(k1 - k3 + kp + q)) -
            (2 * vp * disp(k1 - k3 + kp + q)) - (w * disp(k1 - k3 + kp + q)) -
            (disp(-1k + k1 + k2) * disp(k1 - k3 + kp + q)) +
            2 * disp(k3) * disp(k1 - k3 + kp + q) -
            (disp(k + k1 - k2 + q) * disp(k1 - k3 + kp + q)) +
            disp(k1 - k3 + kp + q) ^ 2 +
            2 * mu * u * fermidist(-1mu + disp(k2), beta) -
            (u * v * fermidist(-1mu + disp(k2), beta)) +
            3 * u * vp * fermidist(-1mu + disp(k2), beta) +
            u * w * fermidist(-1mu + disp(k2), beta) +
            u * disp(k2) * fermidist(-1mu + disp(k2), beta) +
            u * disp(-1k + k1 + k2) * fermidist(-1mu + disp(k2), beta) -
            (3 * u * disp(k3) * fermidist(-1mu + disp(k2), beta)) +
            2 * u * disp(k + k1 - k2 + q) * fermidist(-1mu + disp(k2), beta) -
            (3 * u * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k2), beta)) +
            2 * u * v * fermidist(-1mu + disp(-1k + k1 + k2), beta) -
            (2 * u * vp * fermidist(-1mu + disp(-1k + k1 + k2), beta)) -
            (2 * u * disp(k2) * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
            2 * u * disp(k3) * fermidist(-1mu + disp(-1k + k1 + k2), beta) -
            (2 * u * disp(k + k1 - k2 + q) * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
            2 * u * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(-1k + k1 + k2), beta) -
            (2 * mu * u * fermidist(mu - disp(k + k1 - k2 + q), beta)) -
            (u * v * fermidist(mu - disp(k + k1 - k2 + q), beta)) -
            (u * vp * fermidist(mu - disp(k + k1 - k2 + q), beta)) -
            (u * w * fermidist(mu - disp(k + k1 - k2 + q), beta)) +
            u * disp(k2) * fermidist(mu - disp(k + k1 - k2 + q), beta) -
            (u * disp(-1k + k1 + k2) * fermidist(mu - disp(k + k1 - k2 + q), beta)) +
            u * disp(k3) * fermidist(mu - disp(k + k1 - k2 + q), beta) +
            u * disp(k1 - k3 + kp + q) * fermidist(mu - disp(k + k1 - k2 + q), beta)
        ) *
        (
            2 * mu * u * w * fermidist(-1mu + disp(k3), beta) +
            2 * u * vp * w * fermidist(-1mu + disp(k3), beta) +
            u * w ^ 2 * fermidist(-1mu + disp(k3), beta) +
            2 * mu * u * disp(k3) * fermidist(-1mu + disp(k3), beta) +
            2 * u * vp * disp(k3) * fermidist(-1mu + disp(k3), beta) +
            u * w * disp(k3) * fermidist(-1mu + disp(k3), beta) -
            (u * w * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) -
            (u * disp(k3) * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) -
            (2 * mu * u * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) -
            (2 * u * vp * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) -
            (u * w * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) +
            u * disp(-1k1 + k3 + kp) * disp(k3 + q) * fermidist(-1mu + disp(k3), beta) -
            (u * w * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta)) -
            (u * disp(k3) * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta)) +
            u * disp(k3 + q) * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta) -
            (2 * mu * u * w * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
            (2 * u * vp * w * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
            (u * w ^ 2 * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
            (2 * mu * u * disp(k3) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
            (2 * u * vp * disp(k3) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
            (u * w * disp(k3) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) +
            u * w * disp(-1k1 + k3 + kp) * fermidist(mu - disp(k1 - k3 + kp + q), beta) +
            u *
            disp(k3) *
            disp(-1k1 + k3 + kp) *
            fermidist(mu - disp(k1 - k3 + kp + q), beta) +
            2 * mu * u * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta) +
            2 * u * vp * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta) +
            u * w * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta) - (
                u *
                disp(-1k1 + k3 + kp) *
                disp(k3 + q) *
                fermidist(mu - disp(k1 - k3 + kp + q), beta)
            ) +
            u * w * disp(k1 - k3 + kp + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta) +
            u *
            disp(k3) *
            disp(k1 - k3 + kp + q) *
            fermidist(mu - disp(k1 - k3 + kp + q), beta) - (
                u *
                disp(k3 + q) *
                disp(k1 - k3 + kp + q) *
                fermidist(mu - disp(k1 - k3 + kp + q), beta)
            )
        )
    ) *
    (
        (w + disp(k3) - disp(k3 + q)) *
        (mu + vp + w + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q)) *
        (
            2mu + v + vp + w - disp(k2) + disp(-1k + k1 + k2) - disp(k3) -
            disp(k1 - k3 + kp + q)
        ) *
        (2mu + 2vp + w - disp(-1k1 + k3 + kp) - disp(k1 - k3 + kp + q)) *
        (-1mu - vp + disp(k3) - disp(k1 + q) + disp(k1 - k3 + kp + q)) *
        (v - vp - disp(k2) + disp(k3) - disp(k + k1 - k2 + q) + disp(k1 - k3 + kp + q))
    ) ^ -1
end
