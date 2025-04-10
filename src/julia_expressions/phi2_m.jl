function phi2_m(
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
    -(
        (
            u ^ 2 *
            bosedist(-2mu + disp(k2) + disp(k + k1 - k2 + q), beta) *
            (
                u * fermidist(-1mu + disp(k2), beta) -
                (u * fermidist(mu - disp(k + k1 - k2 + q), beta))
            ) *
            (
                -(v * w) + vp * w + w * disp(k2) - (v * disp(k3)) + vp * disp(k3) -
                (w * disp(k3)) + disp(k2) * disp(k3) - disp(k3) ^ 2 +
                w * disp(k + k1 - k2 + q) +
                disp(k3) * disp(k + k1 - k2 + q) +
                v * disp(k3 + q) - (vp * disp(k3 + q)) - (disp(k2) * disp(k3 + q)) +
                disp(k3) * disp(k3 + q) - (disp(k + k1 - k2 + q) * disp(k3 + q)) -
                (w * disp(k1 - k3 + kp + q)) - (disp(k3) * disp(k1 - k3 + kp + q)) +
                disp(k3 + q) * disp(k1 - k3 + kp + q) +
                u * v * fermidist(-1mu + disp(k3), beta) -
                (u * vp * fermidist(-1mu + disp(k3), beta)) -
                (u * w * fermidist(-1mu + disp(k3), beta)) -
                (u * disp(k2) * fermidist(-1mu + disp(k3), beta)) -
                (u * disp(k + k1 - k2 + q) * fermidist(-1mu + disp(k3), beta)) +
                u * disp(k3 + q) * fermidist(-1mu + disp(k3), beta) +
                u * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta) -
                (u * v * fermidist(-1mu + disp(k3 + q), beta)) +
                u * vp * fermidist(-1mu + disp(k3 + q), beta) +
                u * disp(k2) * fermidist(-1mu + disp(k3 + q), beta) -
                (u * disp(k3) * fermidist(-1mu + disp(k3 + q), beta)) +
                u * disp(k + k1 - k2 + q) * fermidist(-1mu + disp(k3 + q), beta) -
                (u * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3 + q), beta)) +
                u * w * fermidist(mu - disp(k1 - k3 + kp + q), beta) +
                u * disp(k3) * fermidist(mu - disp(k1 - k3 + kp + q), beta) -
                (u * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta))
            )
        ) *
        (
            (mu + v + w + disp(k1) - disp(k2) - disp(k + k1 - k2 + q)) *
            (mu + v - disp(k2) + disp(k1 + q) - disp(k + k1 - k2 + q)) *
            (w + disp(k3) - disp(k3 + q)) *
            (v - vp - disp(k2) + disp(k3) - disp(k + k1 - k2 + q) + disp(k1 - k3 + kp + q))
        ) ^ -1
    ) +
    (
        u ^ 2 *
        bosedist(-2mu + disp(k3) + disp(k1 - k3 + kp + q), beta) *
        (
            -1v + vp + disp(k2) - disp(k3) + disp(k + k1 - k2 + q) -
            disp(k1 - k3 + kp + q) + u * fermidist(-1mu + disp(k2), beta) -
            (u * fermidist(mu - disp(k + k1 - k2 + q), beta))
        ) *
        (
            u * w * fermidist(-1mu + disp(k3), beta) +
            u * disp(k3) * fermidist(-1mu + disp(k3), beta) -
            (u * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) -
            (u * w * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
            (u * disp(k3) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) +
            u * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta)
        )
    ) *
    (
        (w + disp(k3) - disp(k3 + q)) *
        (mu + vp + w + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q)) *
        (-1mu - vp + disp(k3) - disp(k1 + q) + disp(k1 - k3 + kp + q)) *
        (v - vp - disp(k2) + disp(k3) - disp(k + k1 - k2 + q) + disp(k1 - k3 + kp + q))
    ) ^ -1 +
    (
        u ^ 2 *
        fermidist(-1mu + disp(k1), beta) *
        (
            -1mu - v - w - disp(k1) +
            disp(k2) +
            disp(k + k1 - k2 + q) +
            u * fermidist(-1mu + disp(k2), beta) -
            (u * fermidist(mu - disp(k + k1 - k2 + q), beta))
        ) *
        (
            -(mu * w) - (vp * w) - w ^ 2 - (w * disp(k1)) - (mu * disp(k3)) -
            (vp * disp(k3)) - (disp(k1) * disp(k3)) +
            disp(k3) ^ 2 +
            mu * disp(k3 + q) +
            vp * disp(k3 + q) +
            w * disp(k3 + q) +
            disp(k1) * disp(k3 + q) - (disp(k3) * disp(k3 + q)) +
            w * disp(k1 - k3 + kp + q) +
            disp(k3) * disp(k1 - k3 + kp + q) - (disp(k3 + q) * disp(k1 - k3 + kp + q)) +
            mu * u * fermidist(-1mu + disp(k3), beta) +
            u * vp * fermidist(-1mu + disp(k3), beta) +
            2 * u * w * fermidist(-1mu + disp(k3), beta) +
            u * disp(k1) * fermidist(-1mu + disp(k3), beta) -
            (u * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) -
            (u * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta)) -
            (mu * u * fermidist(-1mu + disp(k3 + q), beta)) -
            (u * vp * fermidist(-1mu + disp(k3 + q), beta)) -
            (u * w * fermidist(-1mu + disp(k3 + q), beta)) -
            (u * disp(k1) * fermidist(-1mu + disp(k3 + q), beta)) +
            u * disp(k3) * fermidist(-1mu + disp(k3 + q), beta) +
            u * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3 + q), beta) -
            (u * w * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
            (u * disp(k3) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) +
            u * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta)
        )
    ) *
    (
        (w + disp(k1) - disp(k1 + q)) *
        (mu + v + w + disp(k1) - disp(k2) - disp(k + k1 - k2 + q)) *
        (w + disp(k3) - disp(k3 + q)) *
        (mu + vp + w + disp(k1) - disp(k3) - disp(k1 - k3 + kp + q))
    ) ^ -1 - (
        (
            u ^ 2 *
            fermidist(-1mu + disp(k1 + q), beta) *
            (
                -1mu - v + disp(k2) - disp(k1 + q) +
                disp(k + k1 - k2 + q) +
                u * fermidist(-1mu + disp(k2), beta) -
                (u * fermidist(mu - disp(k + k1 - k2 + q), beta))
            ) *
            (
                -(mu * w) - (vp * w) - (mu * disp(k3)) - (vp * disp(k3)) +
                w * disp(k3) +
                disp(k3) ^ 2 - (w * disp(k1 + q)) - (disp(k3) * disp(k1 + q)) +
                mu * disp(k3 + q) +
                vp * disp(k3 + q) - (disp(k3) * disp(k3 + q)) +
                disp(k1 + q) * disp(k3 + q) +
                w * disp(k1 - k3 + kp + q) +
                disp(k3) * disp(k1 - k3 + kp + q) -
                (disp(k3 + q) * disp(k1 - k3 + kp + q)) +
                mu * u * fermidist(-1mu + disp(k3), beta) +
                u * vp * fermidist(-1mu + disp(k3), beta) +
                u * w * fermidist(-1mu + disp(k3), beta) +
                u * disp(k1 + q) * fermidist(-1mu + disp(k3), beta) -
                (u * disp(k3 + q) * fermidist(-1mu + disp(k3), beta)) -
                (u * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3), beta)) -
                (mu * u * fermidist(-1mu + disp(k3 + q), beta)) -
                (u * vp * fermidist(-1mu + disp(k3 + q), beta)) +
                u * disp(k3) * fermidist(-1mu + disp(k3 + q), beta) -
                (u * disp(k1 + q) * fermidist(-1mu + disp(k3 + q), beta)) +
                u * disp(k1 - k3 + kp + q) * fermidist(-1mu + disp(k3 + q), beta) -
                (u * w * fermidist(mu - disp(k1 - k3 + kp + q), beta)) -
                (u * disp(k3) * fermidist(mu - disp(k1 - k3 + kp + q), beta)) +
                u * disp(k3 + q) * fermidist(mu - disp(k1 - k3 + kp + q), beta)
            )
        ) *
        (
            (w + disp(k1) - disp(k1 + q)) *
            (mu + v - disp(k2) + disp(k1 + q) - disp(k + k1 - k2 + q)) *
            (w + disp(k3) - disp(k3 + q)) *
            (mu + vp - disp(k3) + disp(k1 + q) - disp(k1 - k3 + kp + q))
        ) ^ -1
    )
end
