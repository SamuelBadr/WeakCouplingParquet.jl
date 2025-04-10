function phi2_t(
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
        u ^ 4 *
        bosedist(disp(k3) - disp(-1k1 + k3 + kp), beta) *
        (
            2 * vp * fermidist(-1mu + disp(k3), beta) -
            (w * fermidist(-1mu + disp(k3), beta)) -
            (disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) +
            disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta) -
            (2 * vp * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
            w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
            (disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta))
        ) *
        (
            2 * vp * fermidist(-1mu + disp(k2), beta) -
            (w * fermidist(-1mu + disp(k2), beta)) -
            (disp(-1k + k1 + k2) * fermidist(-1mu + disp(k2), beta)) +
            2 * disp(k3) * fermidist(-1mu + disp(k2), beta) -
            (2 * disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k2), beta)) +
            disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(k2), beta) -
            (v * fermidist(-1mu + disp(-1k + k1 + k2), beta)) -
            (vp * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
            w * fermidist(-1mu + disp(-1k + k1 + k2), beta) +
            disp(k2) * fermidist(-1mu + disp(-1k + k1 + k2), beta) -
            (disp(k3) * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
            disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(-1k + k1 + k2), beta) -
            (disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
            v * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta) -
            (vp * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) -
            (disp(k2) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) +
            disp(-1k + k1 + k2) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta) -
            (disp(k3) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) +
            disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
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
        (2vp - w - disp(-1k1 + k3 + kp) + disp(-1k1 + k3 - kp + q))
    ) ^ -1 - (
        (
            u ^ 4 *
            bosedist(-disp(k2) + disp(-1k + k1 + k2), beta) *
            (
                2 * v * fermidist(-1mu + disp(k2), beta) -
                (w * fermidist(-1mu + disp(k2), beta)) -
                (2 * disp(k2) * fermidist(-1mu + disp(k2), beta)) +
                disp(-1k + k1 + k2) * fermidist(-1mu + disp(k2), beta) +
                disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(k2), beta) -
                (2 * v * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
                w * fermidist(-1mu + disp(-1k + k1 + k2), beta) +
                2 * disp(k2) * fermidist(-1mu + disp(-1k + k1 + k2), beta) -
                (disp(-1k + k1 + k2) * fermidist(-1mu + disp(-1k + k1 + k2), beta)) -
                (disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(-1k + k1 + k2), beta))
            ) *
            (
                -2 * vp * fermidist(-1mu + disp(k3), beta) +
                w * fermidist(-1mu + disp(k3), beta) +
                disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta) -
                (disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta)) +
                v * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
                vp * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
                (w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
                (disp(k2) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
                disp(-1k + k1 + k2) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
                (disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
                disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
                (v * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
                vp * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) +
                disp(k2) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
                (disp(-1k + k1 + k2) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
                disp(k3) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
                (disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta))
            )
        ) *
        (
            2 *
            (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) *
            (v - vp - disp(k2) + disp(-1k + k1 + k2) - disp(k3) + disp(-1k1 + k3 + kp)) *
            (mu - v + w + disp(k2) - disp(-1k + k1 + k2) - disp(-1k1 + q)) *
            (2v - w - (2 * disp(k2)) + disp(-1k + k1 + k2) + disp(-1k - k1 + k2 + q)) *
            (
                v + vp - w - disp(k2) + disp(-1k + k1 + k2) - disp(k3) +
                disp(-1k1 + k3 - kp + q)
            )
        ) ^ -1
    ) +
    (
        u ^ 4 *
        fermidist(-1mu + disp(k1), beta) *
        (
            2 * mu * fermidist(-1mu + disp(k2), beta) +
            w * fermidist(-1mu + disp(k2), beta) -
            (2 * disp(k1) * fermidist(-1mu + disp(k2), beta)) +
            disp(-1k + k1 + k2) * fermidist(-1mu + disp(k2), beta) -
            (disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(k2), beta)) -
            (mu * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
            v * fermidist(-1mu + disp(-1k + k1 + k2), beta) -
            (w * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
            disp(k1) * fermidist(-1mu + disp(-1k + k1 + k2), beta) -
            (disp(k2) * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
            disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(-1k + k1 + k2), beta) -
            (mu * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) -
            (v * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) +
            disp(k1) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta) +
            disp(k2) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta) -
            (disp(-1k + k1 + k2) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta))
        ) *
        (
            2 * vp * fermidist(-1mu + disp(k3), beta) -
            (w * fermidist(-1mu + disp(k3), beta)) -
            (disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) +
            disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta) +
            mu * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
            (vp * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
            w * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
            (disp(k1) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
            disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
            (disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
            (mu * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) -
            (vp * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
            disp(k1) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
            (disp(k3) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
            disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
        )
    ) *
    (
        2 *
        (mu + v - disp(k1) - disp(k2) + disp(-1k + k1 + k2)) *
        (mu + vp - disp(k1) + disp(k3) - disp(-1k1 + k3 + kp)) *
        (2mu + w - disp(k1) - disp(-1k1 + q)) *
        (mu - v + w - disp(k1) + disp(k2) - disp(-1k - k1 + k2 + q)) *
        (mu - vp + w - disp(k1) + disp(k3) - disp(-1k1 + k3 - kp + q))
    ) ^ -1 - (
        (
            u ^ 4 *
            fermidist(mu - disp(-1k1 + q), beta) *
            (
                2 * mu * fermidist(-1mu + disp(k2), beta) +
                w * fermidist(-1mu + disp(k2), beta) -
                (disp(-1k + k1 + k2) * fermidist(-1mu + disp(k2), beta)) -
                (2 * disp(-1k1 + q) * fermidist(-1mu + disp(k2), beta)) +
                disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(k2), beta) -
                (mu * fermidist(-1mu + disp(-1k + k1 + k2), beta)) -
                (v * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
                disp(k2) * fermidist(-1mu + disp(-1k + k1 + k2), beta) +
                disp(-1k1 + q) * fermidist(-1mu + disp(-1k + k1 + k2), beta) -
                (disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(-1k + k1 + k2), beta)) -
                (mu * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) +
                v * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta) -
                (w * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) -
                (disp(k2) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) +
                disp(-1k + k1 + k2) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta) +
                disp(-1k1 + q) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)
            ) *
            (
                -2 * vp * fermidist(-1mu + disp(k3), beta) +
                w * fermidist(-1mu + disp(k3), beta) +
                disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta) -
                (disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta)) +
                mu * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
                vp * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
                (disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
                (disp(-1k1 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
                disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
                (mu * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
                vp * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
                (w * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
                disp(k3) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
                (disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
                disp(-1k1 + q) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)
            )
        ) *
        (
            2 *
            (2mu + w - disp(k1) - disp(-1k1 + q)) *
            (mu - v + w + disp(k2) - disp(-1k + k1 + k2) - disp(-1k1 + q)) *
            (mu - vp + w - disp(k3) + disp(-1k1 + k3 + kp) - disp(-1k1 + q)) *
            (mu + v - disp(k2) - disp(-1k1 + q) + disp(-1k - k1 + k2 + q)) *
            (mu + vp - disp(k3) - disp(-1k1 + q) + disp(-1k1 + k3 - kp + q))
        ) ^ -1
    ) +
    (
        u ^ 4 *
        bosedist(disp(k2) - disp(-1k - k1 + k2 + q), beta) *
        (
            2 * v * fermidist(-1mu + disp(k2), beta) -
            (w * fermidist(-1mu + disp(k2), beta)) -
            (2 * disp(k2) * fermidist(-1mu + disp(k2), beta)) +
            disp(-1k + k1 + k2) * fermidist(-1mu + disp(k2), beta) +
            disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(k2), beta) -
            (2 * v * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) +
            w * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta) +
            2 * disp(k2) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta) -
            (disp(-1k + k1 + k2) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) -
            (disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta))
        ) *
        (
            2 * vp * fermidist(-1mu + disp(k3), beta) -
            (w * fermidist(-1mu + disp(k3), beta)) -
            (disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) +
            disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta) +
            v * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
            (vp * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
            (disp(k2) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) +
            disp(k3) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) +
            disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta) -
            (disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(-1k1 + k3 + kp), beta)) -
            (v * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) -
            (vp * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
            w * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) +
            disp(k2) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
            (disp(k3) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
            disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
            (disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta))
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
        (v - vp - disp(k2) + disp(k3) + disp(-1k - k1 + k2 + q) - disp(-1k1 + k3 - kp + q))
    ) ^ -1 +
    (
        u ^ 4 *
        bosedist(disp(k3) - disp(-1k1 + k3 - kp + q), beta) *
        (
            2 * vp * fermidist(-1mu + disp(k2), beta) -
            (w * fermidist(-1mu + disp(k2), beta)) +
            disp(-1k + k1 + k2) * fermidist(-1mu + disp(k2), beta) -
            (2 * disp(k3) * fermidist(-1mu + disp(k2), beta)) -
            (disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(k2), beta)) +
            2 * disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k2), beta) +
            v * fermidist(-1mu + disp(-1k + k1 + k2), beta) -
            (vp * fermidist(-1mu + disp(-1k + k1 + k2), beta)) -
            (disp(k2) * fermidist(-1mu + disp(-1k + k1 + k2), beta)) +
            disp(k3) * fermidist(-1mu + disp(-1k + k1 + k2), beta) +
            disp(-1k - k1 + k2 + q) * fermidist(-1mu + disp(-1k + k1 + k2), beta) -
            (disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(-1k + k1 + k2), beta)) -
            (v * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) -
            (vp * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) +
            w * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta) +
            disp(k2) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta) -
            (disp(-1k + k1 + k2) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta)) +
            disp(k3) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta) -
            (disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(-1k - k1 + k2 + q), beta))
        ) *
        (
            2 * vp * fermidist(-1mu + disp(k3), beta) -
            (w * fermidist(-1mu + disp(k3), beta)) -
            (disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(k3), beta)) +
            disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(k3), beta) -
            (2 * vp * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta)) +
            w * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) +
            disp(-1k1 + k3 + kp) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta) -
            (disp(-1k1 + k3 - kp + q) * fermidist(-1mu + disp(-1k1 + k3 - kp + q), beta))
        )
    ) *
    (
        2 *
        (mu - vp + w - disp(k1) + disp(k3) - disp(-1k1 + k3 - kp + q)) *
        (
            v + vp - w - disp(k2) + disp(-1k + k1 + k2) - disp(k3) +
            disp(-1k1 + k3 - kp + q)
        ) *
        (2vp - w - disp(-1k1 + k3 + kp) + disp(-1k1 + k3 - kp + q)) *
        (mu + vp - disp(k3) - disp(-1k1 + q) + disp(-1k1 + k3 - kp + q)) *
        (
            -1v + vp + disp(k2) - disp(k3) - disp(-1k - k1 + k2 + q) +
            disp(-1k1 + k3 - kp + q)
        )
    ) ^ -1
end
