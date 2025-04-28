function gamma2_m(u, mu, beta, v, vp, w, k::SVector, kp::SVector, q::SVector, k1::SVector)
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
    ) * (2mu + v + vp + w - disp(k1) - disp(k - k1 + kp + q)) ^ -1
end
