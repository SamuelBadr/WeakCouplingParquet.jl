function gamma2_d(u, mu, beta, v, vp, w, k, kp, q, k1)
    u - 2 * phi1_ph(u, mu, beta, vp - v, kp - k, k1) + 1 / 2 * phi1_s(u, mu, beta, w + v + vp, q + k + kp, k1)
end

function gamma2_m(u, mu, beta, v, vp, w, k, kp, q, k1)
    -u - 1 / 2 * phi1_s(u, mu, beta, w + v + vp, q + k + kp, k1)
end

function gamma2_s(u, mu, beta, v, vp, w, k, kp, q, k1)
    2u - phi1_ph(u, mu, beta, w - vp - v, q - kp - k, k1) - phi1_ph(u, mu, beta, vp - v, kp - k, k1)
end

function gamma2_t(u, mu, beta, v, vp, w, k, kp, q, k1)
    phi1_ph(u, mu, beta, w - vp - v, q - kp - k, k1) - phi1_ph(u, mu, beta, vp - v, kp - k, k1)
end


function full2_d(u, mu, beta, v, vp, w, k, kp, q, k1)
    gamma2_d(u, mu, beta, v, vp, w, k, kp, q, k1) + phi1_ph(u, mu, beta, w, q, k1)
end

function full2_m(u, mu, beta, v, vp, w, k, kp, q, k1)
    gamma2_m(u, mu, beta, v, vp, w, k, kp, q, k1) + phi1_ph(u, mu, beta, w, q, k1)
end

function full2_s(u, mu, beta, v, vp, w, k, kp, q, k1)
    gamma2_s(u, mu, beta, v, vp, w, k, kp, q, k1) + phi1_s(u, mu, beta, w, q, k1)
end

function full2_t(u, mu, beta, v, vp, w, k, kp, q, k1)
    gamma2_t(u, mu, beta, v, vp, w, k, kp, q, k1) + phi1_t(u, mu, beta, w, q, k1)
end


function phi2_d(u, mu, beta, v, vp, w, k, kp, q, k1, k2, k3)
    dim = length(k1)
    cnorm = 1 / (2pi)^dim
    chi0 = (disp(k1) - mu, disp(k1 + q) - mu - w)
    u^2 * matsubara_sum_fermi(beta, chi0...) +
    u^3 * cnorm * if w + disp(k3) - disp(k3 + q) == 0
        dnf(beta, disp(k3) - mu)
    else
        (nf(beta, disp(k3)) - nf(beta, disp(k3 + q) - mu)) / (w + disp(k3) - disp(k3 + q))
    end * matsubara_sum_fermi(beta, chi0...) -
    -2 * u^3 * cnorm * if -v1
end