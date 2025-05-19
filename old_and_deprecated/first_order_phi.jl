function phi1_ph(u, mu, beta, w, q, k1)
    dim = length(k1)
    cnorm = 1 / (2pi)^dim
    prefactor = cnorm * u^2
    matsum = matsubara_sum_fermi(beta, disp(k1) - mu, disp(k1 + q) - mu - w)
    return prefactor * matsum
end

function phi1_s(u, mu, beta, w, q, k1)
    dim = length(k1)
    cnorm = 1 / (2pi)^dim
    prefactor = 2 * cnorm * u^2
    matsum = matsubara_sum_fermi(beta, disp(k1) - mu, -disp(-k1 + q) + mu + w)
    return prefactor * matsum
end

function phi1_t(u, mu, beta, w, q, k1)
    return complex(0.0)
end
