""" Mass from A. """
nucleus_mass(A) = A != 1 ? ATOMIC_MASS * A : PROTON_MASS

""" Lorentz gamma ``1/\\sqrt{1-v^2}`` with speed of light c = 1. """
lorentz_gamma(v) = 1 / sqrt((1 + v) * (1 - v))

""" velocity to kinetic energy """
v_to_Ek(v, m) = (lorentz_gamma(v) - 1) * m

""" kinetic energy to velocity """
Ek_to_v(Ek, m) = sqrt(Ek * (2m + Ek)/((m + Ek)^2))


""" convert flux to f(v) """
function flux_to_vdist(Ek, flux, m)
    v = Ek_to_v(Ek, m)
    return flux * m^2 * lorentz_gamma(v)^3
end

"""Convert 1D velocity distribution f(v) to flux. """
vdist_to_flux(v, fv, m, rho0) = rho0 * fv / (m^2 * lorentz_gamma(v)^3)


"""Convert 3D velocity distribution f(v, theta, phi) to flux. """
vdist_to_flux_3d(v, fv, m, rho0) = vdist_to_flux(v, fv, m, rho0) * v^2


# ==================
# Rigidity to energy
# ==================

R_to_Et(R, Z, A) = sqrt((Z*R)^2 + nucleus_mass(A)^2)
R_to_Ek(R, Z, A) = sqrtm1((Z*R)^2, nucleus_mass(A))
R_to_Ekn(R, Z, A) = R_to_Ek(R, Z, A) / A


dEk_dR(R, Z, A) = Z^2 * R / R_to_Et(R, Z, A)
dEt_dR(R, Z, A) = dEk_dR(R, Z, A)
dEkn_dR(R, Z, A) = dEk_dR(R, Z, A) / A

flux_R_to_Ek(R, flux, Z, A) = flux / dEk_dR(R, Z, A)
flux_R_to_Et(R, flux, Z, A) = flux / dEt_dR(R, Z, A)
flux_R_to_Ekn(R, flux, Z, A) = flux / dEkn_dR(R, Z, A)


spectrum_R_to_Ek(R, flux, Z, A) = (R_to_Ek(R, Z, A), flux_R_to_Ek(R, flux, Z, A))
spectrum_R_to_Ekn(R, flux, Z, A) = (R_to_Ekn(R, Z, A), flux_R_to_Ekn(R, flux, Z, A))
spectrum_R_to_Et(R, flux, Z, A) = (R_to_Et(R, Z, A), flux_R_to_Et(R, flux, Z, A))


# ==================
# energy to rigidity
# ==================

Ek_to_R(Ek, Z, A) = sqrt(Ek*(Ek + 2*nucleus_mass(A))) / Z
Ekn_to_R(Ekn, Z, A) = Ek_to_R(A * Ekn, Z, A)
Et_to_R(Et, Z, A) = sqrt(Et^2 - nucleus_mass(A)^2) / Z


dR_dEk(Ek, Z, A) = (Ek + nucleus_mass(A)) / (Z * Z * Ek_to_R(Ek, Z, A))
dR_dEkn(Ekn, Z, A) = A * dR_dEk(A * Ekn, Z, A)
dR_dEt(Et, Z, A) = Et / (Z^2 * Et_to_R(Et, Z, A))

flux_Ekn_to_R(Ekn, flux, Z, A) = flux / dR_dEkn(Ekn, Z, A)
flux_Ek_to_R(Ek, flux, Z, A) = flux / dR_dEk(Ek, Z, A)
flux_Et_to_R(Et, flux, Z, A) = flux / dR_dEt(Et, Z, A)

spectrum_Ekn_to_R(Ekn, flux, Z, A) = (Ekn_to_R(Ekn, Z, A), flux_Ekn_to_R(Ekn, flux, Z, A))
