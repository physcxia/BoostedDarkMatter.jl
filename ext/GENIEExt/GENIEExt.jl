module GENIEExt

using BoostedDarkMatter
import BoostedDarkMatter as bd
import HEPUnits as units

export XSecGENIE
export dp

include(ENV["DARKPROPJL"])
using .DarkProp; const dp = DarkProp


Base.@kwdef mutable struct XSecGENIE{T <: Number, U <: Number} <: XSec
    mchi::T = 0.1 * units.GeV
    sigma0::U = 1e-32 * units.cm2
    process::String = "DMCEL"
end

function bd.dxsecdT4(xsec::XSecGENIE, T4, T1, m1, m2, target; kwargs...)
    Z = particle_Z(target)
    A = particle_A(target)
    if m2 == dmmass(xsec)
        mchi = m2
        mi = m1
        Tchi = T4
        Ti = T1
        Q2 = 2mchi * Tchi
        γ = (Ti + mi) / mi
        Echi_bar = γ * mchi
        sigma = 2mchi  # Jacobian
        # println(Q2, " ", Echi_bar)
        if xsec.process == "DMEL"
            sigma *= dp.dXSec_dQ2_DMEL(Q2, Echi_bar, Z, A)
        elseif xsec.process == "DMCEL"
            mi > 1 * units.GeV || return zero(sigma)  # exclude free nucleon
            sigma *= dp.dXSec_dQ2_DMCEL(Q2, Echi_bar, Z, A)
        elseif xsec.process == "DMRES"
            sigma *= dp.dXSec_dQ2_DMRES(Q2, Echi_bar, Z, A)
        elseif xsec.process == "DMDIS"
            sigma *= dp.dXSec_dQ2_DMDIS(Q2, Echi_bar, Z, A)
        elseif lowercase(xsec.process) == "all"
            sigma *= (
                dp.dXSec_dQ2_DMEL(Q2, Echi_bar, Z, A)
                + dp.dXSec_dQ2_DMCEL(Q2, Echi_bar, Z, A)
                + dp.dXSec_dQ2_DMRES(Q2, Echi_bar, Z, A)
                + dp.dXSec_dQ2_DMDIS(Q2, Echi_bar, Z, A)
            )
        else
            error("Unknown xsec.process: $(xsec.process)")
        end
        return sigma
    elseif m1 == dmmass(xsec)
        Echi = T1 + m1
        Tchi = Echi - T4 # final state DM kinetic energy
        if xsec.process == "DMEL"
            return dp.dXSec_dTchi_DMEL(Tchi, Echi, Z, A)
        elseif xsec.process == "DMCEL"
            m2 > 1 * units.GeV || return zero(sigma)
            return dp.dXSec_dTchi_DMCEL(Tchi, Echi, Z, A)
        elseif xsec.process == "DMRES"
            return dp.dXSec_dTchi_DMRES(Tchi, Echi, Z, A)
        elseif xsec.process == "DMDIS"
            return dp.dXSec_dTchi_DMDIS(Tchi, Echi, Z, A)
        elseif lowercase(xsec.process) == "all"
            return (
                dp.dXSec_dQ2_DMEL(Tchi, Echi, Z, A)
                + dp.dXSec_dQ2_DMCEL(Tchi, Echi, Z, A)
                + dp.dXSec_dQ2_DMRES(Tchi, Echi, Z, A)
                + dp.dXSec_dQ2_DMDIS(Tchi, Echi, Z, A)
            )
        else
            error("Unknown xsec.process: $(xsec.process)")
        end
    else
        error("DM mass mismatch")
    end
end

function bd.Kinematics.T1_min(xsec::XSecGENIE, T4, m1, m2)
    xsec.process != "DMCEL" || return T1_min(T4, m1, m2)
    mpi = 140*units.MeV
    mN = 0.938  # TODO use constants
    mN2 = mN^2
    Wmin = mN + mpi
    W2 = Wmin^2
    # DM scatter off nucleus
    if dmmass(xsec) == m1
        xsec.process != "DMEL" || return T1_min(T4, m1, mN)
        ν = T4
        mchi = m1
        ν > (W2 - mN2) / (2*(mN - mchi)) || return Inf
        Emin = 1/2 * (ν + sqrt(((mN + ν)^2 - W2) * (4*mchi^2 + mN2 + 2 * mN * ν - W2)
                               / (mN^2+2*mN*ν-W2)))
        return Emin - mchi
    # nucleus scatter off DM
    elseif dmmass(xsec) == m2
        xsec.process != "DMEL" || return m1 / mN * T1_min(T4, mN, m2)
        Q2 = 2m2 * T4
        mchi = m2
        # nucleus rest frame
        ν = (Q2 + W2 - mN2) / (2 * mN)
        Echimin = 1/2 * (ν + sqrt((ν^2 + Q2) * (1 + 4 * mchi^2 / Q2)))
        #  Echimin = ((Q2 * (Q2 + W2 - mN2)
        #              + sqrt(Q2*(Q2+4*mchi^2)*((Q2+W2)^2 + 2*(Q2-W2)*mN2+mN2^2)))
        #             / (4*Q2*mN))
        # boost to DM rest frame
        γ = Echimin / mchi
        E1min = γ * m1  # note here m1 is not mN
        return E1min - m1
    else
        error("DM mass incorrect")
    end
end

function bd.T4_max(xsec::XSecGENIE, T1, m1, m2)
    xsec.process != "DMCEL" || return T4_max(T1, m1, m2)
    mN = 0.9396  # TODO use constants
    Wcut = 1.7units.GeV
    Q2cut = 50*units.GeV^2
    # DM scatters off nucleus
    if dmmass(xsec) == m1
        xsec.process != "DMEL" || return T4_max(T1, m1, mN)
        # A and Z are not so important
        A = div(m2, 0.938units.GeV)
        Z = div(A, 2)
        Echi = T1 + m1
        if xsec.process == "DMRES"
            Echipmin, _ = dp.limits_Echip_DMRES(Echi, Z, A, Wcut, Q2cut)
        elseif xsec.process == "DMDIS"
            Echipmin, _ = dp.limits_Echip_DMDIS(Echi, Z, A)
        end
        return Echi - Echipmin
    # nucleus scatters off DM
    elseif dmmass(xsec) == m2
        xsec.process != "DMEL" || return T4_max(T1, mN, m2)
        mchi = m2
        # nucleus rest frame
        Echi = (T1 + m1) / m1 * mchi
        _, Q2max = dp.limits_Q2(Echi, mN, mchi)
        return Q2max / (2mchi)
    else
        error("DM mass incorrect")
    end
end

end
