using MultiQuad: dblquad, quad
using NumericalTools: sqrtm1

abstract type XSec end
abstract type XSecElastic <: XSec end
abstract type XSecDMElectronElastic <: XSecElastic end

Broadcast.broadcastable(xsec::XSec) = Ref(xsec)

dmmass(xsec::XSec) = xsec.mchi
xsec0(xsec::XSec) = xsec.sigma0
mediatormass(xsec::XSec) = xsec.mmed

dmmass!(xsec::XSec, mchi) = xsec.mchi = mchi
xsec0!(xsec::XSec, sigma0) = xsec.sigma0 = sigma0
mediatormass!(xsec::XSec, mmed) = xsec.mmed = mmed
function set_parameters!(xsec::XSec; params...)
    for (key, value) in params
        hasproperty(xsec, key) && setproperty!(xsec, key, value)
    end
end

Kinematics.T1_min(::XSec, T4, m1, m2) = error("unimplemented")
Kinematics.T1_min(::XSecElastic, T4, m1, m2) = T1_min(T4, m1, m2)

dm_formfactor(::XSec, Q2) = error("unsupproted form factor")
dm_formfactor(::XSec, Q2, T1) = error("unsupproted form factor")


"""
    dxsec(xsec::XSec, T4, T1, m1, m2)

Differential cross section ``\\frac{dσ}{dT_4}(T_4, T_1, m_1, m_2)``.

# Arguments

- `T4`: Recoil energy of the target particle 2.
- `T1`: Kinetic energy of the incident particle 1.
- `m1`: Mass of the incident particle 1.
- `m2`: Mass of the target particle 2.

"""
dxsecdT4(::XSec, T4, T1, m1, m2) = error("unimplemented")

Base.@kwdef mutable struct
XSecDMElectronVectorMediator{T <: Number, U <: Number, V <: Number} <: XSecDMElectronElastic
    sigma0::T = 1e-30
    mchi::U = 1e-6
    mmed::U = 1e-6
    me::U = ELECTRON_MASS
    alpha::V = 1/137
end
function XSecDMElectronVectorMediator(sigma0, mchi, mmed, me, alpha)
    return XSecDMElectronVectorMediator(sigma0, promote(mchi, mmed, me)..., alpha)
end
function dxsecdT4(xsec::XSecDMElectronVectorMediator, T4, T1, m1, m2)
    μχe = reduce_m(xsec.mchi, xsec.me)
    return (xsec.sigma0 * ((xsec.alpha*xsec.me)^2 + xsec.mmed^2)^2 / μχe^2
            * (2*m2*(m1 + T1)^2 - T4*((m1 + m2)^2 + 2*m2*T1) + m2*T4*T4)
            / (4*T1*(2*m1 + T1)*(2*m2*T4+xsec.mmed*xsec.mmed)^2))
end
function dm_formfactor(xsec::XSecDMElectronVectorMediator, Q2, Tchi)
    # return ((xsec.alpha^2 * me^2 + xsec.mmed^2)^2 / (Q2 + xsec.mmed^2)^2
    #     * (8me^2 * (Tchi + xsec.mchi)^2 - 2me*Q2*((xsec.mchi + me)^2 + 2me * Tchi) + Q2^2)
    #     / (8me^2 * xsec.mchi^2))
    # non-relativistic approximation
    return ((xsec.alpha^2 * xsec.me^2 + xsec.mmed^2)^2 / (Q2 + xsec.mmed^2)^2
        * (Tchi + xsec.mchi)^2 / xsec.mchi^2)
end


"""
    dxsecdT4(xsec::XSecDMElectronScalarMediator, T4, T1, m1, m2)

Scalar mediated DM-electron scattering cross section.

"""
Base.@kwdef mutable struct
XSecDMElectronScalarMediator{T <: Number, U <: Number, V <: Number} <: XSecDMElectronElastic
    sigma0::T = 1e-30
    mchi::U = 0.1
    mmed::U = 1e-6
    alpha::V = 1/137
end
function dxsecdT4(xsec::XSecDMElectronScalarMediator, T4, T1, m1, m2)
    me = ELECTRON_MASS
    μχe = reduce_m(xsec.mchi, me)
    Q2 = 2*m2*T4
    return (xsec.sigma0 * (xsec.mmed^2+(xsec.alpha * me)^2)^2 / (32*m2*μχe^2*T1*(T1 + 2*m1))
            * (4*m1^2 + Q2)*(4*m2^2 + Q2)/(xsec.mmed^2 + Q2)^2)
end


struct IronizationFormFactor
    ks::AbstractVector
    qs::AbstractVector
    formfactor::Function
end


Base.@kwdef mutable struct XSecDMElectronBound{T <: XSec} <: XSec
    freexsec::T = XSecDMElectronVectorMediator()
    ironff = (_, _) -> 1
end
# dmmass!(xsec::XSecDMElectronBound, mchi) = dmmass!(xsec.freexsec, mchi)
# xsec0!(xsec::XSecDMElectronBound, sigma0) = xsec0!(xsec.freexsec, sigma0)
function set_parameters!(xsec::XSecDMElectronBound; params...)
    set_parameters!(xsec, params...)
    set_parameters!(xsec.freexsec, params...)
end
function Base.getproperty(xsec::XSecDMElectronBound, name::Symbol)
    hasproperty(xsec, name) && return getfield(xsec, name)
    return getproperty(xsec.freexsec, name)
end
function Base.setproperty!(xsec::XSecDMElectronBound, name::Symbol, value)
    hasproperty(xsec, name) && return setfield!(xsec, name, value)
    return setproperty!(xsec.freexsec, name, value)
end


"""
    recoil_spectrum(xsec::XSec, T4, flux1, T1max, m1, m2; attenuation=nothing)

Recoil spectrum of the target at rest in the scattering.

# Arguments

- `xsec::XSec`: Differential cross section.
- `T4`: Recoil energy of the target particle 2.
- `flux1`: Flux of the incident particle 1.
- `T1max`: Maximal energy of flux1.
- `m1`: Mass of the incident particle 1.
- `m2`: Mass of the target particle 2.

"""
function recoil_spectrum(xsec::XSec, T4, flux1, T1max, m1, m2; attenuation=nothing)
    T1min = T1_min(xsec, T4, m1, m2)
    T1min < T1max || return 0

    if !isnothing(attenuation)
        # ekin_in_min = attenuation.T0Tz(ekin_in_min)
        ekin_in_min0 = attenuation.T0Tz(ekin_in_min)
        ekin_in_min = ekin_in_min0
    end

    Eunit = oneunit(T1min)
    if isnothing(attenuation)
        rate_res = quad(log(T1min / Eunit), log(T1max / Eunit)) do logT1
            T1 = exp(logT1) * Eunit
            dxsecdT4(xsec, T4, T1, m1, m2) * flux1(T1) * T1
        end
    else
        rate_res = quad(log(T1min / Eunit), log(T1max / Eunit)) do logT1
            T1 = exp(logT1) * Eunit
            T1z = attenuation(T1)
            dxsecdT4(xsec, T4, T1z, m1, m2) * flux1(T1) * T1
        end
    end

    return rate_res[1] / m2
end
function recoil_spectrum(
    xsec::XSec,
    T4::AbstractVector,
    T1::AbstractVector,
    flux1::AbstractVector,
    m1,
    m2;
    attenuation=nothing,
)
    flux_in = loginterpolator(T1, flux1)
    return recoil_spectrum.(xsec::XSec, T4, Ref(flux_in), maximum(T1), m1, m2;
                            attenuation=attenuation)
end
function recoil_spectrum(
    xsec::XSecDMElectronBound, Tr, Eb, fluxchi, mchi, me, qmax, Tchimax; unit=1, kwargs...
)
    Eunit = oneunit(Tr)
    ΔE = Tr + Eb
    logqmax = log(qmax / Eunit)
    logqmin = log(ΔE / Eunit) + sqrt(eps(ΔE / Eunit))
    logqmin < logqmax || throw(DomainError("logqmin = $logqmin < logqmax = $logqmax"))
    logTchimin(logq) = log(T1_min_q(exp(logq) * Eunit, ΔE, mchi) / Eunit)
    logtchimax = log(Tchimax / Eunit)
    logTchimax(_) = logtchimax

    ke = sqrt(Tr * (Tr + 2me))
    μχe = reduce_m(mchi, me)
    factor = xsec.sigma0 / (8 * Tr * μχe^2 * me)
    unitzero = zero(unit / oneunit(factor))

    rate_res = dblquad(
        logqmin, logqmax, logTchimin, logTchimax; kwargs...
    ) do logTchi, logq
        logTchi < logtchimax || return unitzero
        Tchi = exp(logTchi) * Eunit
        q = exp(logq) * Eunit
        Q2 = q^2 - ΔE^2
        pchi = sqrt(Tchi * (Tchi + 2mchi))
        return (q * dm_formfactor(xsec.freexsec, Q2, Tchi)
            * xsec.ironff(ke, q)
            * mchi^2 / (pchi*(Tchi+mchi))
            * fluxchi(Tchi)
            * Tchi * q)
    end

    return rate_res[1] * factor
end
