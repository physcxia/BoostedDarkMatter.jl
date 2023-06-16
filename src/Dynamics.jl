using MultiQuad: dblquad, quad
using SpecialFunctions: sphericalbesselj

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
settarget!(::XSec, target) = nothing
gettarget(::XSec) = nothing
gettarget(::XSecDMElectronElastic) = "Electron"

Kinematics.T1_min(::XSec, T4, m1, m2) = error("unimplemented")
Kinematics.T1_min(::XSecElastic, T4, m1, m2) = T1_min(T4, m1, m2)

Kinematics.T4_max(::XSec, T1, m1, m2) = error("unimplemented")
Kinematics.T4_max(::XSecElastic, T1, m1, m2) = T4_max(T1, m1, m2)

dm_formfactor(::XSec, Q2) = error("unsupproted form factor")
dm_formfactor(::XSec, Q2, T1) = error("unsupproted form factor")


"""
    dxsecdT4(xsec::XSec, T4, T1, m1, m2)

Differential cross section ``\\frac{dσ}{dT_4}(T_4, T_1, m_1, m_2)``.

# Arguments

- `T4`: Recoil energy of the target particle 2.
- `T1`: Kinetic energy of the incident particle 1.
- `m1`: Mass of the incident particle 1.
- `m2`: Mass of the target particle 2.

"""
dxsecdT4(::XSec, T4, T1, m1, m2) = error("unimplemented")

function total_xsec(xsec::XSec, T1, m1, m2; Trcut=0*units.eV, options...)
    Tr_max = T4_max(xsec, T1, m1, m2)
    # res, err = quad(Trcut, Tr_max; options...) do T4
    #     return dxsecdT4(xsec, T4, T1, m1, m2)
    # end
    res, err = quad(log(Trcut), log(Tr_max); options...) do lnT4
        T4 = exp(lnT4)
        return dxsecdT4(xsec, T4, T1, m1, m2) * T4
    end
    # res_quad = quad(log(1e-300), log(Tr_max); options...) do lnT4
    #     T4 = exp(lnT4)
    #     return dxsecdT4(xsec, T4, T1, m1, m2) * T4
    # end
    return res, err
end


Base.@kwdef mutable struct
XSecDMElectronVectorMediator{T <: Number, U <: Number, V <: Number} <: XSecDMElectronElastic
    sigma0::T = 1e-30 * units.cm2
    mchi::U = 1 * units.keV
    mmed::U = 1 * units.keV
    me::U = ELECTRON_MASS
    alpha::V = 1/137
    limit_case::String = ""
end
function XSecDMElectronVectorMediator(sigma0, mchi, mmed, me, alpha)
    return XSecDMElectronVectorMediator(sigma0, promote(mchi, mmed, me)..., alpha)
end
function dxsecdT4(xsec::XSecDMElectronVectorMediator, T4, T1, m1, m2)
    μχe = reduce_m(xsec.mchi, xsec.me)
    if xsec.limit_case == "light"
        return (xsec.sigma0 * (xsec.alpha*xsec.me)^4 / μχe^2
                * (2*m2*(m1 + T1)^2 - T4*((m1 + m2)^2 + 2*m2*T1) + m2*T4*T4)
                / (4*T1*(2*m1 + T1)*(2*m2*T4)^2))
    elseif xsec.limit_case == "heavy"
        return (xsec.sigma0 / μχe^2
                * (2*m2*(m1 + T1)^2 - T4*((m1 + m2)^2 + 2*m2*T1) + m2*T4*T4)
                / (4*T1*(2*m1 + T1)))
    elseif isempty(xsec.limit_case)
        return (xsec.sigma0 * ((xsec.alpha*xsec.me)^2 + xsec.mmed^2)^2 / μχe^2
                * (2*m2*(m1 + T1)^2 - T4*((m1 + m2)^2 + 2*m2*T1) + m2*T4*T4)
                / (4*T1*(2*m1 + T1)*(2*m2*T4+xsec.mmed*xsec.mmed)^2))
    end
    throw(KeyError("Unknown limit_case: $(xsec.limit_case)"))
end
""" A Non-relativistic approximation for direct detection. """
function dm_formfactor(xsec::XSecDMElectronVectorMediator, Q2, Tchi)
    # return ((xsec.alpha^2 * me^2 + xsec.mmed^2)^2 / (Q2 + xsec.mmed^2)^2
    #     * (8me^2 * (Tchi + xsec.mchi)^2 - 2me*Q2*((xsec.mchi + me)^2 + 2me * Tchi) + Q2^2)
    #     / (8me^2 * xsec.mchi^2))
    # non-relativistic approximation
    if xsec.limit_case == "light"
        return (xsec.alpha^2 * xsec.me^2)^2 / Q2^2 * (Tchi + xsec.mchi)^2 / xsec.mchi^2
    elseif xsec.limit_case == "heavy"
        return (Tchi + xsec.mchi)^2 / xsec.mchi^2
    end
    return ((xsec.alpha^2 * xsec.me^2 + xsec.mmed^2)^2 / (Q2 + xsec.mmed^2)^2
        * (Tchi + xsec.mchi)^2 / xsec.mchi^2)
end


"""
    dxsecdT4(xsec::XSecDMElectronScalarMediator, T4, T1, m1, m2)

Scalar mediated DM-electron scattering cross section.

"""
Base.@kwdef mutable struct
XSecDMElectronScalarMediator{T <: Number, U <: Number, V <: Number} <: XSecDMElectronElastic
    sigma0::T = 1e-30 * units.cm2
    mchi::U = 0.1 * units.GeV
    mmed::U = 1e-6 * units.GeV
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


Base.@kwdef mutable struct XSecDMElectronBound{X <: XSec, F <: Function, T <: Number} <: XSec
    freexsec::X = XSecDMElectronVectorMediator()
    ionff::F = (_, _) -> 1
    Eb::T = 25.7units.eV
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


function T1_min(xsec::XSecDMElectronBound, T4, m1, m2)
    # FIXME
    # ΔE = T4 + xsec.Eb
    # q = sqrt(
    # p1min = p1_min_q(q, ΔE, m1)
    T1_min(T4, m1, m2)
end
""" using non-relativistic approximation """
function dxsecdT4(xsec::XSecDMElectronBound, T4, T1, m1, m2; kwargs...)
    ΔE = T4 + xsec.Eb
    T1 > ΔE || return zero(xsec0(xsec)/oneunit(T4))
    ke = sqrt(2m2 * T4)
    p1 = sqrt(T1 * (T1 + 2m1))
    T3 = T1 - ΔE
    p3 = sqrt(T3 * (T3 + 2m1))
    qmin = max(p1 - p3, ΔE)
    qmax = p1 + p3
    μχe = reduce_m(dmmass(xsec), xsec.me)
    q_int::typeof(T4), _ = quad(log(qmin), log(qmax); kwargs...) do logq
        q = exp(logq)
        Q2 = q^2 - ΔE^2
        return q^2 * dm_formfactor(xsec.freexsec, Q2, T1) * xsec.ionff(ke, q)
    end
    return xsec0(xsec) * m1^2 / (8 * μχe^2 * T1 * (T1 + 2m1) * T4) * q_int
end


Base.@kwdef mutable struct
XSecDMNucleusConstant{T <: Number, U <: Number} <: XSecElastic
    sigma0::T = 1e-30 * units.cm2
    mchi::U = 0.1 * units.GeV
    nucleus::String = "H"
    ff::String = "helm"
    lambda_dict::Dict = Dict(
        "Proton" => 0.77 * units.GeV^2,
        "He4" => 0.41 * units.GeV^2,
        "H" => 0.77 * units.GeV^2,
        "He" => 0.41 * units.GeV^2,
    )
end
settarget!(xsec::XSecDMNucleusConstant, nuc) = xsec.nucleus = nuc
gettarget(xsec::XSecDMNucleusConstant) = xsec.nucleus
function dxsecdT4(xsec::XSecDMNucleusConstant, T4, T1, m1, m2)
    return dmnucleus_xsec_eff(xsec) * formfactor(xsec, 2m2*T4)^2 / T4_max(xsec, T1, m1, m2)
end
function dmnucleus_xsec_eff(xsec::XSecDMNucleusConstant)
    return (xsec0(xsec) * NUCLEUS_A[xsec.nucleus]^2
        * (reduce_m(dmmass(xsec), NUCLEUS_MASS[xsec.nucleus]) / reduce_m(dmmass(xsec), PROTON_MASS))^2)
end
function formfactor(xsec::XSecDMNucleusConstant, Q2)
    if xsec.ff == "const"
        return one(Q2)
    end
    A = NUCLEUS_A[xsec.nucleus]
    if xsec.ff == "helm" && A > 7
        return ff_helm(Q2, A)
    end
    return 1 / (1 + Q2 / xsec.lambda_dict[xsec.nucleus]^2)^2
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
function recoil_spectrum(
    xsec::XSec, T4, flux1, T1max, m1, m2;
    attenuation=nothing, kwargs...
)
    T1min = T1_min(xsec, T4, m1, m2)
    T1min < T1max || return 0

    if !isnothing(attenuation)
        # ekin_in_min = attenuation.T0Tz(ekin_in_min)
        ekin_in_min0 = attenuation.T0Tz(ekin_in_min)
        ekin_in_min = ekin_in_min0
    end

    Eunit = oneunit(T1min)
    if isnothing(attenuation)
        rate_res = quad(log(T1min / Eunit), log(T1max / Eunit); kwargs...) do logT1
            T1 = exp(logT1) * Eunit
            dxsecdT4(xsec, T4, T1, m1, m2) * flux1(T1) * T1
        end
    else
        rate_res = quad(log(T1min / Eunit), log(T1max / Eunit); kwargs...) do logT1
            T1 = exp(logT1) * Eunit
            T1z = attenuation(T1)
            dxsecdT4(xsec, T4, T1z, m1, m2) * flux1(T1) * T1
        end
    end

    return rate_res[1] / m2
end
function recoil_spectrum(
    xsec::XSec, T4::AbstractVector, T1::AbstractVector, flux1::AbstractVector, m1, m2;
    attenuation=nothing, kwargs...
)
    flux_in = loginterpolator(T1, flux1)
    return recoil_spectrum.(xsec::XSec, T4, Ref(flux_in), maximum(T1), m1, m2;
                            attenuation=attenuation, kwargs...)
end
function recoil_spectrum(
    xsec::XSecDMElectronBound, Tr, fluxchi, mchi, me, qmax, Tchimax; unit=1, kwargs...
)
    Eunit = oneunit(Tr)
    ΔE = Tr + xsec.Eb
    logqmax = log(qmax / Eunit)
    logqmin = log(ΔE / Eunit) + sqrt(eps(ΔE / Eunit))
    logqmin < logqmax || throw(DomainError("logqmin = $logqmin < logqmax = $logqmax"))
    logTchimin(logq) = log(T1_min_q(exp(logq) * Eunit, ΔE, mchi) / Eunit)
    logtchimax = log(Tchimax / Eunit)
    logTchimax(_) = logtchimax

    ke = sqrt(Tr * (Tr + 2me))
    μχe = reduce_m(mchi, me)

    # Cao et al.
    # factor = xsec0(xsec.freexsec) / (8 * Tr * μχe^2 * me)
    # unitzero = zero(unit / oneunit(factor))
    # rate_res::typeof(unitzero), _ = dblquad(
    #     logqmin, logqmax, logTchimin, logTchimax; kwargs...
    # ) do logTchi, logq
    #     logTchi < logtchimax || return unitzero
    #     Tchi = exp(logTchi) * Eunit
    #     q = exp(logq) * Eunit
    #     Q2 = q^2 - ΔE^2
    #     pchi = sqrt(Tchi * (Tchi + 2mchi))
    #     return (q * dm_formfactor(xsec.freexsec, Q2, Tchi)
    #         * xsec.ionff(ke, q)
    #         * mchi^2 / (pchi*(Tchi+mchi))
    #         * fluxchi(Tchi)
    #         * Tchi * q)
    # end

    factor = xsec0(xsec.freexsec) * mchi^2 / (8 * μχe^2 * Tr * me)
    unitzero = zero(unit / oneunit(factor))
    rate_res::typeof(unitzero), _ = dblquad(
        logqmin, logqmax, logTchimin, logTchimax; kwargs...
    ) do logTchi, logq
        logTchi < logtchimax || return unitzero
        Tchi = exp(logTchi) * Eunit
        q = exp(logq) * Eunit
        Q2 = q^2 - ΔE^2
        pchi2 = Tchi * (Tchi + 2mchi)
        return (q * dm_formfactor(xsec.freexsec, Q2, Tchi) * xsec.ionff(ke, q) / pchi2
                * fluxchi(Tchi) * Tchi * q)
    end

    return rate_res * factor
end


@doc raw"""
    ff_helm(Q2, A; s=units.fm)

Helm form factor:

```math
    F_\text{Helm}(Q^2)
=
    \frac{3 j_1(qR_1)}{qR_1} e^{-Q^2 s^2 / 2},
```

where ``j_1`` is the spherical Bessel function of the first kind, ``q = \sqrt{Q^2}``,
``R_A = 1.2 A^{1/3}`` fm, ``R_1 = \sqrt{R_A^2 - 5s^2}``, and ``s = 1`` fm.

# Arguments

- `Q2`: Momentum transfer squared.
- `A`: Number of nucleons.

# Keywords

- `s=units.fm`: The length scale, default to 1 fm.

# References

- Lewin, J.D.D., Smith, P.F.F., 1996. Review of mathematics, numerical factors, and
    corrections for dark matter experiments based on elastic nuclear recoil.
    [Astroparticle Physics 6, 87–112.](https://doi.org/10.1016/S0927-6505(96)00047-3)
- Engel, J., 1991. Nuclear form factors for the scattering of weakly interacting massive
    particles. [Physics Letters B 264, 114–119.]
    (https://doi.org/10.1016/0370-2693(91)90712-Y)

"""
function ff_helm(Q2, A; s=units.fm)
    # Q2 > 2eps(Q2) || return one(Q2)
    q = sqrt(Q2)  # GeV
    RA = 1.2 * A^(1/3) * s
    R1 = sqrt(RA^2 - 5s^2)
    q_R1 = q * R1
    # return 3 * sphericalbesselj(1, q_R1) / (q_R1) * exp(-Q2 * s^2 / 2)
    return _ff_Helm_core(q_R1/oneunit(q_R1)) * exp(-Q2 * s^2 / 2)
end

function _ff_Helm_core(x)
    if x >= 1e-1
        return 3.0 * (sin(x) - x * cos(x)) / x^3
    elseif x > 1e-2
        x2 = x * x;
        return 1.0 + x2 * (
                -0.1 + x2 * (
                 3.5714285714285714286e-3 + x2 * (
                  -6.6137566137566137566e-5 + x2 * 7.5156325156325156325e-7)))
    elseif x > 1e-4
        x2 = x * x;
        return 1.0 + x2 * (-0.1 + x2 * 3.5714285714285714e-3);
    elseif x > 1e-8
        return 1.0 - 0.1 * x * x;
    end
    return one(x);
end

# WARNING temporary using
Kinematics.T4_max(::XSecDMElectronVectorMediator, T1, m1, m2) = T4_max(T1, m1, m2)
Kinematics.T4_max(::XSecDMElectronBound, T1, m1, m2) = T4_max(T1, m1, m2)
