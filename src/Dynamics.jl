using MultiQuad: dblquad, quad
# using SpecialFunctions: sphericalbesselj

abstract type XSec end
abstract type XSecElastic <: XSec end
abstract type XSecDMElectronElastic <: XSecElastic end

export totalxsec_analytic, average_energy_loss, xsec_moment

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

Kinematics.T4_max(::XSec, T1, m1, m2) = error("unimplemented")
Kinematics.T4_max(::XSecElastic, T1, m1, m2) = T4_max(T1, m1, m2)

dm_formfactor(::XSec, Q2) = error("unsupproted form factor")
dm_formfactor(::XSec, Q2, T1) = error("unsupproted form factor")


"""
    dxsecdT4(xsec::XSec, T4, T1, m1, m2, args...; kwargs...)

Differential cross section ``\\frac{dσ}{dT_4}(T_4, T_1, m_1, m_2)``.

# Arguments

- `T4`: Recoil energy of the target particle 2.
- `T1`: Kinetic energy of the incident particle 1.
- `m1`: Mass of the incident particle 1.
- `m2`: Mass of the target particle 2.

"""
dxsecdT4(::XSec, T4, T1, m1, m2, args...; kwargs...) = error("unimplemented")

function xsec_moment(
    xsec::XSec, order::Int, T1, m1, m2, args...;
    Trcut=zero(units.eV), kwargs...
)
    Tr_max = T4_max(xsec, T1, m1, m2)
    Trcut < Tr_max || return zero(units.cm2)
    res, _ = quad(Trcut, Tr_max; kwargs...) do T4
       return T4^order * dxsecdT4(xsec, T4, T1, m1, m2, args...; kwargs...)
    end

    return res
end

function totalxsec(xsec::XSec, T1, m1, m2, args...; Trcut=zero(units.eV), kwargs...)
    Tr_max = T4_max(xsec, T1, m1, m2)
    Trcut < Tr_max || return zero(units.cm2)
    # res, err = quad(Trcut, Tr_max; kwargs...) do T4
    #     return dxsecdT4(xsec, T4, T1, m1, m2, args...; kwargs...)
    # end
    res, _ = quad(log(Trcut), log(Tr_max); kwargs...) do lnT4
        T4 = exp(lnT4)
        return dxsecdT4(xsec, T4, T1, m1, m2, args...; kwargs...) * T4
    end
    return res
end

function average_energy_loss(
    xsec::XSec, T1, m1, m2, args...;
    Trcut=zero(units.eV), kwargs...
)
    return (xsec_moment(xsec, 1, T1, m1, m2, args...; Trcut=Trcut, kwargs...)
            / totalxsec(xsec, T1, m1, m2, args...; Trcut=Trcut, kwargs...))
end


Base.@kwdef mutable struct XSecVectorMediator{T <: Number, U <: Number} <: XSecElastic
    sigma0::T = 1e-30 * units.cm2
    mchi::U = 1 * units.keV
    mmed::U = 1 * units.keV
    mt::U = ELECTRON_MASS
    qref::U = ELECTRON_MASS / 137
    mediator_limit::String = ""
end
function XSecVectorMediator(sigma0, mchi, mmed, mt, qref, limit)
    return XSecVectorMediator(sigma0, promote(mchi, mmed, mt, qref)..., limit)
end
function dxsecdT4(xsec::XSecVectorMediator, T4, T1, m1, m2, args...; kwargs...)
    μχt = reduce_m(xsec.mchi, xsec.mt)
    if xsec.mediator_limit == "light"
        return (xsec.sigma0 * xsec.qref^4 / μχt^2
                * (2*m2*(m1 + T1)^2 - T4*((m1 + m2)^2 + 2*m2*T1) + m2*T4*T4)
                / (4*T1*(2*m1 + T1)*(2*m2*T4)^2))
    elseif xsec.mediator_limit == "heavy"
        return (xsec.sigma0 / μχt^2
                * (2*m2*(m1 + T1)^2 - T4*((m1 + m2)^2 + 2*m2*T1) + m2*T4*T4)
                / (4*T1*(2*m1 + T1)))
    elseif isempty(xsec.mediator_limit)
        return (xsec.sigma0 * (xsec.qref^2 + xsec.mmed^2)^2 / μχt^2
                * (2*m2*(m1 + T1)^2 - T4*((m1 + m2)^2 + 2*m2*T1) + m2*T4*T4)
                / (4*T1*(2*m1 + T1)*(2*m2*T4+xsec.mmed*xsec.mmed)^2))
    end
    throw(KeyError("Unknown mediator_limit: $(xsec.mediator_limit)"))
end
""" A Non-relativistic approximation for direct detection. """
function dm_formfactor(xsec::XSecVectorMediator, Q2, Tchi)
    # return ((xsec.qref^2 + xsec.mmed^2)^2 / (Q2 + xsec.mmed^2)^2
    #     * (8me^2 * (Tchi + xsec.mchi)^2 - 2mt*Q2*((xsec.mchi + mt)^2 + 2mt * Tchi) + Q2^2)
    #     / (8me^2 * xsec.mchi^2))
    # non-relativistic approximation
    if xsec.mediator_limit == "light"
        return xsec.qref^4 / Q2^2 * (Tchi + xsec.mchi)^2 / xsec.mchi^2
    elseif xsec.mediator_limit == "heavy"
        return (Tchi + xsec.mchi)^2 / xsec.mchi^2
    end
    return ((xsec.qref^2 + xsec.mmed^2)^2 / (Q2 + xsec.mmed^2)^2
        * (Tchi + xsec.mchi)^2 / xsec.mchi^2)
end
function totalxsec_analytic(
    xsec::XSecVectorMediator, T1, m1, m2, args...; kwargs...)
    m2 == xsec.mt || throw(KeyError("m2 = $m2 != xsec.mt = $(xsec.mt)"))
    mt = m2
    Trmax = T4_max(xsec, T1, m1, mt)
    return (xsec0(xsec) / Trmax * (m1 + mt)^2 / (2mt * T1 + (m1 + mt)^2)
            * FDMintegral(xsec, T1, Trmax))
end
function FDMintegral(xsec::XSecVectorMediator, T, Tp)
    mchi = dmmass(xsec)
    mt = xsec.mt
    mchi = dmmass(xsec)
    A = 2mt * (mchi + T)^2
    B = -(2mt * T + (mchi + mt)^2)
    C = mt
    fac = 2mt * mchi^2
    if xsec.mediator_limit == "heavy"
        return Tp * (A + Tp * (B / 2 + Tp * C / 3)) / fac
    end
    mmed = mediatormass(xsec)
    mmed > 0 || error("mediator mass <= 0")
    mr = mmed^2 / 2mt
    return (mr + xsec.qref^2 / (2 * mt))^2 / fac * (
            (A / mr - B + C * (Tp + 2mr)) * Tp / (Tp + mr)
             + (2mr * C - B) * log(mr / (Tp + mr)))
    # return (mr + xsec.qref^2 / (2 * mt))^2 / fac * (
    #     A * Tp / (mr * (Tp + mr)) + B * (-Tp / (Tp + mr) + log((Tp + mr) / mr))
    #     + C * (Tp * (Tp + 2mr) / (Tp + mr) + 2mr * log(mr / (Tp + mr))))
end

_f0(t, r) = -1 / (t + r)
_f1(t, r) = r / (t + r) + log(t + r)
_f2(t, r) = -r^2  / (t + r) - 2r * log(t + r) + t


Base.@kwdef mutable struct
XSecScalarMediator{T <: Number, U <: Number, V <: Number} <: XSecDMElectronElastic
    sigma0::T = 1e-30 * units.cm2
    mchi::U = 0.1 * units.GeV
    mmed::U = 1e-6 * units.GeV
    alpha::V = 1/137
end
"""
    dxsecdT4(xsec::XSecScalarMediator, T4, T1, m1, m2, args...; kwargs...)

Scalar mediated DM-electron scattering cross section.

"""
function dxsecdT4(xsec::XSecScalarMediator, T4, T1, m1, m2, args...; kwargs...)
    mt = ELECTRON_MASS
    μχt = reduce_m(xsec.mchi, mt)
    Q2 = 2*m2*T4
    return (xsec.sigma0 * (xsec.mmed^2+(xsec.alpha * mt)^2)^2 / (32*m2*μχt^2*T1*(T1 + 2*m1))
            * (4*m1^2 + Q2)*(4*m2^2 + Q2)/(xsec.mmed^2 + Q2)^2)
end


struct IronizationFormFactor
    ks::AbstractVector
    qs::AbstractVector
    formfactor::Function
end


Base.@kwdef mutable struct XSecDMElectronBound{X<:XSec, F<:Function, T<:Number} <: XSec
    freexsec::X = XSecVectorMediator()
    ionff::F = (_, _) -> 1
    Eb::T = 25.7units.eV  # Xe 5s shell
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


function Kinematics.T1_min(xsec::XSecDMElectronBound, T4, m1, m2)
    # FIXME
    if m1 == dmmass(xsec)
        return T4 + xsec.Eb
    end
    # ΔE = T4 + xsec.Eb
    # p1min = p1_min_q(q, ΔE, m1)
    return T1_min(T4, m1, m2)
end
""" Warning: using non-relativistic approximation """
function dxsecdT4(xsec::XSecDMElectronBound, T4, T1, m1, m2, args...; kwargs...)
    m1 == dmmass(xsec) || error("Only DM scattering off bound electron is implemented")
    ΔE = T4 + xsec.Eb
    T1 > ΔE || return zero(xsec0(xsec)/oneunit(T4))
    ke = sqrt(2 * m2 * T4)
    p1 = sqrt(T1 * (T1 + 2 * m1))
    T3 = T1 - ΔE
    p3 = sqrt(T3 * (T3 + 2 * m1))
    qmin = max(p1 - p3, ΔE)
    qmax = p1 + p3
    μχt = reduce_m(dmmass(xsec), xsec.mt)
    q_int::typeof(T4), _ = quad(log(qmin), log(qmax); kwargs...) do logq
        q = exp(logq)
        Q2 = q^2 - ΔE^2
        Eχ = T1 + m1
        #  return q^2 * dm_formfactor(xsec.freexsec, Q2, T1) * xsec.ionff(ke, q)
        return q^2 * xsec.ionff(ke, q) * (4*Eχ * (Eχ - ΔE) - Q2)  # from Sheng Jie
    end
    #  return xsec0(xsec) * m1^2 / (8 * μχt^2 * T1 * (T1 + 2m1) * T4) * q_int
    return xsec0(xsec) / (32 * μχt^2 * T1 * (T1 + 2*m1) * T4) * q_int
end


Base.@kwdef mutable struct
XSecDMNucleusConstant{T <: Number, U <: Number, U2 <: Number} <: XSecElastic
    sigma0::T = 1e-30 * units.cm2
    mchi::U = 0.1 * units.GeV
    ff::String = "helm"
    lambda_dict::Dict{String, U2} = Dict(
        "Proton" => 0.77 * units.GeV^2,
        "He4" => 0.41 * units.GeV^2,
        "H" => 0.77 * units.GeV^2,
        "He" => 0.41 * units.GeV^2,
    )
end
function dxsecdT4(xsec::XSecDMNucleusConstant, T4, T1, m1, m2, target; kwargs...)
    return (dmnucleus_xsec_eff(xsec, target) * formfactor(xsec, 2m2 * T4, target)^2
            / T4_max(xsec, T1, m1, m2))
end
function dmnucleus_xsec_eff(xsec::XSecDMNucleusConstant, target)
    return (xsec0(xsec) * particle_A(target)^2
            * (reduce_m(dmmass(xsec), particle_mass(target))
               / reduce_m(dmmass(xsec), PROTON_MASS))^2)
end
function formfactor(xsec::XSecDMNucleusConstant, Q2, target)
    if xsec.ff == "const"
        return one(Q2/Q2)
    end
    A = particle_A(target)
    if xsec.ff == "helm" && A > 7
        return ff_helm(Q2, A)
    end
    return 1 / (1 + Q2 / xsec.lambda_dict[target]^2)^2
end


"""
    recoil_spectrum(xsec::XSec, T4, flux1, T1max, m1, m2, args...; attenuation=nothing)

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
    xsec::XSec, T4, flux1, T1max, m1, m2, args...;
    attenuation=nothing, kwargs...
)
    T1min = T1_min(xsec, T4, m1, m2)
    T1min < T1max || return 0

    if !isnothing(attenuation)
        # FIXME
        # ekin_in_min = attenuation.T0Tz(ekin_in_min)
        ekin_in_min0 = attenuation.T0Tz(ekin_in_min)
        ekin_in_min = ekin_in_min0
    end

    Eunit = oneunit(T1min)
    if isnothing(attenuation)
        rate_res = quad(log(T1min / Eunit), log(T1max / Eunit); kwargs...) do logT1
            T1 = exp(logT1) * Eunit
            dxsecdT4(xsec, T4, T1, m1, m2, args...; kwargs...) * flux1(T1) * T1
        end
    else
        rate_res = quad(log(T1min / Eunit), log(T1max / Eunit); kwargs...) do logT1
            T1 = exp(logT1) * Eunit
            T1z = attenuation(T1)
            dxsecdT4(xsec, T4, T1z, m1, m2, args...; kwargs...) * flux1(T1) * T1
        end
    end

    return rate_res[1] / m2
end
function recoil_spectrum(
    xsec::XSec, T4::AbstractVector, T1::AbstractVector, flux1::AbstractVector, m1, m2,
    args...; attenuation=nothing, kwargs...
)
    flux_in = loginterpolator(T1, flux1)
    return recoil_spectrum.(xsec::XSec, T4, Ref(flux_in), maximum(T1), m1, m2, args...;
                            attenuation=attenuation, kwargs...)
end
#=
function recoil_spectrum(
    xsec::XSecDMElectronBound, Tr::Number, fluxchi, mchi, me, qmax, Tchimax;
    unit=1, kwargs...
)
    Eunit = oneunit(Tr)
    ΔE = Tr + xsec.Eb
    qmax = min(qmax, sqrt(2me * Tr) + me - xsec.Eb)
    logqmax = log(qmax / Eunit)
    logqmin = log(ΔE / Eunit) + sqrt(eps(ΔE / Eunit))
    logqmin < logqmax || throw(DomainError("logqmin = $logqmin < logqmax = $logqmax"))
    logTchimin(logq) = log(T1_min_q(exp(logq) * Eunit, ΔE, mchi) / Eunit)
    logtchimax = log(Tchimax / Eunit)
    logTchimax(_) = logtchimax

    ke = sqrt(Tr * (Tr + 2me))
    μχt = reduce_m(mchi, me)

    # Cao et al.
    # factor = xsec0(xsec.freexsec) / (8 * Tr * μχt^2 * me)
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

    # Our result
    factor = xsec0(xsec.freexsec) * mchi^2 / (8 * μχt^2 * Tr * me)
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
function recoil_spectrum(
    xsec::XSecDMElectronBound, Tr::AbstractVector, fluxchi, mchi, me, qmax, Tchimax;
    unit=1, kwargs...
)
    dRdE = similar(Tr, typeof(unit))
    Threads.@threads for i in eachindex(Tr)
        dRdE[i] = recoil_spectrum(
            xsec, Tr[i], fluxchi, mchi, me, qmax, Tchimax; unit=unit, kwargs...
        )
    end
    return dRdE
end
=#


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
Kinematics.T4_max(::XSecVectorMediator, T1, m1, m2) = T4_max(T1, m1, m2)
Kinematics.T4_max(::XSecDMElectronBound, T1, m1, m2) = T4_max(T1, m1, m2)
