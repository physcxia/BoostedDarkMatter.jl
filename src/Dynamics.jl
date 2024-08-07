using MultiQuad: quad
# using SpecialFunctions: sphericalbesselj

abstract type XSec end
abstract type XSecElastic <: XSec end

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
    res, err = quad(Trcut, Tr_max; kwargs...) do T4
       return dxsecdT4(xsec, T4, T1, m1, m2, args...; kwargs...)
    end
    # res, _ = quad(log(Trcut), log(Tr_max); kwargs...) do lnT4
    #     T4 = exp(lnT4)
    #     return dxsecdT4(xsec, T4, T1, m1, m2, args...; kwargs...) * T4
    # end
    return res
end

function average_energy_loss(
    xsec::XSec, T1, m1, m2, args...;
    Trcut=zero(units.eV), kwargs...
)
    return (xsec_moment(xsec, 1, T1, m1, m2, args...; Trcut=Trcut, kwargs...)
            / totalxsec(xsec, T1, m1, m2, args...; Trcut=Trcut, kwargs...))
end


function check_target_mass(xsec::XSec, m1, m2)
    if m1 == dmmass(xsec)
        mt = m2
    elseif m2 == dmmass(xsec)
        mt = m1
    else
        error("DM mass changed.")
    end
    return mt
end


@kwdef mutable struct
XSecVectorMediator{T <: Number, U <: Number, V <: Number} <: XSecElastic
    sigma0::T = 1e-30 * units.cm2
    mchi::U = 1 * units.keV
    mmed::U = 1 * units.keV
    qref::V = ELECTRON_MASS / 137
    mediator_limit::String = ""
end
function XSecVectorMediator(sigma0, mchi, mmed, qref, limit)
    return XSecVectorMediator(sigma0, promote(mchi, mmed)..., qref, limit)
end
function dxsecdT4(xsec::XSecVectorMediator, T4, T1, m1, m2, args...; kwargs...)
    T4 <= T4_max(xsec, T1, m1, m2) || return zero(xsec0(xsec)/T4)
    μχt = reduce_m(xsec.mchi, check_target_mass(xsec, m1, m2))
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
function totalxsec_analytic(
    xsec::XSecVectorMediator, T1, m1, m2, args...; kwargs...)
    m1 == dmmass(xsec) || throw(KeyError("m1 must be DM mass"))
    mt = check_target_mass(xsec, m1, m2)
    Trmax = T4_max(xsec, T1, m1, mt)
    return (xsec0(xsec) / Trmax * (m1 + mt)^2 / (2mt * T1 + (m1 + mt)^2)
            * FDMintegral(xsec, T1, Trmax, mt))
end
function FDMintegral(xsec::XSecVectorMediator, T, Tp, mt)
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


@kwdef mutable struct
XSecAxialVectorMediator{T <: Number, U <: Number, V <: Number} <: XSecElastic
    sigma0::T = 1e-30 * units.cm2
    mchi::U = 1 * units.keV
    mmed::U = 1 * units.keV
    qref::V = ELECTRON_MASS / 137
    mediator_limit::String = ""
end
function XSecAxialVectorMediator(sigma0, mchi, mmed, qref, limit)
    return XSecAxialVectorMediator(sigma0, promote(mchi, mmed)..., qref, limit)
end
function dxsecdT4(xsec::XSecAxialVectorMediator, T4, T1, m1, m2, args...; kwargs...)
    μχt = reduce_m(xsec.mchi, check_target_mass(xsec, m1, m2))
    if xsec.mediator_limit == "light"
        return (xsec.sigma0 * xsec.qref^4 / μχt^2
                * (2*m2*(m1 + T1)^2 - T4*((m1 + m2)^2 + 2*m2*T1) + m2*T4*T4)
                / (4*T1*(2*m1 + T1)*(2*m2*T4)^2))
    elseif xsec.mediator_limit == "heavy"
        return (xsec.sigma0 / μχt^2
                * (2*m2*((m1 + T1)^2 + 2*m1^2) + T4*((m1 - m2)^2 - 2*m2*T1) + m2*T4*T4)
                / (4*T1*(2*m1 + T1)))
    elseif isempty(xsec.mediator_limit)
        return (xsec0(xsec) * (xsec.qref^2 + xsec.mmed^2)^2 / μχt^2
                * (2*m2*((m1 + T1)^2 + 2*m1^2) + T4*((m1 - m2)^2 - 2*m2*T1) + m2*T4*T4)
                / (4*T1*(2*m1 + T1)*(2*m2*T4+xsec.mmed*xsec.mmed)^2))
    end
    throw(KeyError("Unknown mediator_limit: $(xsec.mediator_limit)"))
end


@kwdef mutable struct
XSecScalarMediator{T <: Number, U <: Number, V <: Number} <: XSecElastic
    sigma0::T = 1e-30 * units.cm2
    mchi::U = 0.1 * units.GeV
    mmed::U = 1e-6 * units.GeV
    qref::V = ELECTRON_MASS / 137
    mediator_limit::String = ""
end
function XSecScalarMediator(sigma0, mchi, mmed, qref, limit)
    return XSecScalarMediator(sigma0, promote(mchi, mmed)..., qref, limit)
end
"""
    dxsecdT4(xsec::XSecScalarMediator, T4, T1, m1, m2, args...; kwargs...)

Scalar mediated fermionic DM-fermion scattering cross section.

"""
function dxsecdT4(xsec::XSecScalarMediator, T4, T1, m1, m2, args...; kwargs...)
    μχt = reduce_m(xsec.mchi, check_target_mass(xsec, m1, m2))
    Q2 = 2 * m2 * T4
    if xsec.mediator_limit == "heavy"
        return xsec0(xsec) * (4*m1^2 + Q2)*(4*m2^2 + Q2) / (32*m2*μχt^2*T1*(T1 + 2*m1))
    elseif xsec.mediator_limit == "light"
        return (xsec0(xsec) * xsec.qref^4 / (32*m2*μχt^2*T1*(T1 + 2*m1))
                * (4*m1^2 + Q2)*(4*m2^2 + Q2) / Q2^2)
    elseif isempty(xsec.mediator_limit)
        return (xsec0(xsec) * (xsec.mmed^2+xsec.qref^2)^2 / (32*m2*μχt^2*T1*(T1 + 2*m1))
                * (4*m1^2 + Q2)*(4*m2^2 + Q2)/(xsec.mmed^2 + Q2)^2)
    end
    throw(KeyError("Unknown mediator_limit: $(xsec.mediator_limit)"))
end

@kwdef mutable struct
XSecPseudoScalarMediator{T <: Number, U <: Number, V <: Number} <: XSecElastic
    sigma0::T = 1e-30 * units.cm2
    mchi::U = 0.1 * units.GeV
    mmed::U = 1e-6 * units.GeV
    qref::V = ELECTRON_MASS / 137
    mediator_limit::String = ""
end
function XSecPseudoScalarMediator(sigma0, mchi, mmed, qref, limit)
    return XSecPseudoScalarMediator(sigma0, promote(mchi, mmed)..., qref, limit)
end
"""
    dxsecdT4(xsec::XSecPseudoScalarMediator, T4, T1, m1, m2, args...; kwargs...)

Pseudo scalar mediated fermionic DM-fermion scattering cross section.

"""
function dxsecdT4(xsec::XSecPseudoScalarMediator, T4, T1, m1, m2, args...; kwargs...)
    μχt = reduce_m(xsec.mchi, check_target_mass(xsec, m1, m2))
    Q2 = 2 * m2 * T4
    if xsec.mediator_limit == "heavy"
        return xsec0(xsec) * Q2^2 / (32*m2*μχt^2*T1*(T1 + 2*m1))
    elseif xsec.mediator_limit == "light"
        return xsec0(xsec) * xsec.qref^4 / (32*m2*μχt^2*T1*(T1 + 2*m1))
    elseif isempty(xsec.mediator_limit)
        return (xsec.sigma0 * (xsec.mmed^2+xsec.qref^2)^2 / (32*m2*μχt^2*T1*(T1 + 2*m1))
            * Q2^2 / (xsec.mmed^2 + Q2)^2)
    end
    throw(KeyError("Unknown mediator_limit: $(xsec.mediator_limit)"))
end


@kwdef mutable struct XSecDMElectronBound{X<:XSec, D1<:Dict, D2<:Dict} <: XSec
    freexsec::X = XSecVectorMediator()
    ionff::D1 = Dict("5s" => (_, _) -> 1)
    Eb::D2 = Dict("5s" => 25.7units.eV)  # Xe 5s shell
    shell::Vector{String} = ["5s"]
end
dmmass!(xsec::XSecDMElectronBound, mchi) = dmmass!(xsec.freexsec, mchi)
xsec0!(xsec::XSecDMElectronBound, sigma0) = xsec0!(xsec.freexsec, sigma0)
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
        return T4 + minimum([xsec.Eb[s] for s in xsec.shell])
    end
    return T1_min(T4, m1, m2)
end
# WARNING temporary using
# Kinematics.T4_max(::XSecDMElectronBound, T1, m1, m2) = T4_max(T1, m1, m2)
Kinematics.T4_max(::XSecDMElectronBound, T1, m1, m2) = T1

""" Warning: using non-relativistic approximation """
function dxsecdT4(xsec::XSecDMElectronBound, T4, T1, m1, m2, args...; kwargs...)
    m1 == dmmass(xsec) || error("Only DM scattering off bound electron is implemented")
    q_int = zero(T4)
    μχt = reduce_m(dmmass(xsec), check_target_mass(xsec, m1, m2))
    ke = sqrt(2 * m2 * T4)
    p1 = sqrt(T1 * (T1 + 2 * m1))
    for s in xsec.shell
        ΔE = T4 + xsec.Eb[s]
        T1 > ΔE || continue
        # T1 > ΔE || return zero(xsec0(xsec)/oneunit(T4))
        T3 = T1 - ΔE
        p3 = sqrt(T3 * (T3 + 2 * m1))
        qmin = max(p1 - p3, ΔE)
        qmax = p1 + p3
        q_int_s::typeof(T4), _ = quad(log(qmin), log(qmax); kwargs...) do logq
            q = exp(logq)
            Eχ = T1 + m1
            if xsec.freexsec isa XSecVectorMediator
                factor = (4*Eχ * (Eχ - ΔE) - (q^2 - ΔE^2)) / T4
            elseif xsec.freexsec isa XSecAxialVectorMediator
                factor = (2*Eχ * (Eχ - ΔE) - (q^2 - ΔE^2) / 2 - 2 * m1^2) / m2
            elseif xsec.freexsec isa XSecScalarMediator
                factor = (4 * m1^2 + q^2 - ΔE^2) / T4
            elseif xsec.freexsec isa XSecPseudoScalarMediator
                factor = ((q^2 - ΔE^2) / 2) / m2
            else
                error("unimplemented free xsec for bound electron scattering")
            end
            return q^2 * xsec.ionff[s](ke, q) * factor
        end
        q_int += q_int_s
    end
    return xsec0(xsec) / (32 * μχt^2 * p1^2) * q_int
end


@kwdef mutable struct
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
- `T1max`: Maximal energy of `flux1`.
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

@kwdef mutable struct XSecCombine{X1 <: XSec, X2 <: XSec, T <: Number} <: XSec
    xsec1::X1
    xsec2::X2
    Tbr::T
end
function dxsecdT4(xsec::XSecCombine, T4, T1, m1, m2, args...; kwargs...)
    if T4 < xsec.Tbr
        return dxsecdT4(xsec.xsec1, T4, T1, m1, m2, args...; kwargs...)
    end
    return dxsecdT4(xsec.xsec2, T4, T1, m1, m2, args...; kwargs...)
end
function Kinematics.T4_max(xsec::XSecCombine, T1, m1, m2)
    T4max1 = T4_max(xsec.xsec1, T1, m1, m2)
    T4max2 = T4_max(xsec.xsec2, T1, m1, m2)
    return max(T4max1, T4max2)
end
function dmmass(xsec::XSecCombine)
    m = dmmass(xsec.xsec1)
    m == dmmass(xsec.xsec2) || error("m1 != m2")
    return m
end
function xsec0(xsec::XSecCombine)
    σ₀ = xsec0(xsec.xsec1)
    σ₀ == xsec0(xsec.xsec2) || error("σ1 != σ2")
    return σ₀
end
