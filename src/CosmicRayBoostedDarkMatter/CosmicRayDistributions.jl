abstract type CRDistribution end
abstract type CRFlux end

Broadcast.broadcastable(cr::CRDistribution) = Ref(cr)
Broadcast.broadcastable(cr::CRFlux) = Ref(cr)

@enum EnergyType ETR ETEk ETEkn ETEt

function crflux(cr::CRFlux, coordinates; etype::EnergyType=ETEk)
    if etype == ETEk
        return crflux_Ekn(cr, coordinates)
    end
end

@doc raw"""
    crflux_Ekn(cr::CRFluxSBPLElectron, coordinates)

Intensity in kinetic energy per nucleon ``\frac{dI}{dE_\text{kn}}``
[m``^{-2}``s``^{-1}``sr``^{-1}``GeV^``{-1}``].

# Arguments

- `T`: kinetic energy of CR electron

"""
function crflux_Ekn(cr::CRFlux, coordinates::Tuple)
    _, __, Ekn = coordinates
    return crflux_Ekn(cr, Ekn)
end
function crflux_Ekn(cr::CRFlux, r, b, l, Ekn)
    return crflux_Ekn(cr, Ekn)
end

@doc raw"""
    CRFluxSBPLElectron{T <: Number, U <: Number, V <: Number} <: CRFlux

Smooth broken power law parameterization of local CR electron flux:

```math
    Φ_\text{CRE}
=
    Φ_0 \left[ 1 + \left( \frac{E_\text{br1}}{E} \right)^k \right]^{(γ_1 - γ_2)/k}
    \left( \frac{E}{E_\text{br}} \right)^{-γ_2}
    \left[ 1 + \left( \frac{E}{E_\text{br2}} \right)^k \right]^{(γ_3 - γ_2)/k}.
```

# Fields

- `phi0::T<:Number=5.02e-6`: ``Φ₀`` [m``^{-2}``sr``^{-1}``s``^{-1}``GeV``^{-1}``].
- `Ebr1::U<:Number=46.0`: ``E_\text{br1}`` [GeV].
- `Ebr2::U<:Number=987.8`: ``E_\text{br2}`` [GeV].
- `Ebr::U<:Number=300.0`: ``E_\text{br}`` [GeV].
- `k::V<:Number=10.0`: ``k``.
- `gamma1::V<:Number=3.24`: ``γ_1``.
- `gamma2::V<:Number=3.1`: ``γ_2``.
- `gamma3::V<:Number=3.89`: ``γ_3``.

# References

- Ding, Y. C., Li, N., Wei, C. C., Wu, Y. L., & Zhou, Y. F. (2021). Implications of a
    possible TeV break in the cosmic-ray electron and positron flux.
    [Physical Review D, 103(11), 115010.](https://doi.org/10.1103/PhysRevD.103.115010)

"""
Base.@kwdef struct CRFluxSBPLElectron{T <: Number, U <: Number, V <: Number} <: CRFlux
    phi0::T = 5.02e-6
    Ebr1::U = 46.0
    Ebr2::U = 987.8
    Ebr::U = 300.0
    k::V = 10.0
    gamma1::V = 3.24
    gamma2::V = 3.1
    gamma3::V = 3.89
end

function crflux_Ekn(cr::CRFluxSBPLElectron, T::Number)
    return (cr.phi0 *
            (1 + (cr.Ebr1 / T)^cr.k)^((cr.gamma1 - cr.gamma2) / cr.k)
            * (T / cr.Ebr)^-cr.gamma2
            * (1 + (T / cr.Ebr2)^cr.k)^((cr.gamma2 - cr.gamma3) / cr.k))
end


@doc raw"""
    CRFluxLISElectron{T <: Number, U <: Number} <: CRFlux

A parameterization of local CR electron flux:

```math
    F(T)
=
    \begin{cases}
    \frac{1.181×10^{11} T^{-12.061}}{1 + 4.307 × 10^8 T^{-9.269} + 3.125×10^8 T^{-10.697}},
    & T < 6.88 \text{GeV} \\
    995.598 T^{-3.505} + 4.423 T^{-2.620},
    & T ≥ 6.88 \text{GeV}
    \end{cases},
```

where the units are m``^{-2}``sr``^{-1}``s``^{-1}``GeV``^{-1}``. The parameterization is
valid from ``E_\text{cut1} = 2`` MeV up to ``E_\text{cut2} = 90`` GeV. Nevertheless, the two
cutoff energies are set as free parameters for convenience.

# Fields

- `Ecut1::T<:Number=0.002`: Lower cutoff kinetic energy [GeV].
- `Ebr::T<:Number=6.88`: Kinetic energy break [GeV].
- `Ecut2::T<:Number=90.0`: Upper cutoff kinetic energy [GeV].
- `unitflux::U<:Number=1.0`: Unit of the flux (for Unitful.jl support).

# References

- Boschini, M. J., et al. "HelMod in the works: from direct observations to the local
    interstellar spectrum of cosmic-ray electrons."
    [The Astrophysical Journal 854.2 (2018): 94.](https://doi.org/10.3847/1538-4357/aaa75e)

"""
Base.@kwdef struct CRFluxLISElectron{T <: Number, U <: Number} <: CRFlux
    Tcut1::T = 0.002
    Tbr::T = 6.88
    Tcut2::T = 90.0
    unitflux::U = 1.0
end
function crflux_Ekn(cr::CRFluxLISElectron, Ek::Number)
        if Ek < cr.Tcut1 || Ek > cr.Tcut2
            return zero(cr.unitflux)
        end
        T = Ek / oneunit(cr.Tbr)
        if Ek < cr.Tbr
            return cr.unitflux*1.181e11*T^-12.061/(1 + 4.307e8*T^-9.269 + 3.125e8*T^-10.697)
        end
        return cr.unitflux*(995.598 * T^-3.505 + 4.423 * T^-2.62)
end
