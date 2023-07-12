include("FluxConversion.jl")
using Interpolations: linear_interpolation

using HDF5: h5open

abstract type CRDistribution end
abstract type CRFlux end
abstract type CRFluxElectron <: CRFlux end

Broadcast.broadcastable(cr::CRDistribution) = Ref(cr)
Broadcast.broadcastable(cr::CRFlux) = Ref(cr)


@enum EnergyType ETR ETEk ETEkn ETEt

function crflux(cr::CRFlux, coordinates...; etype::EnergyType=ETEk, kwargs...)
    if etype == ETEk
        return crflux_Ek(cr, coordinates..., kwargs...)
    elseif etype == ETEkn
        return crflux_Ekn(cr, coordinates..., kwargs...)
    end
end

function (cr::CRFlux)(coordinates...; etype::EnergyType=ETEk, kwargs...)
    return crflux(cr, coordinates...; etype=etype, kwargs...)
end

crflux_Ek(cr::CRFluxElectron, T::Number) = crflux_Ekn(cr, T)
crflux_Ek(cr::CRFluxElectron, args...) = crflux_Ekn(cr, args...)
# crflux_Ek(cr::CRFluxElectron, r, b, l, Ekn) = crflux_Ekn(cr, Ekn)

@doc raw"""
    crflux_Ekn(cr::CRFluxSBPLElectron, coordinates::Tuple, args...; kwargs...)
    crflux_Ekn(cr::CRFlux, r, b, l, Ekn, args...; kwargs...)

Intensity in kinetic energy per nucleon ``\frac{dI}{dE_\text{kn}}``
[m``^{-2}``s``^{-1}``sr``^{-1}``GeV``^{-1}``].

# Arguments

- `T`: kinetic energy of CR electron

"""
function crflux_Ekn(
    cr::CRFlux,
    coordinates::Tuple{T, T, U},
    args...;
    kwargs...
) where {T <: Number, U <: Number}
    _, __, Ekn = coordinates
    return crflux_Ekn(cr, Ekn, args..., kwargs...)
end
function crflux_Ekn(cr::CRFlux, r, b, l, Ekn, args...; kwargs...)
    return crflux_Ekn(cr, Ekn, args...; kwargs...)
end

@doc raw"""
    CRFluxSBPLElectron{T <: Number, U <: Number, V <: Number} <: CRFluxElectron

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
Base.@kwdef struct CRFluxSBPLElectron{
    T <: Number, U <: Number, V <: Number
} <: CRFluxElectron
    phi0::T = 5.02e-6 * units.GeV^-1 * units.m^-2 * units.sec^-1
    Ebr1::U = 46.0 * units.GeV
    Ebr2::U = 987.8 * units.GeV
    Ebr::U = 300.0 * units.GeV
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
    CRFluxLISElectron{T <: Number, U <: Number} <: CRFluxElectron

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
Base.@kwdef struct CRFluxLISElectron{T <: Number, U <: Number} <: CRFluxElectron
    Tcut1::T = 0.002 * units.GeV
    Tbr::T = 6.88 * units.GeV
    Tcut2::T = 90.0 * units.GeV
    unitflux::U = units.GeV^-1 * units.m^-2 * units.sec^-1
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


@doc raw"""
    CRFluxLISHelMod2017{T <: Number, U <: Number} <: CRFlux

A parametrization of local interstellar CR nuclei flux:

```math
    F(R) × R^{2.7}
=
    \begin{cases}
    ∑_{i=0}^5 a_i R^i, & R ≤ 1 \text{GV}, \\
    b + \frac{c}{R} + \frac{d_1}{d_2 + R} + \frac{e_1}{e_2 + R} + \frac{f_1}{f_2 + R} + gR,
    & R ≥ 1 \text{GV}
    \end{cases},
```


where ``R`` is rigidity and the unit is [m``^{-2}`` s``^{-1}`` sr``^{-1}`` GV``^{-1}``].

# References

- M. J. Boschini et al., 2017. Solution of Heliospheric Propagation: Unveiling the Local
    Interstellar Spectra of Cosmic-ray Species.
    [Astrophys. J. 840, 115.](https://doi.org/10.3847/1538-4357/aa6e4f)

"""
Base.@kwdef struct CRFluxLISHelMod2017{T <: Number, U <: Number} <: CRFlux
    params::Dict{String, Dict{Char, Vector{Float64}}} = Dict(
        "Proton" => Dict(
            'a' => [94.1, -831., 0., 16700., -10200., 0.],
            'b' => [10800.],
            'c' => [8590.],
            'd' => [-4230000., 3190.],
            'e' => [274000., 17.4],
            'f' => [-39400., 0.464],
            'g' => [0.]
        ),
        "He" => Dict(
            'a' => [1.14, 0., -118., 578., 0., -87.],
            'b' => [3120.],
            'c' => [-5530.],
            'd' => [3370., 1.29],
            'e' => [134000., 88.5],
            'f' => [-1170000., 861.],
            'g' => [0.03]
        )
    )
    A::Dict{String, Int} = Dict("Proton" => 1, "He" => 4)
    Z::Dict{String, Int} = Dict("Proton" => 1, "He" => 2)
    Rbr::Dict{String, T} = Dict("Proton" => 1.0, "He" => 2.0)
    Rcut::T = 0.2
    unitflux::U = units.GV^-1 * units.m^-2 * units.sec^-1
end
function crflux_R(cr::CRFluxLISHelMod2017, R, particle="Proton")
    df = cr.params[particle]

    Rbr = cr.Rbr[particle]
    if R > Rbr
        return (df['b'][1]
                + df['c'][1] / R
                + df['d'][1] / (df['d'][2] + R)
                + df['e'][1] / (df['e'][2] + R)
                + df['f'][1] / (df['f'][2] + R)
                + df['g'][1] * R) * R^-2.7 * cr.unitflux
    elseif R > cr.Rcut
        return sum([a * R^(i-1) for (i, a) in enumerate(df['a'])]) * R^-2.7 * cr.unitflux
    end
    return zero(cr.unitflux)
end
function crflux_Ek(cr::CRFluxLISHelMod2017, Ek, particle="Proton")
    A = cr.A[particle]
    Z = cr.Z[particle]
    R = Ek_to_R(Ek, Z, A)
    return flux_R_to_Ek(R, crflux_R(cr, R, particle), Z, A)
end
function crflux_Ekn(cr::CRFluxLISHelMod2017, Ekn, particle="Proton")
    A = cr.A[particle]
    Z = cr.Z[particle]
    R = Ekn_to_R(Ekn, Z, A)
    return flux_R_to_Ekn(R, crflux_R(cr, R, particle), Z, A)
end

abstract type CRDGalprop <: CRDistribution end

function crflux_Ekn(cr::CRDGalprop, coordinates::Tuple; particle="Proton")
    return crflux_Ekn(cr, coordinates...; particle=particle)
end
function crflux_Ekn(cr::CRDGalprop, coordinates...; particle="Proton")
    if particle isa String
        if lowercase(particle) == "all"
            return sum([crflux(coordinates...) for (_, crflux) in cr.crflux])
        elseif lowercase(particle) == "proton"
            return sum([cr.crflux[p](coordinates...)
                        for p in ["Hydrogen_1", "secondary_protons"]
                        if p in keys(cr.crflux)])
        elseif lowercase(particle) == "electron"
            return sum([cr.crflux[p](coordinates...)
                        for p in ["primary_electrons", "secondary_electrons"]
                        if p in keys(cr.crflux)])
        elseif particle in keys(NUCLEUS_NAME)
            return sum([cr.crflux[p](coordinates...)
                        for p in _find_isotopes(cr, NUCLEUS_NAME[particle])])
        else
            return cr.crflux[particle](coordinates...)
        end
    elseif particle isa AbstractVector
        return sum([crflux_Ekn(cr, coordinates...; particle=p) for p in particle])
    end
    error("unknown particle: $particle")
end

function make_crflux_dict_galactic(
    cr::CRDGalprop, particles::Union{String, AbstractVector}, rsun=8.5*units.kpc
)
    particles = particles isa String ? [particles] : particles
    return Dict(
        (p => (r, b, l, Ek) -> crflux_Ek(cr, r, b, l, Ek; rsun=rsun, particle=p))
        for p in particles
    )
end
function make_crflux_dict_galactic(cr::CRDGalprop)
    particles = replace!(collect(keys(NUCLEUS_NAME)), "H" => "Proton")
    return make_crflux_dict_galactic(cr::CRDGalprop, particles)
end

function _find_isotopes(cr::CRDGalprop, particle)
    iso_names = [key for key in keys(cr.crflux) if occursin(particle, key)]
    length(iso_names) > 0 || error("$particle not found")
    return iso_names
end

struct CRDGalpropCylindrical{T <: Number, U <: Number, V <: Number} <: CRDGalprop
    r::Vector{Float64}
    z::Vector{Float64}
    Ekn::Vector{Float64}
    crflux::Dict{String, Function}
    Eunit::T
    coorunit::U
    fluxunit::V
end
function CRDGalpropCylindrical(
    crfilename::String;
    Eunit=units.GeV,
    coorunit=units.kpc,
    fluxunit=units.cm^-2 * units.sec^-1 * units.GeV^-1,
)
    fid = h5open(crfilename)
    r = read(fid, "r")
    z = read(fid, "z")
    Ekn = read(fid, "Ekin")
    lnEkn = log.(Ekn)
    crflux = Dict{String, Function}()
    for p_name in keys(fid)
        p_name != "r" && p_name != "z" && p_name != "Ekin" || continue
        dset = read(fid, p_name * "/flux")
        lndset = @. ifelse(dset <= 0, -Inf, log(dset))
        itp = linear_interpolation(
            (r, z, lnEkn),
            permutedims(lndset, (3, 2, 1)),
            extrapolation_bc=-Inf
        )

        function wrapper(r, z, Ekn)
            flux = exp(itp(r / coorunit, z / coorunit, log(Ekn / Eunit)))
            return isnan(flux) ? zero(fluxunit) : flux * fluxunit
        end

        crflux[p_name] = wrapper
    end
    close(fid)
    return CRDGalpropCylindrical(r, z, Ekn, crflux, Eunit, coorunit, fluxunit)
end
function crflux_Ek(
    crd::CRDGalpropCylindrical, r, b, l, Ek;
    rsun=8.5*units.kpc, particle="Proton"
)
    z = r * sin(b)
    rho = sqrt(rsun^2 + r * cos(b) * (r * cos(b) - 2 * rsun*cos(l)))
    A = lowercase(particle) == "electron" ? 1 : particle_A(particle)
    return crflux_Ekn(crd, rho, z, Ek / A; particle=particle) / A
end
