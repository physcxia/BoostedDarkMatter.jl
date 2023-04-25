include("CosmicRayDistributions.jl")

using .Kinematics: T1_min

using DarkMatterProfiles: DMProfile, DMPNFW, dmdensity_galactic
using MultiQuad: tplquad, quad


"""
    struct CRDM{T <: Number} <: BDM

Cosmic ray boosted dark matter (CRDM) model.

# Fields

- `kfactor::Dict`: Dictionary of the ``K`` factor as a function of total kinetic energy for
                   each CR species. See Ref.[3] for definition.
- `crdist::Dict`: Dictionary of CR distributions. The CR distributions have to be callable
                  with the signature (`r`, `b`, `l`, `Ek`), where (`r`, `b`, `l`) is
                  galactic coordinate and `Ek` is the kinetic energy of the CR particle.
- `Ekn_max::Dict`: Dictionary of the maximal *kinetic energy per nucleon* for CR species.
- `dmprofile::DMProfile`: DM profile.
- `xsec::XSec`: DM-CR scattering cross section.
- `selected_cr::Vector{String}`: The selected CR species for calculation.
- `rsun::T<:Number`: The distance between the Sun and the Galactic Center.

# References

- [1] Bringmann, T., Pospelov, M., 2019. Novel Direct Detection Constraints on Light Dark
      Matter. [Phys. Rev. Lett. 122, 171801.](https://doi.org/10.1103/PhysRevLett.122.171801)
- [2] Ema, Y., Sala, F., Sato, R., 2018. Light Dark Matter at Neutrino Experiments.
      [Phys. Rev. Lett. 122, 181802.](https://doi.org/10.1103/PhysRevLett.122.181802)
- [3] Xia, C., Xu, Y.-H., Zhou, Y.-F., 2022. Production and attenuation of cosmic-ray boosted
      dark matter. [JCAP 02, 028.](https://doi.org/10.1088/1475-7516/2022/02/028)

"""
Base.@kwdef mutable struct CRDM{T <: Number} <: BDM
    crdist::Dict = Dict("Electron" => CRFluxLISElectron())
    kfactor::Dict{String, Function} = Dict{String, Function}()
    Ekn_max::Dict = Dict((p => 100*units.GeV) for p in keys(crdist))
    dmprofile::DMProfile = DMPNFW(
        rho0=0.3*units.GeV/units.cm^3,
        rsun=8.5*units.kpc,
        rs=20*units.kpc,
        rcut=4e-10*units.kpc
    )
    xsec::XSec = XSecDMElectronVectorMediator()
    selected_cr::Vector{String} = collect(keys(crdist))
    rsun::T = 8.5 * units.kpc
end
function crdist!(bdm::CRDM, crdist::Dict)
    bdm.crdist = crdist
end
dmmass!(bdm::CRDM, mchi) = dmmass!(bdm.xsec, mchi)
xsec0!(bdm::CRDM, sigma0) = xsec0!(bdm.xsec, sigma0)
mediatormass!(bdm::CRDM, m) = mediatormass!(bdm.xsec, m)
dmmass(bdm::CRDM) = dmmass(bdm.xsec)
xsec0(bdm::CRDM) = xsec0(bdm.xsec)
mediatormass(bdm::CRDM) = mediatormass(bdm.xsec)
select_cr!(bdm::CRDM, crs) = bdm.selected_cr = crs

"""
    dmflux(bdm::CRDM, Tchi; Ekn_cutoff=nothing, fluxunit=1.0, options...)

Return the total flux of CRDM on the Earth.

Use the ``K`` factor and the local cosmic ray spectrum to calculate the flux of CRDM
produced by selected CR species.

# Arguments

- `bdm::CRDM`: Cosmic ray boosted dark matter model.
- `Tchi`: Kinetic energy of DM particle.

# Keywords

- `Ekn_cutoff`: Cutoff of CR kinetic energy per nucleon. If `nothing`,
                `bdm.Ekn_max` will be used (the default is `nothing`).
- `fluxunit=1.0`: The unit of the flux.
- `options...`: Other keyword arguments are passed to the
                [`quad`](https://github.com/aurelio-amerio/MultiQuad.jl#quad) integrator.

"""
function dmflux(bdm::CRDM, Tchi; Ekn_cutoff=nothing, fluxunit=1.0, options...)
    flux = 0 * fluxunit
    old_target = gettarget(bdm.xsec)
    for particle in bdm.selected_cr
        A = NUCLEUS_A[particle]
        if isnothing(Ekn_cutoff)
            cutoff = bdm.Ekn_max[particle] * A
        elseif Ekn_cutoff isa Number
            cutoff = Ekn_cutoff * A
        elseif Ekn_cutoff isa Dict
            cutoff = get(Ekn_cutoff, particle, bdm.Ekn_max[particle]) * A
        else
            throw(KeyError("Unknown Ekn_cutoff type."))
        end

        settarget!(bdm.xsec, particle)
        mi = particle_mass(particle)
        Ti_min = T1_min(bdm.xsec, Tchi, mi, dmmass(bdm))

        Ti_min < cutoff || continue

        flux_res = quad(Ti_min, cutoff, options...) do Ti
            ((kfactor(bdm, Ti, particle)
                * _los_integrand(bdm, zero(bdm.rsun), 0, 0, Ti, bdm.rsun, particle)
                * dxsecdT4(bdm.xsec, Tchi, Ti, mi, dmmass(bdm))
                / dmmass(bdm)) * 4π)  # 4 π total flux
        end
        flux += flux_res[1]
    end
    settarget!(bdm.xsec, old_target)
    return flux
end

"""
    kfactor!(bdm::CRDM, kfactor::Number, particles="All")
    kfactor!(bdm::CRDM, kfactor::Function, particles="All")

Set the angular averaged effective distance. Since the calculation of ``D_\\text{eff}`` is
time-consuming, it is better to manually set the Deff function calculated in advance.

# Arguments

- `bdm`: Deff function, given a single number means constant Deff.
particles : str or tuple of str, optional
    The names of particles to be set (the default is "All").

"""
function kfactor!(bdm::CRDM, kfactor::Number, particles="All")
    kfactor!(bdm, (_ -> kfactor), particles)
end
function kfactor!(bdm::CRDM, func::Function, particles="All")
    particles = _parse_particles(bdm, particles)
    for p in particles
        bdm.kfactor[p] = func
    end
end
function kfactor(bdm::CRDM, Ek, particle)
    bdm.kfactor[particle](Ek)
end
function kfactor(dmp::DMProfile, rmax, zmax=zero(rmax), rsun=8.5*units.kpc)
    local_value = dmdensity_galactic(dmp, zero(rmax), zero(rmax), zero(rmax), rsun=rsun)

    # upper limit of r as a function of r and b
    rfun = (iszero(zmax) ? (_, _) -> rmax : (l, b) -> _los_bound(b, l, rmax, zmax, rsun))

    kfactor_res = tplquad(
        0, 2π,
        b -> -π/2, b -> π/2,
        (l, b) -> 0, rfun
    ) do r, b, l
        dmdensity_galactic(dmp, r, b, l, rsun=rsun) * cos(b) / (4π * local_value)
    end

    return kfactor_res[1]
end

function kfactor!(
    bdm::CRDM,
    Ti::AbstractVector,
    rmax,
    zmax=zero(rmax),
    particles="all";
    kwargs...
)
    for p in _parse_particles(bdm, particles)
        # upper limit of r as a function of r and b
        rfun = (iszero(zmax) ?
                (_, _) -> rmax :
                (l, b) -> _los_bound(b, l, rmax, zmax, bdm.rsun))

        kfacs = similar(Ti)
        Threads.@threads for i in 1:length(Ti)
            local_value = _los_integrand(bdm, zero(bdm.rsun), 0, 0, Ti[i], bdm.rsun, p)
            kfacs[i], _ = tplquad(
                0.0, 2π, b -> -π/2, b -> π/2, (l, b) -> 0.0, rfun; kwargs...
            ) do r, b, l
                _los_integrand(bdm, r, b, l, Ti[i], bdm.rsun, p) * cos(b) / local_value
            end
            kfacs[i] /= 4π
        end
        bdm.kfactor[p] = loginterpolator(Ti, kfacs, method="xlog")
    end
    return nothing
end


"""Line of sight boundary. """
function _los_bound(b, l, rmax, zmax, rsun)
    rsun_cosl = rsun * cos(l)

    # boundary is on the z = +/- zmax planes
    if abs(tan(b)) > (zmax / (rsun_cosl + sqrt(rsun_cosl^2 + rmax^2 - rsun^2)))
        return zmax / abs(sin(b))
    end

    return (rsun_cosl + sqrt(rsun_cosl^2 + rmax^2 - rsun^2)) / cos(b)
end


"""Dark matter profile times cosmic ray flux. """
function _los_integrand(bdm::CRDM, r, b, l, Ek, rsun, particle)
    return (bdm.crdist[particle](r, b, l, Ek)
            * dmdensity_galactic(bdm.dmprofile, r, b, l, rsun=rsun))
end


function _parse_particles(bdm::CRDM, particles::String)
    if lowercase(particles) == "all"
        return [key for key in keys(bdm.crdist)]
    end
    return [particles]
end
_parse_particles(::CRDM, particles::Vector{String}) = particles
