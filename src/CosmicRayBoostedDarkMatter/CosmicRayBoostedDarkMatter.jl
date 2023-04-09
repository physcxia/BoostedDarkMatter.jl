include("CosmicRayDistributions.jl")

using .Kinematics: T1_min

using DarkMatterProfiles: DMProfile, DMPNFW, dmdensity_galactic
using MultiQuad: tplquad, quad
using NumericalTools: loginterpolator


Base.@kwdef mutable struct CRDM <: BDM
    kfactor::Dict = Dict("Electron" => _ -> 1)
    crdist::Dict = Dict("Electron" => CRFluxLISElectron())
    crmass::Dict = Dict("Electron" => ELECTRON_MASS)
    Ekn_max::Dict = Dict("Electron" => 100)
    dmprofile::DMProfile = DMPNFW()
    xsec::XSec = XSecDMElectronVectorMediator()
end
function crdist!(bdm::CRDM, crdist::Dict)
    bdm.crdist = crdist
end
dmmass!(bdm::CRDM, mchi) = dmmass!(bdm.xsec, mchi)
xsec0!(bdm::CRDM, sigma0) = xsec0!(bdm.xsec, sigma0)
mediatormass!(bdm::CRDM, m) = mediatormass!(bdm.xsec, m)
dmmass(bdm::CRDM) = dmmass(bdm.xsec)
xsec0(bdm::CRDM) = xsec0(bdm.xsec)

"""
    dmflux(bdm::CRDM, Tchi; Ekn_cutoff=nothing, fluxunit=1.0, options...)

Return the total flux of CRDM on the Earth.

Use angular averaged effective distance and local cosmic ray spectrum
to calculate the flux of CRDM produced by all particles set.

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
    for particle in keys(bdm.crdist)
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

        mi = bdm.crmass[particle]
        Ti_min = T1_min(bdm.xsec, Tchi, mi, dmmass(bdm))

        Ti_min < cutoff || continue

        flux_res = quad(Ti_min, cutoff, options...) do Ti
            ((kfactor(bdm, Ti, particle)
                * _los_integrand(bdm, zero(bdm.dmprofile.rsun), 0, 0, Ti, particle)
                * dxsecdT4(bdm.xsec, Tchi, Ti, mi, dmmass(bdm))
                / dmmass(bdm)) * 4π)  # 4 pi total flux
        end
        flux += flux_res[1]
    end
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
function kfactor(dmp::DMProfile, rmax, zmax=zero(rmax), rsun=8.5*oneunit(rmax))
    local_value = dmdensity_galactic(dmp, zero(rmax), zero(rmax), zero(rmax))

    # upper limit of r as a function of r and b
    rfun = (iszero(zmax) ? (_, _) -> rmax : (l, b) -> _los_bound(b, l, rmax, zmax, rsun))

    kfactor_res = tplquad(
        0, 2π,
        b -> -π/2, b -> π/2,
        (l, b) -> 0, rfun
    ) do (r, b, l)
        dmdensity_galactic(dmp, r, b, l, rsun=dmp.rsun) * cos(b) / (4π * local_value)
    end

    return kfactor_res[1]
end


"""Line of sight boundary. """
function _los_bound(b, l, rmax, zmax, rsun)
    rsun_cosl = rsun * cos(l)

    # boundary is on the z = +/- zmax planes
    on_plane = abs(tan(b)) > (zmax / (rsun_cosl + sqrt(rsun_cosl^2 + rmax^2 - rsun^2)))
    if on_plane
        return zmax / abs(sin(b))
    end

    return (rsun_cosl + sqrt(rsun_cosl^2 + rmax^2 - rsun^2)) / cos(b)
end


"""Dark matter profile times cosmic ray flux. """
function _los_integrand(bdm::CRDM, r, l, b, Ek, particle)
    A = NUCLEUS_A[particle]
    return (bdm.crdist[particle](r, l, b, Ek / A) / A
            * dmdensity_galactic(bdm.dmprofile, r, l, b, rsun=bdm.dmprofile.rsun))
end


function _parse_particles(bdm, particles::String)
    if lowercase(particles) == "all"
        return [key for key in keys(bdm.crdist)]
    end
    return [particles]
end
