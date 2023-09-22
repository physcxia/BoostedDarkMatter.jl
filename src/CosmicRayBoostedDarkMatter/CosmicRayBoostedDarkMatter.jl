include("CosmicRayDistributions.jl")

import JSON
using HDF5: h5open, create_group, create_dataset, write, read, attributes

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
                  with the signature (`r`, `b`, `l`, `Ek`), where (`r`, `b`, `l`) is the
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
    kfactor::Dict{String, Function} = Dict{String, Function}(
        p => (_ -> 1units.kpc) for p in keys(crdist)
    )
    Ekn_max::Dict = Dict((p => 1000*units.GeV) for p in keys(crdist))
    r_max::Dict = Dict((p => 20*units.kpc) for p in keys(crdist))
    z_max::Dict = Dict((p => zero(units.kpc)) for p in keys(crdist))
    dmprofile::DMProfile = DMPNFW(
        rho0=0.3*units.GeV/units.cm^3,
        rsun=8.5*units.kpc,
        rs=20*units.kpc,
        rcut=4e-10*units.kpc
    )
    xsec::XSec = XSecVectorMediator()
    selected_cr::Vector{String} = collect(keys(crdist))
    rsun::T = 8.5 * units.kpc
end
function CRDM(xsec::XSec, crd::CRDGalpropCylindrical, particles)
    crdist = make_crflux_dict_galactic(crd, particles)
    r_max = maximum(crd.r) * crd.coorunit
    z_max = maximum(crd.z) * crd.coorunit
    Ekn_max = maximum(crd.Ekn) * crd.Eunit
    return CRDM(
        crdist = crdist,
        Ekn_max = Dict((p => Ekn_max) for p in keys(crdist)),
        r_max = Dict((p => r_max) for p in keys(crdist)),
        z_max = Dict((p => z_max) for p in keys(crdist)),
        xsec = xsec,
    )
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
selectcr!(bdm::CRDM, crs) = bdm.selected_cr = _parse_particles(bdm, crs)

"""
    dmflux(bdm::CRDM, Tchi; Ekn_cutoff=nothing, fluxunit=1.0, kwargs...)

Return the total flux of CRDM on the Earth.

Use the ``K`` factor and the local cosmic ray spectrum to calculate the flux (per solid
angle) of CRDM produced by selected CR species.

# Arguments

- `bdm::CRDM`: Cosmic ray boosted dark matter model.
- `Tchi`: Kinetic energy of DM particle.

# Keywords

- `Ekn_cutoff`: Cutoff of CR kinetic energy per nucleon. If `nothing`,
                `bdm.Ekn_max` will be used (the default is `nothing`).
- `fluxunit=1.0`: The unit of the flux.
- `kwargs...`: Other keyword arguments are passed to the
               [`quad`](https://github.com/aurelio-amerio/MultiQuad.jl#quad) integrator.

"""
function dmflux(
    bdm::CRDM, Tchi::Number, angles...;
    Ekn_cutoff=nothing, fluxunit=1.0, kwargs...
)
    flux = zero(fluxunit)
    for particle in bdm.selected_cr

        cutoff = _determine_Ekn_cutoff(bdm, particle, Ekn_cutoff)

        mi = particle_mass(particle)
        mchi = dmmass(bdm)
        Ti_min = T1_min(bdm.xsec, Tchi, mi, mchi)

        Ti_min < cutoff || continue

        flux_res::typeof(flux), _ = quad(Ti_min, cutoff; atol=zero(flux), kwargs...) do Ti
            (kfactor(bdm, Ti, angles..., particle; kwargs...)
             * _los_integrand(bdm, zero(bdm.rsun), 0, 0, Ti, particle) / mchi
             * dxsecdT4(bdm.xsec, Tchi, Ti, mi, mchi, particle; kwargs...))
        end
        flux += flux_res
    end
    return flux
end
function _determine_Ekn_cutoff(bdm::CRDM, particle, Ekn_cutoff)
    if lowercase(particle) == "electron"
        return bdm.Ekn_max[particle]
    end
    A = particle_A(particle)
    if isnothing(Ekn_cutoff)
        cutoff = bdm.Ekn_max[particle] * A
    elseif Ekn_cutoff isa Number
        cutoff = Ekn_cutoff * A
    elseif Ekn_cutoff isa Dict
        cutoff = get(Ekn_cutoff, particle, bdm.Ekn_max[particle]) * A
    else
        throw(KeyError("Unknown Ekn_cutoff type."))
    end
    return cutoff
end

"""
    kfactor!(bdm::CRDM, kfactor::Number)
    kfactor!(bdm::CRDM, kfactor::Function)
    kfactor!(bdm::CRDM, jsonfile::String)

Set the angular averaged effective distance. Since the calculation of the ``K`` factor is
time-consuming, it is better to manually set it in advance.

The particles selected in the `bdm.selected_cr` field will be set.

# Arguments

- `bdm::CRDM`: For CRDM model.
- `kfactor`: K factor function, given a single number means constant K factor.

"""
function kfactor!(bdm::CRDM, kfactor::Number)
    kfactor!(bdm, (_ -> kfactor))
end
function kfactor!(bdm::CRDM, func::Function)
    for p in bdm.selected_cr
        bdm.kfactor[p] = func
    end
    return nothing
end
function kfactor!(bdm::CRDM, jsonfile::String)
    data = JSON.parsefile(jsonfile)
    for (p, kfac) in data
        bdm.kfactor[p] = loginterpolator(
            Vector{Float64}(kfac[1]) .* units.GeV,
            Vector{Float64}(kfac[2]) .* units.kpc;
            method="xlog"
        )
    end
    return nothing
end
"""
    kfactor!(bdm::CRDM, Ti::AbstractVector; save_to_file=nothing, kwargs...)

Calculate and set the angular averaged effective distance.

# Arguments

- `bdm::CRDM`: For CRDM model.
- `Ti`: Energy grid of CR kinetic energy.

# Keywords

- `save_to_file=nothing`: Switch whether save the results to a json file, default is no.
                          Giving a file name means save to the file.

"""
function kfactor!(bdm::CRDM, Ti::AbstractVector; save_to_file=nothing, kwargs...)
    for p in bdm.selected_cr
        # upper limit of r as a function of r and b
        rfun = (iszero(bdm.z_max[p]) ?
                (_, _) -> bdm.r_max[p] :
                (l, b) -> _los_bound(b, l, bdm.r_max[p], bdm.z_max[p], bdm.rsun))

        kfacs = similar(Ti)
        Threads.@threads for i in 1:length(Ti)
            local_value = _los_integrand(bdm, zero(bdm.rsun), 0, 0, Ti[i], p)
            kfacs[i], _ = tplquad(
                0.0, 2π, b -> -π/2, b -> π/2, (l, b) -> 0.0, rfun;
                atol=0, rtol=1e-3, kwargs...
            ) do r, b, l
                _los_integrand(bdm, r, b, l, Ti[i], p) * cos(b) / local_value
            end
            kfacs[i] /= 4π
        end
        bdm.kfactor[p] = loginterpolator(Ti, kfacs, method="xlog")
    end
    if !isnothing(save_to_file)
        jsondump_kfactor(bdm, save_to_file, Ti)
    end
    return nothing
end
function kfactor(bdm::CRDM, Ek, particle; kwargs...)
    # TODO handle the kwargs
    bdm.kfactor[particle](Ek)
end
function kfactor(bdm::CRDM, Ek, b, l, particle; kwargs...)
    rzero = zero(bdm.rsun)
    loc_value = _los_integrand(bdm, rzero, 0, 0, Ek, particle)
    loc_value > zero(loc_value) || return zero(loc_value/loc_value)
    los_max = _los_bound(b, l, bdm.r_max[particle], bdm.z_max[particle], bdm.rsun)
    losi, _ = quad(rzero, los_max; atol=rzero, kwargs...) do r
        _los_integrand(bdm, r, b, l, Ek, particle) / loc_value
    end
    return losi
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

function jsondump_kfactor(bdm::BDM, filename, Ti::AbstractVector)
    data = Dict(p => [Ti./units.GeV bdm.kfactor[p].(Ti)./units.kpc]
                for p in keys(bdm.kfactor))
    open(filename, "w") do f
        write(f, JSON.json(data, 2))
    end
    return nothing
end

function make_intensity_file(
    filename::String, mchis::AbstractVector, shape::Tuple{Int, Int, Int};
    dtype=Float64, kwargs...
)
    shape = reverse(shape)
    f = h5open(filename, "cw"; kwargs...)
    f["mchi"] = mchis ./ units.GeV
    attributes(f["mchi"])["unit"] = "GeV"
    dset = create_dataset(f, "sigma", dtype, ())
    attributes(dset)["unit"] = "cm^2"

    for (i, mchi) in enumerate(mchis)
        grp = create_group(f, "m_$(i-1)")
        dset = create_dataset(grp, "b", dtype, (shape[3],))
        attributes(dset)["unit"] = "rad"
        dset = create_dataset(grp, "l", dtype, (shape[2],))
        attributes(dset)["unit"] = "rad"
        dset = create_dataset(grp, "T", dtype, (shape[1],))
        attributes(dset)["unit"] = "GeV"
        dset = create_dataset(grp, "intensity", dtype, shape)
        attributes(dset)["unit"] = "GeV^-1 cm^-2 s^-1 sr^-1"
        attributes(dset)["dimension"] = "(b, l, T)"
        grp["mchi"] = mchi / units.GeV
        attributes(grp["mchi"])["unit"] = "GeV"
    end
    close(f)
end

"""Dump 3-D intensity to a HDF5 file. """
function h5dumpflux(
    bdm::BDM, filename::String, b::AbstractVector, l::AbstractVector, Tchi::AbstractVector;
    groupname::Union{Nothing, String}=nothing, kwargs...
)
    f = h5open(filename, "cw")
    _require_dataset(f, "sigma", xsec0(bdm) / units.cm2)

    grp = isnothing(groupname) ? f : f[groupname]
    _require_dataset(grp, "mchi", dmmass(bdm) / units.GeV)
    _require_dataset(grp, "b", b)
    _require_dataset(grp, "l", l)
    _require_dataset(grp, "T", Tchi / units.GeV)

    # With reverse order in Julia since HDF5 library is C ordered
    intensity = zeros(length(Tchi), length(l), length(b))
    Threads.@threads for ib in eachindex(b)
        for il in eachindex(l), iT in eachindex(Tchi)
            intensity[iT, il, ib] = dmflux(bdm, Tchi[iT], b[ib], l[il]; rtol=0.01, kwargs...)
        end
    end
    fluxunit = units.GeV^-1 * units.sec^-1 * units.cm^-2
    _require_dataset(grp, "intensity", intensity / fluxunit)
    close(f)
end

"""Dump 2D (b, l) distribution to a hdf5 file. """
function h5dump_bl_distribution(
    bdm::BDM, filename::String, b::AbstractVector, l::AbstractVector, Tchimin, Tchimax;
    groupname::Union{Nothing, String}=nothing, kwargs...
)
    f = h5open(filename, "cw")
    grp = isnothing(groupname) ? f : create_group(f, groupname)

    bldist = zeros(length(l), length(b))
    fluxunit = units.sec^-1 * units.cm^-2

    Threads.@threads for ib in eachindex(b)
        for il in eachindex(l)
            res::typeof(fluxunit), _ = quad(
                log(Tchimin), log(Tchimax); atol=0, rtol=0.01, kwargs...
            ) do lnT
                T = exp(lnT)
                return dmflux(bdm, T, b[ib], l[il]; rtol=0.01, kwargs...) * T
            end
            bldist[il, ib] = res
        end
    end

    grp["sigma"] = xsec0(bdm) / units.cm2
    attributes(grp["sigma"])["unit"] = "cm^2"
    grp["mchi"] = dmmass(bdm) / units.GeV
    attributes(grp["mchi"])["unit"] = "GeV"
    grp["b"] = b
    attributes(grp["b"])["unit"] = "rad"
    grp["l"] = l
    attributes(grp["l"])["unit"] = "rad"
    grp["Tmin"] = Tchimin / units.GeV
    attributes(grp["Tmin"])["unit"] = "GeV"
    grp["Tmax"] = Tchimax / units.GeV
    attributes(grp["Tmax"])["unit"] = "GeV"
    grp["distribution"] = bldist / fluxunit
    attributes(grp["distribution"])["unit"] = "cm^-2 s^-1 sr^-1"
    attributes(grp["distribution"])["dimension"] = "(b, l)"
    close(f)
end

"""Load and interpolate the 3-D intensity dumped by `h5loadflux`. """
function h5loadflux(filename, groupname=nothing, log_scale=false)
    f = h5open(filename, "r")
    grp = isnothing(groupname) ? f : f[groupname]
    b = read(grp, "b")
    l = read(grp, "l")
    T = read(grp, "T")
    flux = read(grp, "intensity")
    lnflux = @. ifelse(flux <= 0, log(1e-300), log(flux))
    logspectrum = linear_interpolation(
        (b, l, log.(T / units.GeV)), permutedims(lnflux, (3, 2, 1)), extrapolation_bc=-Inf
    )

    local wrapper::Function
    if log_scale
        # log probability for MCMC
        wrapper = (b, l, Tchi) -> logspectrum((b, l, log(Tchi / units.GeV)))
    else
        wrapper = (b, l, Tchi) -> let fluxunit = units.GeV^-1 * units.sec^-1 * units.cm^-2
            exp(logspectrum(b, l, log(Tchi / units.GeV))) * fluxunit
        end
    end

    blobs = (b, l, T, read(f, "sigma") * units.cm2, read(f, "mchi") * units.GeV)
    close(f)
    return wrapper, blobs
end

"""Load and interpolate the 2-D distribution dumped by `h5dump_bl_distribution`. """
function h5load_bl_distribution(filename, groupname=nothing, log_scale=false)
    f = h5open(filename, "r")
    grp = isnothing(groupname) ? f : f[groupname]
    b = read(grp, "b")
    l = read(grp, "l")

    flux = read(grp, "distribution")
    lnflux = @. ifelse(flux <= 0, -Inf, log(flux))
    logspectrum = linear_interpolation(
        (b, l), permutedims(lnflux, (2, 1)), extrapolation_bc=-Inf
    )

    local wrapper::Function
    if log_scale
        # log probability for MCMC
        wrapper = (b, l) -> logspectrum(b, l)
    else
        wrapper = (b, l) -> let fluxunit = units.sec^-1 * units.cm^-2
            exp(logspectrum(b, l)) * fluxunit
        end
    end

    blobs = (b, l, read(grp, "sigma") * units.cm2, read(grp, "mchi") * units.GeV)
    return wrapper, blobs
end

"""Line of sight boundary. """
function _los_bound(b, l, rmax, zmax, rsun)
    rsun_cosl = rsun * cos(l)
    d = rsun_cosl + sqrt(rsun_cosl^2 + rmax^2 - rsun^2)

    # boundary is on the z = +/- zmax planes
    if abs(tan(b)) * d > zmax
        return zmax / abs(sin(b))
    end

    # TODO check 0
    return d / cos(b)
end


"""Dark matter profile times cosmic ray flux. """
function _los_integrand(bdm::CRDM, r, b, l, Ek, particle)
    return (bdm.crdist[particle](r, b, l, Ek)
            * dmdensity_galactic(bdm.dmprofile, r, b, l, rsun=bdm.rsun))
end


function _parse_particles(bdm::CRDM, particles::String)
    if lowercase(particles) == "all"
        return [key for key in keys(bdm.crdist)]
    end
    return [particles]
end
_parse_particles(::CRDM, particles::Vector{String}) = particles

function _require_dataset(group, name, data)
    if haskey(group, name)
        dset = group[name]
    else
        dset = create_dataset(group, name, eltype(data), size(data))
    end
    write(dset, data)
    return nothing
end
