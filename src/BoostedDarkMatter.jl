module BoostedDarkMatter

export units
export dmflux

# Particles.jl
export Particle, particle_A, particle_Z, particle_mass, ZA_to_pdgcode, pdgcode

# Kinematics.jl
export
    reduce_m,
    T4_max,
    T4_max_nonrel,
    T1_min,
    T1_min_nonrel,
    p1_min_q,
    T1_min_q

# Attenuation.jl
export Attenuation, dmspectrum_z, cache_dTdz!, mean_free_path
export TzT0, T0Tz, dTzdT0, dT0dTz, T_z

# Dynamics.jl
export XSec, XSecElastic, XSecDMElectronElastic, XSecDMElectronScalarMediator
export XSecDMElectronVectorMediator, XSecDMElectronBound
export XSecDMNucleusConstant
export dmmass, dmmass!, xsec0, xsec0!, dxsecdT4, recoil_spectrum, set_parameters!
export mediatormass, mediatormass!, total_xsec

# CosmicRayBoostedDarkMatter.jl
export BDM, CRDM, kfactor, kfactor!, crdist!, selectcr!
export CRDistribution, CRFlux, EnergyType, crflux, crflux_Ekn, crflux_Ek
export CRFluxSBPLElectron, CRFluxLISElectron
export CRFluxLISHelMod2017
export CRDGalprop, CRDGalpropCylindrical
export make_crflux_dict_galactic
export h5dumpflux, h5loadflux, h5dump_bl_distribution, h5load_bl_distribution



abstract type BDM end

Broadcast.broadcastable(bdm::BDM) = Ref(bdm)

dmflux(::BDM, Tchi::Number) = error("unimplemented")
function dmflux(bdm::BDM, Tchi::AbstractVector; kwargs...)
    flux = zero(Tchi)  # TODO support unitful type
    Threads.@threads for i in eachindex(Tchi)
        flux[i] = dmflux(bdm, Tchi[i]; kwargs...)
    end
    return flux
end

dmmass(::BDM) = error("unimplemented")

using NumericalTools: geomspace, loginterpolator, sqrtm1

include("Units.jl")
import .Units
const units = Units
include("Kinematics.jl")
using .Kinematics
include("Constants.jl")
include("Particles.jl")
include("Dynamics.jl")
include("CosmicRayBoostedDarkMatter/CosmicRayBoostedDarkMatter.jl")
include("Attenuation.jl")


end
