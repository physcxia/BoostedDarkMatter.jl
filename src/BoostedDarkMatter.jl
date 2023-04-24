module BoostedDarkMatter

export units
export dmflux

# Kinematics.jl
export
    reduce_m,
    T4_max,
    T4_max_nonrel,
    T1_min,
    T1_min_nonrel,
    p1_min_q,
    T1_min_q

# Dynamics.jl
export XSec, XSecElastic, XSecDMElectronElastic, XSecDMElectronScalarMediator
export XSecDMElectronVectorMediator, XSecDMElectronBound
export XSecDMNucleusConstant
export settarget!, gettarget
export dmmass, dmmass!, xsec0, xsec0!, dxsecdT4, recoil_spectrum, set_parameters!
export mediatormass, mediatormass!, total_xsec

# CosmicRayBoostedDarkMatter.jl
export CRDM, kfactor, kfactor!, crdist!, select_cr!
export CRDistribution, CRFlux, EnergyType, crflux, crflux_Ekn
export CRFluxSBPLElectron, CRFluxLISElectron
export CRFluxLISHelMod2017

# Attenuation.jl
export Attenuation, dmspectrum_z, cache_dTdz!


abstract type BDM end

Broadcast.broadcastable(bdm::BDM) = Ref(bdm)

dmflux(::BDM, Tchi) = error("unimplemented")
dmmass(::BDM) = error("unimplemented")

using NumericalTools: geomspace, loginterpolator, sqrtm1

include("Units.jl")
import .Units
const units = Units
include("Kinematics.jl")
using .Kinematics
include("Constants.jl")
include("Dynamics.jl")
include("CosmicRayBoostedDarkMatter/CosmicRayBoostedDarkMatter.jl")
include("Attenuation.jl")


end
