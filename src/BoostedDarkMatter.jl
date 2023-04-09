module BoostedDarkMatter

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
export dmmass, dmmass!, xsec0, xsec0!, dxsecdT4, recoil_spectrum, set_parameters!
export mediatormass, mediatormass!

# CosmicRayBoostedDarkMatter.jl
export CRDM, kfactor, kfactor!, crdist!
export CRDistribution, CRFlux, EnergyType, crflux, crflux_Ekn
export CRFluxSBPLElectron, CRFluxLISElectron


abstract type BDM end

Broadcast.broadcastable(bdm::BDM) = Ref(bdm)

dmflux(::BDM, Tchi) = error("unimplemented")
dmmass(::BDM) = error("unimplemented")

include("Kinematics.jl")
using .Kinematics
include("Constants.jl")
include("Dynamics.jl")
include("CosmicRayBoostedDarkMatter/CosmicRayBoostedDarkMatter.jl")


end
