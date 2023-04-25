function particle_Z(pdgc::Int)
    d = digits(pdgc)
    Z = sum(d[i+4] * 10^(i-1) for i in 1:3)
    return Z
end
particle_Z(name::String) = particle_Z(pdgcode(name))
function particle_A(pdgc::Int)
    d = digits(pdgc)
    A = sum(d[i+1] * 10^(i-1) for i in 1:3)
    return A
end
particle_A(name::String) = particle_A(pdgcode(name))
particle_ZA(pdgc::Int) = (particle_Z(pdgc), particle_A(pdgc))
particle_ZA(name::String) = particle_ZA(pdgcode(name))
ZA_to_pdgcode(Z::Int, A::Int) = 1000000000 + Z * 10000 + A * 10
function pdgcode(name::String)
    if lowercase(name) == "proton"
        return 1000010010
    end
    return ZA_to_pdgcode(NUCLEUS_Z[name], NUCLEUS_A[name])
end

struct Particle{T}
    Z::Int
    A::Int
    mass::T
end
Particle(Z::Int, A::Int) = Particle(Z, A, A == 1 ? PROTON_MASS : A * ATOMIC_MASS)

particle_mass(name::String) = NUCLEUS_MASS[name]
