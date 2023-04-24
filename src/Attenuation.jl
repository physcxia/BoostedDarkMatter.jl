using Roots: find_zero, ITP
using FiniteDifferences: central_fdm
# using OrdinaryDiffEq: ODEProblem, AutoTsit5, Rosenbrock23, solve
using OrdinaryDiffEq: solve, ODEProblem, Rodas4, AutoVern7
# using OrdinaryDiffEq: Rodas5P, AutoVern9, RadauIIA3

Base.@kwdef mutable struct Attenuation{T <: Number, U <: XSec}
    target_density::Dict = Dict(
        nuc => r * 2.7 * units.g_cm3 / BoostedDarkMatter.NUCLEUS_MASS[nuc] for (nuc, r) in
        Dict(
            "O"  => 46.6e-2 / 0.985,
            "Si" => 27.7e-2 / 0.985,
            "Al" => 8.1e-2  / 0.985,
            "Fe" => 5.0e-2  / 0.985,
            "Ca" => 3.6e-2  / 0.985,
            "K"  => 2.8e-2  / 0.985,
            "Na" => 2.6e-2  / 0.985,
            "Mg" => 2.1e-2  / 0.985
        )
    )
    xsec::U = XSecDMNucleusConstant()
    Tcut::T = 1e-2 * units.eV
    _dTdz = nothing
end

Broadcast.broadcastable(atten::Attenuation) = Ref(atten)

""" Cache dTdz """
function cache_dTdz!(
    atten::Attenuation,
    Tchi::AbstractVector=geomspace(1e-2*units.eV, 1e10*units.GeV, 1000);
    kwargs...
)
    dtdz = -dTdz.(atten::Attenuation, Tchi; usecache=false, kwargs...)
    atten._dTdz = loginterpolator(Tchi, dtdz)
    return nothing
end

function dTdz(atten::Attenuation, T; usecache=true, atol=0, rtol=1e-6, kwargs...)
    T > atten.Tcut || return 0
    !usecache || isnothing(atten._dTdz) || return -atten._dTdz(T)

    Eunit = oneunit(T)

    old_target = gettarget(atten.xsec)
    res = 0
    for (target, density) in atten.target_density
        settarget!(atten.xsec, target)
        m_in = dmmass(atten.xsec)
        m_t = BoostedDarkMatter.NUCLEUS_MASS[target]
        Trmax = T4_max(atten.xsec, T, m_in, m_t)

        # res_quad, _ = quad(zero(Trmax), Trmax; atol=atol, rtol=rtol, kwargs...) do Tr
        #     dxsecdT4(atten.xsec, Tr, T, m_in, m_t) * Tr * density
        # end

        Trmax > atten.Tcut || continue
        res_quad, _ = quad(
            log(atten.Tcut / Eunit),
            log(Trmax / Eunit);
            atol=atol, rtol=rtol, kwargs...
        ) do lnTr
            Tr = exp(lnTr)
            dxsecdT4(atten.xsec, Tr, T, m_in, m_t) * Tr^2 * density
        end

        res += res_quad
        if res_quad <= 0
            @error "res_quad = $res_quad <= 0, T = $T, Trmax = $Trmax"
        end
    end
    settarget!(atten.xsec, old_target)
    return -res
end

function dlnTdz(atten::Attenuation, lnT; atol=0, rtol=1e-5)
    T = exp(lnT)
    return dTdz(atten, T; atol=atol, rtol=rtol) / T
end


function lnTzlnT0(atten::Attenuation, lnT0, z; options...)
    f(lntz, _, _) = (lnT0 >= lntz > log(atten.Tcut)) ? dlnTdz(atten, lntz) : zero(lnT0/z)
    prob = ODEProblem(f, lnT0, (zero(z), z))
    # sol = solve(prob, AutoDP5(Vern9()); save_everystep=false, options...)
    # sol = solve(prob, AutoVern9(Rodas4P()); save_everystep=false, options...)
    # sol = solve(prob, RadauIIA3(); save_everystep=false, abstol=0, reltol=1e-2, options...)
    # alg = Rodas5P()
    alg = AutoVern7(Rodas4())
    sol = solve(prob, alg;
        dt=1*units.mm,
        dtmin=1*units.mm,
        force_dtmin=true,
        save_everystep=false,
        abstol=0,
        reltol=1e-2,
        options...
    )
    return max(log(atten.Tcut), sol.u[end])
end
# function lnTzlnT0(atten::Attenuation, lnT0, z; Tzcut=1e-2units.eV, options...)
#     return log(TzT0(atten, exp(lnT0), z; Tzcut=Tzcut, options...))
# end


"""
    TzT0(atten::Attenuation, T0, z; options...)

Parameters
----------

- `atten`: Attenuation object.
- `T0`: Kinetic energy at the surface (`z` = 0).
- `z`: The depth.

"""
function TzT0(atten::Attenuation, T0, z; options...)
    lnT0 = log(T0)
    return exp(lnTzlnT0(atten, lnT0, z; options...))
end
# function TzT0(atten::Attenuation, T0, z; Tzcut=1e-2units.eV, options...)
#     f(tz, _, _) = (T0 >= tz > Tzcut) ? dTdz(atten, tz) : zero(T0/z)
#     prob = ODEProblem(f, T0, (zero(z), z))
#     sol = solve(prob, AutoDP5(Vern9()); save_everystep=false, abstol=0, reltol=1e-6, options...)
#     # sol = solve(prob, RadauIIA3(); save_everystep=false, abstol=0, reltol=1e-2, options...)
#     # sol = solve(prob, AutoVern9(Rodas4P()); save_everystep=false, abstol=0, reltol=1e-6, options...)
#     return max(Tzcut, sol.u[end])
# end

"""
    T0Tz(atten::Attenuation, Tz, z; T0_max=10, options...)

Inverse of [`TzT0`](@ref).

Parameters
----------

- `atten`: Attenuation object.
- `Tz`: Kinetic energy at depth `z`.
- `z`: The depth.

"""
# function T0Tz(atten::Attenuation, Tz, z; T0_max=10, options...)
#     Tz > zero(Tz) || return zero(Tz)
#     lnTz = log(Tz)
#     lnT0max = log(T0_max)
#     f(lnT0) = lnTzlnT0(atten, lnT0, z; options...) - lnTz

#     if f(lnTz) * f(lnT0max) > zero(lnTz)
#         # @warn "no root found"
#         return T0_max
#     end
#     return exp(find_zero(f, (lnTz, lnT0max), ITP()))
# end
# function T0Tz(atten::Attenuation, Tz, z; T0_max=10, options...)
#     Tz > zero(Tz) || return zero(Tz)
#     f(T0) = TzT0(atten, T0, z; options...) - Tz

#     if f(Tz) * f(T0_max) > zero(Tz)
#         # @warn "no root found"
#         return T0_max
#     end
#     return find_zero(f, (Tz, T0_max), ITP())
# end
function T0Tz(atten::Attenuation, Tz, z; T0_max=10, options...)
    return exp(lnT0lnTz(atten, log(Tz), z; T0_max, options...))
end


# function lnT0lnTz(atten::Attenuation, lnTz, z; T0_max, options...)
#     return log(T0Tz(atten, exp(lnTz), z; T0_max=T0_max, options...))
# end


function lnT0lnTz(atten::Attenuation, lnTz, z; T0_max=10, options...)
    lnT0max = log(T0_max)
    f(lnT0) = lnTzlnT0(atten, lnT0, z; options...) - lnTz

    # no energy loss
    if f(lnTz) * f(lnT0max) > zero(lnTz)
        @warn "no root found"
        return lnT0max
    end
    return find_zero(f, (lnTz, lnT0max), ITP())
end


function dT0dTz(atten::Attenuation, Tz, T0, z; T0_max, options...)
    T0 < T0_max || return zero(T0/Tz)
    deriv = central_fdm(5, 1)(
        lntz -> lnT0lnTz(atten, lntz, z; T0_max=T0_max, options...), log(Tz)
    ) * T0 / Tz
    return max(zero(deriv), deriv)
end


function dmflux(atten::Attenuation, Tz, z, flux0; T0_max=10, options...)
    T0 = T0Tz(atten, Tz, z; T0_max=T0_max, options...)
    (zero(T0) < T0 < T0_max) || return zero(flux0(T0_max))
    return flux0(T0) * dT0dTz(atten, Tz, T0, z; T0_max=T0_max, options...)
end


function dmflux(
    atten::Attenuation,
    Tz::AbstractArray,
    z,
    T0::AbstractArray,
    flux0::AbstractArray;
    options...
)
    flux_0 = loginterpolator(T0, flux0)
    return dmflux.(atten, Tz, z, Ref(flux_0); T0_max=maximum(T0), options...)
end


function dTzdT0(atten::Attenuation, T0, z)
    Tz = TzT0(atten, T0, z)
    return dTzdT0(atten, T0, Tz, z)
end
function dTzdT0(atten::Attenuation, T0, Tz, z)
    T0 > zero(T0) || return zero(T0/T0)
    deriv = central_fdm(5, 1)(
        lnt0 -> lnTzlnT0(atten, lnt0, z), log(T0)
    ) * Tz / T0
    return max(zero(deriv), deriv)
end


function dmspectrum_z(atten::Attenuation, T0, flux0, z)
    Tz = TzT0(atten, T0, z)
    fluxz = flux0 / dTzdT0(atten, T0, z)
    return Tz, fluxz
end

function dmspectrum_z(atten::Attenuation, T0::AbstractArray, flux0::AbstractArray, z)
    Tz = TzT0.(atten, T0, z)
    fluxz = flux0 ./ dTzdT0.(atten, T0, z)
    return Tz, fluxz
end
