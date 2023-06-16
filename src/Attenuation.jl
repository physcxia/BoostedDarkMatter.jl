using Roots: find_zero, ITP
using FiniteDifferences: central_fdm, extrapolate_fdm
using OrdinaryDiffEq: solve, ODEProblem, Rodas4, AutoVern7

Base.@kwdef mutable struct Attenuation{D <: Dict, T <: Number, U <: XSec}
    target_density::D = Dict(
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
    Tmin::T = 1e-2 * units.keV
    Tmax::T = 1e2 * units.GeV
    Tcut::T = 4 * units.eV
    _dTdz = nothing
end

Broadcast.broadcastable(atten::Attenuation) = Ref(atten)

xsec0!(atten::Attenuation, sigma0) = xsec0!(atten.xsec, sigma0)
xsec0(atten::Attenuation) = xsec0(atten.xsec)
dmmass!(atten::Attenuation, mchi) = dmmass!(atten.xsec, mchi)
dmmass(atten::Attenuation) = dmmass(atten.xsec)

""" Cache dTdz """
function cache_dTdz!(
    atten::Attenuation,
    Tchi::AbstractVector=geomspace(atten.Tmin, atten.Tmax, 10000);
    kwargs...
)
    atten.Tmin = minimum(Tchi)
    atten.Tmax = maximum(Tchi)
    dtdz = -dTdz.(atten::Attenuation, Tchi; usecache=false, kwargs...)
    atten._dTdz = loginterpolator(Tchi, dtdz)
    return nothing
end
function cache_dTdz!(atten::Attenuation, Tmin, Tmax, n; kwargs...)
    cache_dTdz!(atten, geomspace(Tmin, Tmax, n); kwargs...)
end

@doc raw"""
    dTdz(atten::Attenuation, T; usecache=true, atol=0, rtol=1e-6, kwargs...)

Calculate ``\frac{dT}{dz}`` at `T` with the attenuation model `atten`,

```math
\frac{dT}{dz} = -\int_0^{T_\text{r}^{\max} n \frac{dσ}{dT_\text{r}} T_\text{r} dT_\text{r}.
```

# Arguments

- `atten::Attenuation`: The Attenuation object.
- `T`: Kinetic energy.

# Keywords

- `usecache`: If `true` and the cache exists, then the cached `_dTdz` will be used.
- `atol=0`: Absolute tolerence.
- `rtol=1e-6`: Relative tolerence.

"""
function dTdz(atten::Attenuation, T; usecache=true, atol=0, rtol=1e-6, kwargs...)
    zunit = oneunit(1 / (oneunit(valtype(atten.target_density)) * xsec0(atten.xsec)))
    dTdz_unit = oneunit(T / zunit)
    T > atten.Tmin || return zero(dTdz_unit)
    !usecache || isnothing(atten._dTdz) || return -atten._dTdz(T)::typeof(dTdz_unit)

    Eunit = oneunit(T)

    old_target = gettarget(atten.xsec)
    res = zero(dTdz_unit)
    for (target, density) in atten.target_density
        settarget!(atten.xsec, target)
        m_in = dmmass(atten.xsec)
        m_t = NUCLEUS_MASS[target]
        Trmax = T4_max(atten.xsec, T, m_in, m_t)
        Trmax > atten.Tcut || continue

        res_quad::typeof(dTdz_unit) = zero(dTdz_unit)
        if iszero(atten.Tmin)
            res_quad, _ = quad(zero(Trmax), Trmax; atol=atol, rtol=rtol, kwargs...) do Tr
                dxsecdT4(atten.xsec, Tr, T, m_in, m_t) * Tr * density
            end
        else
            res_quad, _ = quad(
                log(atten.Tcut / Eunit),
                log(Trmax / Eunit);
                atol=atol, rtol=rtol, kwargs...
            ) do lnTr
                Tr = exp(lnTr)
                dxsecdT4(atten.xsec, Tr, T, m_in, m_t) * Tr^2 * density
            end
        end

        res += res_quad
        if res_quad <= 0
            @error "res_quad = $res_quad <= 0, T = $T, Trmax = $Trmax"
        end
    end
    settarget!(atten.xsec, old_target)
    return -res
end


function T_z(atten::Attenuation, T0, z; kwargs...)
    f(tz, _, _) = (T0 >= tz > atten.Tmin) ? dTdz(atten, tz) : zero(T0/z)
    prob = ODEProblem(f, T0, (zero(z), z))
    alg = AutoVern7(Rodas4())
    sol = solve(prob, alg;
        dt=1e-10*units.m,
        dtmin=1e-10*units.m,
        force_dtmin=true,
        abstol=0,
        reltol=1e-8,
        kwargs...
    )
    return sol.t, sol.u
end


"""
    TzT0(atten::Attenuation, T0, z; kwargs...)

Parameters
----------

- `atten`: Attenuation object.
- `T0`: Kinetic energy at the surface (`z` = 0).
- `z`: The depth.

"""
function TzT0(atten::Attenuation, T0, z; kwargs...)
    f(tz, _, _) = (T0 >= tz > atten.Tcut) ? dTdz(atten, tz) : zero(T0/z)
    prob = ODEProblem(f, T0, (zero(z), z))
    alg = AutoVern7(Rodas4())
    sol = solve(
        prob, alg;
        dt=1e-12*units.m,
        dtmin=1e-12*units.m,
        force_dtmin=true,
        save_everystep=false,
        abstol=0,
        reltol=1e-6,
        kwargs...
    )
    return max(atten.Tcut, sol.u[end])::typeof(T0)
end


"""
    T0Tz(atten::Attenuation, Tz, z; kwargs...)

Inverse of [`TzT0`](@ref).

Parameters
----------

- `atten`: Attenuation object.
- `Tz`: Kinetic energy at depth `z`.
- `z`: The depth.

"""
function T0Tz(atten::Attenuation, Tz, z; kwargs...)
    Tz > zero(Tz) || return zero(Tz)
    f(T0) = TzT0(atten, T0, z; kwargs...) - Tz

    f(Tz) < zero(Tz) || return Tz  # no energy loss
    f(atten.Tmax) > zero(Tz) || return typemax(Tz)  # not in the range
    return find_zero(f, (Tz, atten.Tmax), ITP())
end


function dT0dTz(atten::Attenuation, Tz, T0, z; kwargs...)
    atten.Tmin < Tz <= T0 < atten.Tmax || return zero(T0/Tz)
    function f(tz)
        T0Tz(atten, tz, z; kwargs...)
    end
    max_range = min(Tz - atten.Tmin, atten.Tmax - Tz)
    deriv = central_fdm(5, 1, factor=1, max_range=max_range)(f, Tz)
    return max(zero(deriv), deriv)::typeof(T0/Tz)
end


function dmflux(atten::Attenuation, Tz, z, flux0; kwargs...)
    T0 = T0Tz(atten, Tz, z; kwargs...)
    (zero(T0) < Tz <= T0 < atten.Tmax) || return zero(flux0(atten.Tmax))
    return flux0(T0) * dT0dTz(atten, Tz, T0, z; kwargs...)
end


function dmflux(
    atten::Attenuation,
    Tz::AbstractArray,
    z,
    T0::AbstractArray,
    flux0::AbstractArray;
    kwargs...
)
    flux_0 = loginterpolator(T0, flux0)
    return dmflux.(atten, Tz, z, Ref(flux_0); kwargs...)
end


function dTzdT0(atten::Attenuation, T0, z)
    Tz = TzT0(atten, T0, z)
    return dTzdT0(atten, T0, Tz, z)
end
function dTzdT0(atten::Attenuation, T0, Tz, z; kwargs...)
    T0 > zero(T0) || return zero(T0/T0)
    deriv = central_fdm(5, 1)(t0 -> TzT0(atten, t0, z; kwargs), T0)
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


function mean_free_path(atten::Attenuation, Tchi; kwargs...)
    old_target = gettarget(atten.xsec)
    λ = 0
    for (target, density) in atten.target_density
        settarget!(atten.xsec, target)
        σtot, _ = total_xsec(
            atten.xsec, Tchi, dmmass(atten.xsec), NUCLEUS_MASS[target];
            Trcut=atten.Tcut, kwargs...)
        λ += density * σtot
    end
    settarget!(atten.xsec, old_target)
    return 1 / λ
end
