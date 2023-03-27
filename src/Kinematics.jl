"""
The Kinematics sub-module contains kinematics function for elastic two body scattering,

```math
p_1 + p_2 → p_3 + p_4.
```

The convention is ``1 → 3`` and ``2 → 4``. Currently, the functions in this module are given
for the laboratory frame where the particle 2 is at rest.

"""
module Kinematics

export
    reduce_m,
    T4_max,
    T4_max_nonrel,
    T1_min,
    T1_min_nonrel,
    p1_min

using NumericalTools: sqrtm1


@doc raw"""
    reduce_m(m1, m2)

Reduced mass ``μ_{12} = \frac{m_1 m_2}{m_1 + m_2}``.

"""
function reduce_m(m1, m2)
    return m1 * m2 / (m1 + m2)
end


@doc raw"""
    T4_max(T1, m1, m2)

Maximum recoil energy of the stationary target particle 2 in elastic scattering,

```math
T_4^{\max} = \frac{2m_2 T_1 (T_1 + 2m_1)}{2m_2 T_1 + (m_1 + m_2)^2}.
```

# Arguments

- `T1`: Kinetic energy of the incident particle 1.
- `m1`: Mass of the incident particle 1.
- `m2`: Mass of the target particle 2.

"""
function T4_max(T1, m1, m2)
    return (2 * m2 * T1 * (T1 + 2 * m1)) / (2 * m2 * T1 + (m1 + m2)^2)
end


@doc raw"""
    T4_max_nonrel(T1, m1, m2)

Non-relativistic version of [`T4_max`](@ref),

```math
T_4^{\max} = \frac{4 μ_{12} T_1}{m_1 m_2} = \frac{4 m_1 m_2 T_1}{(m_1 + m_2)^2}.
```

"""
function T4_max_nonrel(T1, m1, m2)
    return 4 * reduce_m(m1, m2)^2 / (m1 * m2) * T1
end


@doc raw"""
    T1_min(T4, m1, m2)

Minimal kinetic energy of particle 1 required to obtain a recoil energy of `T4` for the
stationary target 2 in elastic scattering,

```math
\begin{aligned}
    T_1^{\min}
&=
    \left(\frac{T_4}{2} - m_1\right)
    \left(1 ± \sqrt{1 + \frac{2T_4}{m_2}\frac{(m_1+m_2)^2}{(2m_1-T_4)^2}}\right)\\
&=
    \frac{T_4}{2} - m_1
    +
    \sqrt{\left(\frac{T_4}{2} - m_1\right)^2 + \frac{T_4(m_1+m_2)^2}{2m_2}}.
\end{aligned}
```

In the case of
``T_4 < 2m_1`` and ``\left(\frac{T_4}{2} - m_1\right)^2 ≫ \frac{T_4(m_1+m_2)^2}{2m_2}``,
we need to expand the square root to avoid precision loss, and the result is
```math
T_4^{\min} ≈ \frac{T_4(m_1+m_2)^2}{2m_2(2m_1 - T_4)}.
```

# Arguments

- `T4`: Recoil kinetic energy of the target particle 2.
- `m1`: Mass of the incident particle 1.
- `m2`: Mass of the target particle 2.

"""
function T1_min(T4, m1, m2)
    fac1 = T4 / 2 - m1
    fac2 = 2T4 / m2 * ((m1 + m2) / (2m1 - T4))^2
    if fac1 < zero(T4)
        return -fac1 * sqrtm1(fac2)
    end
    return fac1 * (1 + sqrt(1 + fac2))
end


@doc raw"""
    T1_min_nonrel(T4, m1, m2)

Non-relativistic version of [`T1_min`](@ref),

```math
    T_1^{\min} = \frac{m_1 m_2 T_4}{4 μ_{12}^2} = \frac{(m_1 + m_2)^2 T_4}{4m_1 m_2}.
```

"""
function T1_min_nonrel(T4, m1, m2)
    return m1 * m2 * T4 / (4 * reduce_m(m1, m2)^2)
end


function p1_min(q, ΔE, m1, m2)
    if ΔE >= q
        @warn "Kinematics Invalid: ΔE = $ΔE >= q = $q, Inf returned"
        return Inf
    end
    return q / 2 + ΔE / 2 * sqrt(1 + 4m1^2 / ((q + ΔE) * (q - ΔE)))
end


function T1_min(q, ΔE, m1, m2)
    return m1 * sqrtm1((p1_min(q, ΔE, m1, m2) / m1)^2)
end


end
