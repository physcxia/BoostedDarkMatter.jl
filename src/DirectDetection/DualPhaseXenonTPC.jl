using SpecialFunctions: loggamma
using Distributions: Poisson, cdf

abstract type DirectDetectionExperiment end

abstract type DualPhaseXenonTPC <: DirectDetectionExperiment end

xsec_limits(::DirectDetectionExperiment, dRdE) = error("unimplemented")

struct Xenon1T <: DualPhaseXenonTPC
    efficiency::Function
    events_bg
    events_dt
    Er_bins
    exposure
end

struct PandaX4T <: DualPhaseXenonTPC
end

function _find_sigma(
    xenon::DualPhaseXenonTPC, events_sig, sigmas, method, alpha, alpha_bin,
    chi2_min, delta_chi2
)
    event_bg = sum(xenon.events_bg)
    event_dt = sum(xenon.events_dt)
    for sig in sigmas
        if method == "binned" && any(
            @. cdf(Poisson(xenon.events_bg + events_sig), xenon.events_dt) < alpha_bin)
            return sig
        end
        if method == "onebin" && cdf(Poisson(event_bg + sum(events_sig)), event_dt) < alpha
            return sig
        end
        if method == "binned_chi2"
            chi2 = -2 * sum(log_poisson.(events_sig + xenon.events_bg, xenon.events_dt))
            if chi2 - chi2_min > delta_chi2
                return sig
            end
        end
        if method == "number" && sum(events_sig) > event_dt
            return sig
        end
    end

    return nothing
end

function xsec_limits(
    xenon::Xenon1T, dRdE, xsecs::AbstractVector;
    confidence_level::Real=0.9, delta_chi2::Real=2.71, method::String="binned"
)
    xsecs = sort(xsecs)
    alpha = 1 - confidence_level
    alpha_bin = 1 - confidence_level^(1/length(xenon.events_dt))
    chi2_min = -2 * sum(log_poisson.(xenon.events_bg, xenon.events_dt))
    events_sig = dRdE.(xenon.Er_bins)  # FIXME
    lower = _find_sigma(xenon, events_sig, xsecs, method, alpha, alpha_bin, chi2_min, delta_chi2)
    if isnothing(lower)
        return (NaN, NaN)
    end
    upper = _find_sigma(xenon, events_sig, reverse(xsecs), method, alpha, alpha_bin, chi2_min, delta_chi2)
    return (lower, upper)
end


"""Log Poisson distribution.

Parameters
----------
n : array
    Number of events.
mean : array
    Mean of Poisson.

Returns
-------
array
    Log poisson PDF.

"""
function log_poisson(n, mean)
    if n > 150
        return n * log(mean) - mean - n * (log(n) - 1)
    end

    return n * log(mean) - mean - loggamma(n + 1)
end
