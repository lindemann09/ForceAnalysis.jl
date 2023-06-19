struct OnsetCriterion
    min_increase::Float64  # min change between two samples
    min_increase_in_window::Float64
    window_size::Int
    function OnsetCriterion(min_increase, min_increase_in_window, window_size)
        (min_increase >= 0 && min_increase_in_window >= 0
         && window_size >= 0) || throw(ArgumentError(
            "No negative parameter allowed for OnsetCriterion"))
        return new(min_increase, min_increase_in_window, window_size)
    end
end;

struct ResponseBounds
    onset::Int
    offset::Int
    function ResponseBounds(a, b)
        a > b && b > 0 && throw(ArgumentError("Onset must be before offset."))
        return new(a, b)
    end
end

function OnsetCriterion(;
    min_increase::Real,
    min_increase_in_window::Real,
    window_size::Integer
)
    return OnsetCriterion(min_increase, min_increase_in_window, window_size)
end

function detect_responses(
    fp::ForceProfiles,
    criterion::OnsetCriterion
)::Vector{ResponseBounds}
    # returns onset, offset and peak force
    #   takes into account zero sample. That is, onset and off are samples after
    #       zero sample
    rtn = ResponseBounds[]
    for f in eachrow(fp.dat)
        onset = response_onset(f, criterion)
        if onset > 0
            offset = response_offset(f, criterion; onset)
            push!(rtn, ResponseBounds(onset, offset))
        else
            push!(rtn, ResponseBounds(-1, -1))
        end
    end
    return rtn
end

function response_onset(
    force_vector::AbstractVector{<:FloatOrMissing},
    criterion::OnsetCriterion
)::Int
    # returns sample in vector with onset

    # function for onset detection, finds first onset
    #
    # detects if  difference threshold (current-previous) is exceeded and force has
    # INCREASED after n trials for at least the increase_criterion
    #
    # returns sample of onset or -1

    force_diff = diff(force_vector)
    l = length(force_vector)
    l - 1 == length(force_diff) || throw(ArgumentError(
        "length(force_difference) must be length(force_vector)-1 ."))

    for (c, d) in enumerate(force_diff)
        if !ismissing(d) && d > criterion.min_increase
            x = c + 1 # first difference refers to second force
            last = x + criterion.window_size
            if last > l
                last = l
            end
            init_force = force_vector[x]
            for f in skipmissing(force_vector[(x + 1):last])
                if f > init_force + criterion.min_increase_in_window
                    return x
                end
            end
        end
    end
    return -1
end

function response_offset(
    force_vector::AbstractVector{<:FloatOrMissing},
    criterion::OnsetCriterion;
    onset::Integer,
)::Int
    # returns sample in vector with offset

    if onset < 0
        return -1
    else
        # find approaching first crossing initial force
        l = length(force_vector)
        from = onset + criterion.window_size + 1
        initial_force = force_vector[onset]
        for x in range(from, l)
            f = force_vector[x]
            if !ismissing(f) && f <= initial_force
                return x
            end
        end
        # reverse data and find onset
        rev_onset = response_onset(force_vector[l:-1:from], criterion)
        if rev_onset >= 0
            return l - rev_onset
        else
            return -1
        end
    end
end

function latency(
    rb::ResponseBounds;
    sampling_rate::Int,
    zero_sample::Int,
)::FloatOrMissing
    # returns response latency in millisecond
    if rb.onset < 0
        return missing
    end
    return duration(rb.onset - zero_sample; sampling_rate)
end

function latency(
    fp::ForceProfiles,
    vrb::AbstractVector{ResponseBounds}
)::Vector{FloatOrMissing}
    n_profiles(fp) == length(vrb) || throw(ArgumentError(
        "Number of profiles and ResponseBounds don't match!"))

    @unpack zero_sample, sr = fp
    return [latency(rb; zero_sample, sampling_rate=sr) for rb in vrb]
end

function duration(
    rb::ResponseBounds;
    sampling_rate::Int,
)::FloatOrMissing
    # returns response latency in millisecond
    if rb.onset < 0 || rb.offset < 0
        return missing
    end
    return duration(rb.offset - rb.onset; sampling_rate)
end

function duration(
    fp::ForceProfiles,
    vrb::AbstractVector{ResponseBounds}
)::Vector{FloatOrMissing}
    n_profiles(fp) == length(vrb) || throw(ArgumentError(
        "Number of profiles and ResponseBounds don't match!"))
    return [duration(rb; sampling_rate=fp.sr) for rb in vrb]
end

function peak_force(
    force_vector::AbstractVector{<:FloatOrMissing},
    rb::ResponseBounds,
)
    # returns (peak force, n sample to peak)
    resp = _extract_response(force_vector, rb)
    if isnothing(resp)
        rtn = (missing, -1)
    else
        rtn = findmax(skipmissing(resp))
    end
    return (; peak=rtn[1], sample_to_peak=rtn[2])
end

function peak_force(
    fp::ForceProfiles,
    rb::AbstractVector{ResponseBounds},
)
    n_profiles(fp) == length(rb) || throw(ArgumentError(
        "Number of profiles and ResponseBounds don't match!"))

    return [peak_force(f, b) for (f, b) in zip(eachrow(fp.dat), rb)]
end


function impulse_size(
    force_vector::AbstractVector{<:FloatOrMissing},
    rb::ResponseBounds,
)   # TODO not yet tested
    resp = _extract_response(force_vector, rb)
    if !isnothing(resp)
        return sum(skipmissing(resp .- resp[1]))
    end
end

function impulse_size(
    fp::ForceProfiles,
    rb::AbstractVector{ResponseBounds},
)
    n_profiles(fp) == length(rb) || throw(ArgumentError(
        "Number of profiles and ResponseBounds don't match!"))

    return [impulse_size(f, b) for (f, b) in zip(eachrow(fp.dat), rb)]
end

# helper
function _extract_response(
    force_vector::AbstractVector{<:FloatOrMissing},
    rb::ResponseBounds
)
    if rb.offset < 0 || rb.offset < 0
        return nothing
    end
    return force_vector[(rb.onset):(rb.offset)]
end;