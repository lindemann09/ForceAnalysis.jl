const ColumnIndex = Union{Symbol,AbstractString}
# const MultiColumnIndex = Union{
#     ColumnIndex, AbstractVector{<:ColumnIndex},Tuple{<:ColumnIndex}
# }

function lowpass_filter(dat::Vector{<:AbstractFloat};
    sampling_rate::Integer,
    cutoff_feq::Integer,
    butterworth_order::Integer,
)
    responsetype = Lowpass(cutoff_feq; fs=sampling_rate)
    myfilter = digitalfilter(responsetype, Butterworth(butterworth_order))
    return filtfilt(myfilter, dat .- dat[1]) .+ dat[1]  # filter centered data
end;

function lowpass_filter(force_data::ForceData;
    cutoff_feq::Integer,
    butterworth_order::Integer=4
)
    @unpack dat, sr = force_data
    flt = lowpass_filter(dat; sampling_rate=sr, cutoff_feq, butterworth_order)
    return ForceData(flt, force_data.ts, sr, force_data.meta)
end;

function lowpass_filter(force_data::MultiForceData;
    cutoff_feq::Integer,
    butterworth_order::Integer=4
)
    @unpack dat, sr = force_data
    flt = similar(dat)
    for (i, d) in enumerate(eachcol(dat))
        flt[i,] = lowpass_filter(d; sampling_rate=sr, cutoff_feq, butterworth_order)
    end
    return MultiForceData(flt, force_data.ts, sr, force_data.ids, force_data.meta)
end;

function force_profile_matrix(
    force_data::ForceData{T};
    zero_times::AbstractVector{<:Integer},
    n_samples::Integer,
    n_samples_before::Integer,
) where T<:AbstractFloat
    @unpack dat, ts = force_data
    len_force = length(dat)
    nrow = length(zero_times)
    rtn = Matrix{T}(undef, nrow, n_samples_before + n_samples)
    for r in 1:nrow
        i = _find_larger_or_equal(zero_times[r], ts)
        if i !== nothing
            from = (i - n_samples_before)
            to = (i + n_samples - 1)
            if from < len_force
                if to > len_force
                    to = len_force
                end
                rtn[r, 1:(to - from + 1)] .= dat[from:to]
            end
        end
    end
    return rtn
end;

function force_data_preprocess(force_data::ForceData;
    profiles_zero_times::AbstractVector{<:Integer},
    n_samples::Integer,
    n_samples_before::Integer,
    baseline_sample_range::UnitRange{<:Integer},
    scale_forces::AbstractFloat=1,
    filter_cutoff_feq::Integer=15,
    butterworth_order::Integer=4,
)
    ## filter and scale data
    flt =
        lowpass_filter(force_data.dat; sampling_rate=force_data.sr,
            cutoff_feq=filter_cutoff_feq, butterworth_order) .* scale_forces
    tmp = ForceData(flt, force_data.ts, force_data.sr, force_data.meta)
    # extract force profile per trial
    force_mtx = force_profile_matrix(tmp;
        zero_times=profiles_zero_times, n_samples, n_samples_before)
    #baseline adjustment (row_mean)
    bsl = vec(mean(force_mtx[:, baseline_sample_range], dims=2))
    force_mtx = force_mtx .- bsl
    return ForceProfiles(
        force_mtx, sampling_rate(force_data), DataFrame(), bsl, n_samples_before + 1)
end;


### helper functions
function _find_larger_or_equal(needle::T, sorted_array::AbstractVector{T}) where {T<:Real}
    cnt::Int = 0
    for x in sorted_array
        cnt = cnt + 1
        if x >= needle
            return cnt
        end
    end
    return nothing
end;
