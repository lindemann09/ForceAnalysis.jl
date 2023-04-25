function lowpass_filter(force_data::AbstractVector{T};
                    sampling_rate::Integer=1000,
                    cutoff_feq::Integer= 15,
                    butterworth_order::Integer = 4) where T<:FloatOrMissing
    responsetype = Lowpass(cutoff_feq, fs=sampling_rate)
    myfilter = digitalfilter(responsetype, Butterworth(butterworth_order))
    return filtfilt(myfilter, force_data .- force_data[1]) .+ force_data[1]  # filter centered data
end;

function force_profile_matrix(;force::AbstractVector{T},
                    time_stamps::AbstractVector{<:Integer},
                    zero_times::AbstractVector{<:Integer},
                    n_samples::Integer,
                    n_samples_before::Integer) where {T<:FloatOrMissing}

    if length(time_stamps) != length(force)
        throw(error("forces and times_stamps have to have the same length"))
    end
    nrow = length(zero_times)
    rtn = Matrix{FloatOrMissing}(missing, nrow, n_samples_before + n_samples)
    for r in 1:nrow
        i = _find_larger_or_equal(zero_times[r], time_stamps)
        if i !== nothing
            rtn[r,:] .= force[(i-n_samples_before):(i+n_samples-1)]
        end
    end
    return rtn
end;

function force_data_preprocess(;
                    force::AbstractVector{T},
                    time_stamps::AbstractVector{<:Integer},
                    sampling_rate::Integer,
                    profiles_zero_times::AbstractVector{<:Integer},
                    n_samples::Integer,
                    n_samples_before::Integer,
                    baseline_sample_range::UnitRange{<:Integer},
                    scale_forces::AbstractFloat = 1,
                    filter_cutoff_feq::Integer = 15,
                    butterworth_order::Integer = 4) where T<:FloatOrMissing
    ## convenience function
    force = lowpass_filter(force;
            sampling_rate = sampling_rate,
            cutoff_feq = filter_cutoff_feq,
            butterworth_order) .* scale_forces; # filter and scale data

    # extract force profile per trial
    force_mtx =  force_profile_matrix(; force, time_stamps,
                        zero_times=profiles_zero_times,
                        n_samples, n_samples_before)
    #baseline adjustment
    bsl = row_mean(force_mtx[:, baseline_sample_range]);
    force_mtx = force_mtx .- bsl
    return ForceProfiles(force_mtx, DataFrame(), bsl, n_samples_before+1)
end;

function peak_difference(force_mtx::Matrix{<:FloatOrMissing};
                            window_size::Integer = 100)
    # peak difference per row
    (nr, nc) = size(force_mtx)
    peak = Vector{FloatOrMissing}(missing, nr)
    for r in 1:nr
        for i in 1:nc-window_size
            diff = abs(force_mtx[r, i+window_size] - force_mtx[r, i])
            if ismissing(peak[r])
                peak[r] = diff
            elseif peak[r] < diff
                peak[r] = diff
            end
        end
    end
    return peak
end;

function profile_parameter(force_profiles::ForceProfiles;
                force_range::UnitRange = -400:400, # criteria for good trial
                max_difference = 200, # criteria for good trial
                max_diff_windows_size = 100)
    return profile_parameter(force_profiles.force; force_range, max_difference,
                max_diff_windows_size)
end;

function profile_parameter(
                force_profile_matrix::Matrix{<:FloatOrMissing};
                force_range::UnitRange = -400:400,# criteria for good trial
                max_difference = 200, # criteria for good trial
                max_diff_windows_size = 100)
    # force profile quality parameter
    (min, max) = row_minmax(force_profile_matrix)
    df = DataFrame(min = min[:],
        max = max[:],
        peak_differences = peak_difference(force_profile_matrix;
                                    window_size=max_diff_windows_size))

    df.good_trial = df.min.>force_range.start .&&
                    df.max.<force_range.stop .&&
                    abs.(df.peak_differences) .< max_difference
    return df
end;


function aggregate_force_profiles(fp::ForceProfiles; dv::Symbol,
                        var_sub_id = :subject_id,
                        var_row = :row)
    # aggregate per subject
    agg_forces = Matrix{Float64}(undef, 0, size(fp.force, 2))
    agg_df = DataFrame()
    design = fp.design
    for sid in unique(design[:, var_sub_id])
        for condition in unique(design[:, dv])
            append!(agg_df, DataFrame(; condition, subject_id=sid))
            ids = findall(design[:, var_sub_id] .== sid .&&
                          design[:, dv] .== condition)
            rows = design[ids, var_row]
            m = column_mean(fp.force; rows)
            agg_forces = vcat(agg_forces, transpose(m))
        end
    end;
    #agg_df.row = 1:nrow(agg_df)
    return ForceProfiles(agg_forces, agg_df, Float64[],  fp.zero_sample)
end;


### helper functions
function _find_larger_or_equal(needle::T, sorted_array::AbstractVector{T}) where T<:Real
    cnt::UInt64 = 0
    for x in sorted_array
        cnt = cnt + 1
        if x >= needle
            return cnt
        end
    end
    return nothing
end;
