function lowpass_filter(force_data::AbstractVector{T};
                    sampling_rate::Integer=1000,
                    cutoff_feq::Integer= 15,
                    butterworth_order::Integer = 4) where T<:Float64OrMissing
    responsetype = Lowpass(cutoff_feq, fs=sampling_rate)
    myfilter = digitalfilter(responsetype, Butterworth(butterworth_order))
    return filtfilt(myfilter, force_data .- force_data[1]) .+ force_data[1]  # filter centered data
end;

function force_profile_matrix(;force::AbstractVector{T},
                    time_stamps::AbstractVector{<:Integer},
                    zero_times::AbstractVector{<:Integer},
                    n_samples::Integer,
                    n_samples_before::Integer) where {T<:Float64OrMissing}

    len_force = length(force)
    if length(time_stamps) != len_force
        throw(error("forces and times_stamps have to have the same length"))
    end
    nrow = length(zero_times)
    rtn = Matrix{Float64OrMissing}(missing, nrow, n_samples_before + n_samples)
    for r in 1:nrow
        i = _find_larger_or_equal(zero_times[r], time_stamps)
        if i !== nothing
            from = (i-n_samples_before)
            to = (i+n_samples-1)
            if from<len_force
                if to>len_force
                    to = len_force
                end
                rtn[r, 1:(to-from+1)] .= force[from:to]
            end
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
                    butterworth_order::Integer = 4) where T<:Float64OrMissing
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
    rtn = ForceProfiles(force_mtx, DataFrame(), bsl, n_samples_before+1)
    return rtn
end;

function peak_difference(force_mtx::Matrix{<:Float64OrMissing};
                            window_size::Integer = 100)
    # peak difference per row
    (nr, nc) = size(force_mtx)
    peak = Vector{Float64OrMissing}(missing, nr)
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
                force_profile_matrix::Matrix{<:Float64OrMissing};
                force_range::UnitRange = -400:400,# criteria for good trial
                max_difference = 200, # criteria for good trial
                max_diff_windows_size = 100)
    # force profile quality parameter
    (min, max) = row_minmax(force_profile_matrix)
    df = DataFrame(min = min[:],
        max = max[:],
        peak_differences = peak_difference(force_profile_matrix;
                                    window_size=max_diff_windows_size))

    tmp = df.min.>force_range.start .&&
                    df.max.<force_range.stop .&&
                    abs.(df.peak_differences) .< max_difference
    df.good_trial = convert(Vector{Bool}, tmp) # TODO why convert
    return df
end;


function aggregate(fp::ForceProfiles; dv::Symbol,
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
    return ForceProfiles(agg_forces, agg_df, eltype(agg_forces)[],  fp.zero_sample)
end;


function aggregate(aggregate_fnc::Function, fp::ForceProfiles,
                    grouped_design::GroupedDataFrame; var_row = :row)
    # aggregate per subject
    agg_forces = Matrix{Float64}(undef, 0, size(fp.force, 2))
    agg_df = DataFrame()
    cols = groupcols(grouped_design)
    for sub_df in grouped_design
        append!(agg_df, DataFrame(sub_df[1, cols]))
        rows = sub_df[:, var_row]
        m = aggregate_fnc(fp.force[rows, :])
        agg_forces = vcat(agg_forces, transpose(m)) # TODO hcat transpose later?
    end
    return ForceProfiles(agg_forces, agg_df, eltype(agg_forces)[],  fp.zero_sample)
end;

aggregate(fp::ForceProfiles, grouped_design::GroupedDataFrame; var_row = :row) =
    aggregate(column_mean, fp, grouped_design; var_row)


function subset(fp::ForceProfiles, subset_design::DataFrame)
    i = subset_design.row
    mtx = fp.force[i,:]
    subset_design.row = 1:nrow(subset_design) # renumber
    return ForceProfiles(mtx, subset_design, eltype(mtx)[], fp.zero_sample)
end

### helper functions
function _find_larger_or_equal(needle::T, sorted_array::AbstractVector{T}) where T<:Real
    cnt::Int = 0
    for x in sorted_array
        cnt = cnt + 1
        if x >= needle
            return cnt
        end
    end
    return nothing
end;
