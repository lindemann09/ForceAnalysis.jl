
function lowpass_filter(dat::Vector{<:FloatOrMissing};
    sampling_rate::Integer,
    cutoff_feq::Integer,
    butterworth_order::Integer
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
    force_data::ForceData;
    zero_times::AbstractVector{<:Integer},
    n_samples::Integer,
    n_samples_before::Integer
)
    @unpack dat, ts = force_data
    len_force = length(dat)
    nrow = length(zero_times)
    rtn = Matrix{FloatOrMissing}(missing, nrow, n_samples_before + n_samples)
    for r in 1:nrow
        i = _find_larger_or_equal(zero_times[r], ts)
        if i !== nothing
            from = (i - n_samples_before)
            to = (i + n_samples - 1)
            if from < len_force
                if to > len_force
                    to = len_force
                end
                rtn[r, 1:(to-from+1)] .= dat[from:to]
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
    butterworth_order::Integer=4
)
    ## filter and scale data
    flt = lowpass_filter(force_data.dat; sampling_rate=force_data.sr,
        cutoff_feq=filter_cutoff_feq, butterworth_order) .* scale_forces
    tmp = ForceData(flt, force_data.ts, force_data.sr, force_data.meta)
    # extract force profile per trial
    force_mtx = force_profile_matrix(tmp;
        zero_times=profiles_zero_times, n_samples, n_samples_before)
    #baseline adjustment
    bsl = row_mean(force_mtx[:, baseline_sample_range])
    force_mtx = force_mtx .- bsl
    return ForceProfiles(
        force_mtx, sampling_rate(force_data), DataFrame(), bsl, n_samples_before + 1)
end;

function peak_difference(force_mtx::Matrix{<:FloatOrMissing};
    window_size::Integer=100
)
    # peak difference per row
    (nr, nc) = size(force_mtx)
    peak = Vector{FloatOrMissing}(missing, nr)
    for r in 1:nr
        for i in 1:(nc-window_size)
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

function profile_parameter(fp::ForceProfiles;
    force_range::UnitRange=-400:400, # criteria for good trial
    max_difference=200, # criteria for good trial
    max_diff_windows_size=100
)
    force_profile_matrix = force(fp)
    # force profile quality parameter
    (min, max) = row_minmax(force_profile_matrix)
    df = DataFrame(;
        min=min[:],
        max=max[:],
        peak_differences=peak_difference(force_profile_matrix;
            window_size=max_diff_windows_size)
    )
    tmp = df.min .> force_range.start .&&
        df.max .< force_range.stop .&&
        abs.(df.peak_differences) .< max_difference
    df.good_trial = convert(Vector{Bool}, tmp)
    return df
end;

function aggregate(
    fp::ForceProfiles;
    iv::Symbol,
    subject_id=:subject_id,
    row_column=:row,
    agg_fnc=column_mean
)
    # aggregate per subject
    agg_forces = Matrix{Float64}(undef, 0, size(fp.force, 2))
    agg_baseline = Float64[]
    rtn = Dict(iv => [], subject_id => [])
    design = fp.design
    bsln = hcat(fp.baseline) # convert to nx1 matrix
    for sid in unique(design[:, subject_id])
        for condition in unique(design[:, iv])
            append!(rtn[iv], condition)
            append!(rtn[subject_id], sid)
            ids = findall(design[:, subject_id] .== sid .&& design[:, iv] .== condition)
            rows = design[ids, row_column]
            m = agg_fnc(fp.force; rows) #<===============
            agg_forces = vcat(agg_forces, transpose(m))
            m = agg_fnc(bsln; rows) #<===============
            append!(agg_baseline, m)

        end
    end
    rtn[row_column] = 1:length(rtn[iv])
    return ForceProfiles(
        agg_forces, fp.sr, DataFrame(rtn), agg_baseline, fp.zero_sample)
end;

# function aggregate(aggregate_fnc::Function, fp::ForceProfiles,
#                     grouped_design::GroupedDataFrame; var_row = :row)
#     # aggregate per subject
#     agg_forces = Matrix{Float64}(undef, 0, size(fp.force, 2))
#     agg_df = DataFrame()
#     cols = groupcols(grouped_design)
#     for sub_df in grouped_design
#         append!(agg_df, DataFrame(sub_df[1, cols]))
#         rows = sub_df[:, var_row]
#         m = aggregate_fnc(fp.force[rows, :])
#         agg_forces = vcat(agg_forces, transpose(m)) # TODO hcat transpose later?
#     end
#     return ForceProfiles(agg_forces, agg_df, eltype(agg_forces)[],  fp.zero_sample)
# end;

# aggregate(fp::ForceProfiles, grouped_design::GroupedDataFrame; var_row = :row) =
#     aggregate(column_mean, fp, grouped_design; var_row)

function subset(fp::ForceProfiles, subset_design::DataFrame)
    i = subset_design.row
    subset_design.row = 1:nrow(subset_design) # renumber
    return ForceProfiles(fp.force[i, :], fp.sr, subset_design,
                fp.baseline[i], fp.zero_sample)
end

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
