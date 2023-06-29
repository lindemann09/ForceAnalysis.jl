function peak_difference(force_mtx::Matrix{<:FloatOrMissing};
    window_size::Integer=100,
)
    # peak difference per row
    (nr, nc) = size(force_mtx)
    peak = Vector{FloatOrMissing}(missing, nr)
    for r in 1:nr
        for i in 1:(nc - window_size)
            diff = abs(force_mtx[r, i + window_size] - force_mtx[r, i])
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
    max_diff_windows_size=100,
)
    force_profile_matrix = force(fp)
    # force profile quality parameter
    (min, max) = row_minmax(force_profile_matrix)
    df = DataFrame(;
        min=min[:],
        max=max[:],
        peak_differences=peak_difference(force_profile_matrix;
            window_size=max_diff_windows_size),
    )
    tmp =
        df.min .> force_range.start .&&
        df.max .< force_range.stop .&&
        abs.(df.peak_differences) .< max_difference
    df.good_trial = convert(Vector{Bool}, tmp)
    return df
end;

function aggregate( # FIXME generate methods with multiple variables
    fp::ForceProfiles;
    iv::ColumnIndex,
    subject_id::Union{Nothing,ColumnIndex}=nothing,
    row_idx_column=:row,
    agg_fnc=column_mean,
)
    # aggregate per subject
    agg_forces = Matrix{Float64}(undef, 0, size(fp.dat, 2))
    agg_baseline = Float64[]
    Tiv = eltype(fp.design[:, iv])
    design = fp.design
    bsln = hcat(fp.baseline) # convert to nx1 matrix
    if isnothing(subject_id)
        rtn = Dict(iv => Tiv[])
        for condition in unique(design[:, iv])
            push!(rtn[iv], condition)
            ids = findall(design[:, subject_id] .== sid .&& design[:, iv] .== condition)
            rows = design[ids, row_idx_column]
            m = agg_fnc(fp.dat; rows) #<===============
            agg_forces = vcat(agg_forces, transpose(m))
            m = agg_fnc(bsln; rows) #<===============
            append!(agg_baseline, m)
        end
    else
        Tsid = eltype(fp.design[:, subject_id])
        rtn = Dict(iv => Tiv[], subject_id => Tsid[])
        for sid in unique(design[:, subject_id])
            for condition in unique(design[:, iv])
                push!(rtn[iv], condition)
                push!(rtn[subject_id], sid)
                ids = findall(design[:, subject_id] .== sid .&& design[:, iv] .== condition)
                rows = design[ids, row_idx_column]
                m = agg_fnc(fp.dat; rows) #<===============
                agg_forces = vcat(agg_forces, transpose(m))
                m = agg_fnc(bsln; rows) #<===============
                append!(agg_baseline, m)
            end
        end
    end
    rtn[row_idx_column] = 1:length(rtn[iv])
    return ForceProfiles(
        agg_forces, fp.sr, DataFrame(rtn), agg_baseline, fp.zero_sample)
end;

function subset(fp::ForceProfiles, rows::Base.AbstractVecOrTuple{Integer}; row_idx_column::String="row")
    force = fp.dat[rows, :]
    bsln = fp.baseline[rows]
    subset_design = fp.design[rows, :]
    subset_design[:, row_idx_column] = 1:nrow(subset_design) # renumber
    return ForceProfiles(force, fp.sr, subset_design, bsln, fp.zero_sample)
end

function subset(fp::ForceProfiles, subset_design::DataFrame; row_idx_column::String="row")
    rows = subset_design[:, row_idx_column]
    return subset(fp, rows; row_idx_column)
end

function subset(fp::ForceProfiles, row::Integer; row_idx_column::String="row")
    return subset(fp, [row]; row_idx_column)
end


