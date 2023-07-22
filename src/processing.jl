function peak_difference(force_mtx::Matrix{<:T};
	window_size::Integer = 100,
) where T<:AbstractFloat
	# peak difference per row
	(nr, nc) = size(force_mtx)
	peak = Vector{T}(undef, nr)
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
	force_range::UnitRange = -400:400, # criteria for good trial
	max_difference = 200, # criteria for good trial
	max_diff_windows_size = 100,
)
	force_profile_matrix = fp.dat
	# force profile quality parameter
	df = DataFrame(;
		min = minimum(fp),
		max = maximum(fp),
		peak_differences = peak_difference(force_profile_matrix;
			window_size = max_diff_windows_size),
	)
	tmp =
		df.min .> force_range.start .&&
		df.max .< force_range.stop .&&
		abs.(df.peak_differences) .< max_difference
	df.good_trial = convert(Vector{Bool}, tmp)
	return df
end;

function aggregate(
	# TODO generate methods with multiple IVs
	fp::ForceProfiles{T};
	condition::ColumnIndex = :all,
	subject_id::Union{Nothing, ColumnIndex} = nothing,
	row_idx_column = :row,
	agg_fnc = mean,
) where T<:AbstractFloat

	# aggregate per subject
	agg_forces = Matrix{T}(undef, 0, size(fp.dat, 2))
	agg_baseline = T[]
	rows = fp.design[:, row_idx_column]
	if condition == :all
		conditions = repeat([true], nrow(fp.design))
	else
		conditions = fp.design[:, condition]
	end
	Tiv = eltype(unique(conditions)) # unique require to deal with CategoricalArrays
	bsln = hcat(fp.baseline) # convert to nx1 matrix

	if isnothing(subject_id)
		dsgn = Dict(condition => Tiv[], row_idx_column => Int64[])
		for cond in unique(conditions)
			push!(dsgn[condition], cond)
			ids = findall(conditions .== cond)
			agg_fp = agg_fnc(fp; rows=rows[ids])
            agg_forces = vcat(agg_forces, agg_fp.dat) #FIXME check transpose?
            append!(agg_baseline, agg_fp.bsl)

		end
	else
		subject_ids = fp.design[:, subject_id]
		Tsid = eltype(unique(subject_ids))
		dsgn = Dict(condition => Tiv[], subject_id => Tsid[],
			row_idx_column => Int64[])
		for sid in unique(subject_ids)
			for cond in unique(conditions)
				push!(dsgn[condition], cond)
				push!(dsgn[subject_id], sid)
				ids = findall(subject_ids .== sid .&& conditions .== cond)
				agg_fp = agg_fnc(fp; rows=rows[ids])
				agg_forces = vcat(agg_forces, agg_fp.dat)
				append!(agg_baseline, agg_fp.baseline)
				end
		end
	end

	dsgn[row_idx_column] = 1:length(dsgn[condition])
	delete!(dsgn, :all)
	return ForceProfiles(
		agg_forces, fp.sr, DataFrame(dsgn), agg_baseline, fp.zero_sample)
end;

function subset(fp::ForceProfiles, rows::Base.AbstractVecOrTuple{Integer}; row_idx_column::String = "row")
	force = fp.dat[rows, :]
	bsln = fp.baseline[rows]
	subset_design = fp.design[rows, :]
	subset_design[:, row_idx_column] = 1:nrow(subset_design) # renumber
	return ForceProfiles(force, fp.sr, subset_design, bsln, fp.zero_sample)
end

function subset(fp::ForceProfiles, subset_design::DataFrame; row_idx_column::String = "row")
	rows = subset_design[:, row_idx_column]
	return subset(fp, convert(Vector{Int64}, rows); row_idx_column)
end

function subset(fp::ForceProfiles, args...; row_idx_column::String = "row")
	df = subset(fp.design, args...)
	rows = convert(Vector{Int64}, df[:, row_idx_column])
	return subset(fp, rows; row_idx_column)
end

