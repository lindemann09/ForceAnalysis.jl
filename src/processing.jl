peak_differences(fe::ForceEpochs, window_size::Integer) = peak_differences(fe.dat, window_size)

function peak_differences(force_mtx::Matrix{T}, window_size::Integer) where T <: AbstractFloat
	# peak difference per row
	nc = size(force_mtx, 2)
	rtn = T[]
	for row in eachrow(force_mtx)
		peak = 0.0
		for i in 1:(nc-window_size)
			@inbounds diff = abs(row[i+window_size] - row[i])
			if peak < diff
				peak = diff
			end
		end
		push!(rtn, peak)
	end
	return rtn
end;

function epoch_rejection(fe::ForceEpochs;
	force_range::UnitRange, # criteria for good trial
	max_difference, # criteria for good trial
	max_diff_windows_size
)
	min = minimum(fe)
	max = maximum(fe)
	peak_diff = peak_differences(fe.dat, max_diff_windows_size)
	good = min .>= force_range.start .&&
		  max .<= force_range.stop .&&
		  abs.(peak_diff) .<= max_difference
	return subset(fe, good)
end


function aggregate(
	# TODO generate methods with multiple IVs
	fe::ForceEpochs{T};
	condition::ColumnIndex = :all,
	subject_id::Union{Nothing, ColumnIndex} = nothing,
	agg_fnc = mean,
) where T <: AbstractFloat

	# aggregate per subject
	agg_forces = Matrix{T}(undef, 0, size(fe.dat, 2))
	agg_baseline = T[]

	if condition == :all
		conditions = repeat([true], nrow(fe.design))
	else
		conditions = fe.design[:, condition]
	end
	Tiv = eltype(unique(conditions)) # unique require to deal with CategoricalArrays

	if isnothing(subject_id)
		dsgn = Dict(condition => Tiv[])
		for cond in unique(conditions)
			push!(dsgn[condition], cond)
			ids = findall(conditions .== cond)
			agg_fe = agg_fnc(fe; rows = ids)
			agg_forces = vcat(agg_forces, agg_fe.dat)
			append!(agg_baseline, agg_fe.baseline)

		end
	else
		subject_ids = fe.design[:, subject_id]
		Tsid = eltype(unique(subject_ids))
		dsgn = Dict(condition => Tiv[], subject_id => Tsid[])
		for sid in unique(subject_ids)
			for cond in unique(conditions)
				push!(dsgn[condition], cond)
				push!(dsgn[subject_id], sid)
				ids = findall(subject_ids .== sid .&& conditions .== cond)
				agg_fe = agg_fnc(fe; rows = ids)
				agg_forces = vcat(agg_forces, agg_fe.dat)
				append!(agg_baseline, agg_fe.baseline)
			end
		end
	end

	delete!(dsgn, :all)
	return ForceEpochs(
		agg_forces, fe.sr, DataFrame(dsgn), agg_baseline, fe.zero_sample)
end;

function subset(fe::ForceEpochs, rows::Base.AbstractVecOrTuple{Integer})
	force = fe.dat[rows, :]
	bsln = fe.baseline[rows]
	subset_design = fe.design[rows, :]
	return ForceEpochs(force, fe.sr, subset_design, bsln, fe.zero_sample)
end


function subset(fe::ForceEpochs, args...)
	df = copy(fe.design)
	df.row_xxx .= 1:nrow(df)
	df = subset(df, args...)
	return subset(fe, df[:, :row_xxx])
end

