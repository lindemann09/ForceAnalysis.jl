const ColumnIndex = Union{Symbol, AbstractString}

peak_differences(fe::BeForEpochs, window_size::Integer) = peak_differences(fe.dat, window_size)

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

"""
	epoch_rejection_ids(fe::BeForEpochs; force_range::UnitRange, max_difference::Integer,
		max_diff_windows_size::Integer)

Return BitVector indicating 'good' epochs
"""
function epoch_rejection_ids(fe::BeForEpochs;
	force_range::UnitRange, # criteria for good trial
	max_difference::Integer, # criteria for good trial
	max_diff_windows_size::Integer,
)
	min = minimum(fe)
	max = maximum(fe)
	peak_diff = peak_differences(fe.dat, max_diff_windows_size)
	good = min .>= force_range.start .&&
		   max .<= force_range.stop .&&
		   abs.(peak_diff) .<= max_difference
	return good
end

"""
	epoch_rejection(fe::BeForEpochs; force_range::UnitRange, max_difference::Integer,
		max_diff_windows_size::Integer)

Return BeForEpochs fitting the criteria for 'good' epochs
"""
function epoch_rejection(fe::BeForEpochs;
	force_range::UnitRange, # criteria for good trial
	max_difference::Integer, # criteria for good trial
	max_diff_windows_size::Integer,
)
	good = epoch_rejection_ids(fe; force_range, max_difference, max_diff_windows_size)
	return subset(fe, good)
end


"""
	aggregate(fe::BeForEpochs; condition::ColumnIndex = :all,
		subject_id::Union{Nothing, ColumnIndex} = nothing, agg_fnc::Function = mean)
TODO
"""
function aggregate(
	# TODO generate methods with multiple IVs
	fe::BeForEpochs;
	condition::ColumnIndex = :all,
	subject_id::Union{Nothing, ColumnIndex} = nothing,
	agg_fnc::Function = mean,
)
	# aggregate per subject
	agg_forces = Matrix{Float64}(undef, 0, size(fe.dat, 2))
	agg_baseline = Float64[]

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
	return BeForEpochs(
		agg_forces, fe.sampling_rate, DataFrame(dsgn), agg_baseline, fe.zero_sample)
end;


"""
	subset(fe::BeForEpochs, rows::Base.AbstractVecOrTuple{Integer})
	subset(fe::BeForEpochs, args...)

TODO
"""
function DataFrames.subset(fe::BeForEpochs, rows::Base.AbstractVecOrTuple{Integer})
	force = fe.dat[rows, :]
	bsln = fe.baseline[rows]
	subset_design = fe.design[rows, :]
	return BeForEpochs(force, fe.sampling_rate, subset_design, bsln, fe.zero_sample)
end

function DataFrames.subset(fe::BeForEpochs, args...)
	df = copy(fe.design)
	df.row_xxx .= 1:nrow(df)
	df = subset(df, args...)
	return subset(fe, df[:, :row_xxx])
end