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

TOptionalRowSelect = Union{Nothing, BitVector, Vector{<:Integer}}
function aggregate(fe::BeForEpochs,
	agg_fnc::Function,
	rows::TOptionalRowSelect,
	design::Union{DataFrame, DataFrameRow})

	if isnothing(rows)
		dat = agg_fnc(fe.dat, dims = 1)
		bsl = agg_fnc(fe.baseline)
	else
		dat = agg_fnc(fe.dat[rows, :], dims = 1)
		bsl = agg_fnc(fe.baseline[rows, :])
	end

	if design isa DataFrameRow
		design = DataFrame(design)
	elseif nrow(design) > 1
		throw(ArgumentError("Design has to be a DataFrame with one or no rows or DataFrameRow."))
	end

	if !(bsl isa Vector)
		bsl = [bsl]
	end
	meta = copy(fe.meta)
	meta["agg_fnc"] = string(agg_fnc)
	return BeForEpochs(dat, fe.sampling_rate;
		baseline = bsl,
		zero_sample = fe.zero_sample,
		design, meta)
end

aggregate(fe::BeForEpochs, agg_fnc::Function, rows::TOptionalRowSelect) = aggregate(fe, agg_fnc, rows, DataFrame())
aggregate(fe::BeForEpochs, agg_fnc::Function, design::Union{DataFrame, DataFrameRow}) = aggregate(fe, agg_fnc, nothing, design)

"""
	aggregate(fe::BeForEpochs, agg_fnc::Function rows::Union{Nothing, BitVector, Vector{<:Integer}}, design::Union{DataFrame, DataFrameRow})
	aggregate(fe::BeForEpochs, agg_fnc::Function, rows::Union{Nothing, BitVector, Vector{<:Integer}})
	aggregate(fe::BeForEpochs, agg_fnc::Function, design::Union{DataFrame, DataFrameRow})
	aggregate(fe::BeForEpochs, agg_fnc::Function;
			condition::ColumnIndex = :all,
			subject_id::Union{Nothing, ColumnIndex} = nothing)
TODO
"""
function aggregate(
	fe::BeForEpochs,
	agg_fnc::Function;
	condition::Union{Nothing, ColumnIndex} = :all, # condition column name
	subject_id::Union{Nothing, ColumnIndex} = nothing,
)
	if isnothing(condition)
		condition = :all
	end

	if condition == :all
		if isnothing(subject_id)
			return aggregate(fe, agg_fnc, nothing)
		end
		conditions = nothing
	else
		conditions = fe.design[:, condition]
	end
	rtn_array = BeForEpochs[]
	# aggregate per subject
	if isnothing(subject_id)
		for cond in unique(conditions)
			ids = findall(conditions .== cond)
			design = DataFrame(condition => cond)
			push!(rtn_array, aggregate(fe, agg_fnc, ids, design))
		end
	else
		subject_ids = fe.design[:, subject_id]
		for sid in unique(subject_ids)
			subj_rows = subject_ids .== sid
			if isnothing(conditions)
				design = DataFrame(:subject_id => sid)
				push!(rtn_array, aggregate(fe, agg_fnc, subj_rows, design))
			else
				for cond in unique(conditions)
					design = DataFrame(:subject_id => sid, condition => cond)
					ids = findall(subj_rows .&& conditions .== cond)
					push!(rtn_array, aggregate(fe, agg_fnc, ids, design))
				end
			end
		end
	end
	rtn = reduce(vcat, rtn_array)
	rtn.meta["agg_fnc"] = string(agg_fnc)
	return rtn
end;
