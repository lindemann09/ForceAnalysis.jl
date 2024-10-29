const ColumnIndex = Union{Symbol, AbstractString}
# const MultiColumnIndex = Union{
#     ColumnIndex, AbstractVector{<:ColumnIndex},Tuple{<:ColumnIndex}
# }

"""
	epochs(fd::ForceData{<: AbstractFloat};
		zero_times::AbstractVector{<:Integer},
		n_samples::Integer,	n_samples_before::Integer)

TODO
"""
function epochs(
	fd::ForceData{T};
	zero_times::AbstractVector{<:Integer},
	n_samples::Integer,
	n_samples_before::Integer,
	design::Union{Nothing, DataFrame}=nothing
) where T <: AbstractFloat
	@unpack dat, ts = fd
	samples_fd = fd.n_samples # samples for data
	n_epochs = length(zero_times)
	ncol = n_samples_before + n_samples
	force_mtx = Matrix{T}(undef, n_epochs, ncol)
	for (r, zt) in enumerate(zero_times)
		i = _find_larger_or_equal(zt, ts)
		if i !== nothing
			from = i - n_samples_before
			if from < samples_fd
				to = i + n_samples - 1
				if to > samples_fd
					@warn string("extract_force_epochs: last force epoch is incomplete, ",
								to-samples_fd, " samples missing.")
					force_mtx[r, :] .= vcat(dat[from:samples_fd], zeros(T, to - samples_fd))
				else
					force_mtx[r, :] .= dat[from:to]
				end
			end
		end
	end

	if isnothing(design)
		design = DataFrame()
	end
	return ForceEpochs(force_mtx, fd.sampling_rate, design,
		zeros(T, n_epochs), n_samples_before + 1)
end;

"""
	adjust_baseline!(fe::ForceEpochs, baseline_window::UnitRange{<:Integer})

TODO
"""
function adjust_baseline!(fe::ForceEpochs, baseline_window::UnitRange{<:Integer})
	dat = fe.dat .+ fe.baseline
	bsl = mean(dat[:, baseline_window], dims = 2)
	fe.dat[:, :] .= dat .- bsl
	fe.baseline[:] .= bsl
	return fe
end


### helper functions
function _find_larger_or_equal(needle::T, sorted_array::AbstractVector{T}) where {T <: Real}
	cnt::Int = 0
	for x in sorted_array
		cnt = cnt + 1
		if x >= needle
			return cnt
		end
	end
	return nothing
end;
