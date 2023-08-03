const ColumnIndex = Union{Symbol, AbstractString}
# const MultiColumnIndex = Union{
#     ColumnIndex, AbstractVector{<:ColumnIndex},Tuple{<:ColumnIndex}
# }

function lowpass_filter(dat::AbstractVector{<:AbstractFloat};
	sampling_rate::Real,
	cutoff_freq::Real,
	butterworth_order::Integer,
)
	responsetype = Lowpass(cutoff_freq; fs = sampling_rate)
	myfilter = digitalfilter(responsetype, Butterworth(butterworth_order))
	return filtfilt(myfilter, dat .- dat[1]) .+ dat[1]  # filter centered data
end;

function lowpass_filter!(fd::ForceData;
	cutoff_freq::Real,
	butterworth_order::Integer = 4,
)
	fd.dat[:] = lowpass_filter(fd.dat;
		sampling_rate = fd.sr, cutoff_freq, butterworth_order)
	return fd
end;

function lowpass_filter!(fe::ForceEpochs{T};
	cutoff_freq::Real,
	butterworth_order::Integer = 4,
) where T <: AbstractFloat
	for i in 1:fe.n_epochs
		fe.dat[i, :] = lowpass_filter(vec(fe.dat[i, :]);
			sampling_rate = fe.sr, cutoff_freq, butterworth_order)
	end
	return fe
end;

function lowpass_filter!(fd::MultiForceData;
	cutoff_freq::Real,
	butterworth_order::Integer = 4,
)
	for (i, d) in enumerate(eachcol(fd.dat))
		fd.dat[i, :] = lowpass_filter(d;
			sampling_rate = fd.sr, cutoff_freq, butterworth_order)
	end
	return fd
end;

function extract_force_epochs(
	force_data::ForceData{T};
	zero_times::AbstractVector{<:Integer},
	n_samples::Integer,
	n_samples_before::Integer,
) where T <: AbstractFloat
	@unpack dat, ts = force_data
	len_force = length(dat)
	nrow = length(zero_times)
	ncol = n_samples_before + n_samples
	force_mtx = Matrix{T}(undef, nrow, ncol)
	for r in 1:nrow
		i = _find_larger_or_equal(zero_times[r], ts)
		if i !== nothing
			from = (i - n_samples_before)
			to = (i + n_samples - 1)
			if from < len_force
				if to > len_force
					@warn "extract_force_epochs: last force epoch is incomplete"
					force_mtx[r, :] .= vcat(dat[from:len_force], zeros(T, to - len_force))
				else
					force_mtx[r, :] .= dat[from:to]
				end
			end
		end
	end
	return ForceEpochs(force_mtx, force_data.sampling_rate,
		DataFrame(), zeros(T, nrow), n_samples_before + 1)
end;

function scale_force!(fd::ForceData, factor::AbstractFloat)
	fd.dat[:] = fd.dat .* factor
	return fd
end

function scale_force!(fe::ForceEpochs, factor::AbstractFloat)
	fe.dat[:, :] = fe.dat .* factor
	fe.baseline[:] = fe.baseline .* factor
	return fe
end

function adjust_baseline!(fe::ForceEpochs; sample_range::UnitRange{<:Integer})
	dat = fe.dat .+ fe.baseline
	bsl = mean(dat[:, sample_range], dims = 2)
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
