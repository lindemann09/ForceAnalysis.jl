#### ForceData
"""
	ForceData{T <: AbstractFloat}

TODO
"""
struct ForceData{T <: AbstractFloat}
	dat::Vector{T} # force data
	ts::Vector{Int} # timestamps
	sampling_rate::Float64
	meta::Dict

	function ForceData(dat::Vector{T}, ts::Vector{Int}, sampling_rate::Real, meta::Dict) where {T}
		check_timestamps(ts, length(dat))
		return new{T}(dat, ts, sampling_rate, meta)
	end
end;


ForceData(bfr::BeForRecord, force_column::String) =
	return ForceData(bfr.dat[!, force_column],  collect(bfr.time_stamps), bfr.sampling_rate, bfr.meta)

Base.propertynames(::ForceData) = (:dat, :timestamps, :sampling_rate, :meta, :n_samples)
function Base.getproperty(fd::ForceData, s::Symbol)
	if s === :timestamps
		return fd.ts
	elseif s === :n_samples
		return length(fd.dat)
	else
		return getfield(fd, s)
	end
end

"""
	force(fd::ForceData)
	force(fe::BeForEpochs)

Returns the force data as `Matrix` or `Vector`.
"""
force(fd::ForceData) = fd.dat

function Base.show(io::IO, mime::MIME"text/plain", x::ForceData)
	println(io, "ForceData")
	print(io, " $(x.n_samples) samples, sampling rate: $(x.sampling_rate)")
end;


force(fe::BeForEpochs) = fe.dat


# helper
function check_timestamps(timestamps::Vector, required_length::Int)
	lt = length(timestamps)
	return lt == 0 || lt == required_length ||
		   throw(
			   ArgumentError(
				   "if timestamps are defined, they must have the same length as samples of one sensor.",
			   ),
		   )
end
