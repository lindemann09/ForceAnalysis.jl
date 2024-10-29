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
	force(fd::MultiForceData, id::Union{Int, String, Symbol})
	force(fe::ForceEpochs)

Returns the force data as `Matrix` or `Vector`.
"""
force(fd::ForceData) = fd.dat

function Base.show(io::IO, mime::MIME"text/plain", x::ForceData)
	println(io, "ForceData")
	print(io, " $(x.n_samples) samples, sampling rate: $(x.sampling_rate)")
end;



#### MultiForceData
"""
	MultiForceData{T <: AbstractFloat}

TODO
"""
struct MultiForceData{N, T <: AbstractFloat}
	dat::Matrix{T} # force data
	ts::Vector{Int} # timestamps
	sr::Float64
	ids::NTuple{N, Symbol} # ids/labels
	meta::Dict

	function MultiForceData(dat::Matrix{T}, ts::Vector{Int}, sr::Real,
		ids::Vector, meta::Dict) where {T}
		check_timestamps(ts, size(dat, 1))
		N = size(dat, 2)
		ls = length(ids)
		if ls == 0
			ids = [Symbol(Char(x)) for x in 97:(97+N-1)]
		else
			ls == size(dat, 2) || throw(
				ArgumentError(
					"if ids are defined, there must be one label for each sensor."),
			)
			ids = Symbol.(ids)
		end
		return new{N, T}(dat, ts, sr, NTuple{N, Symbol}(ids), meta)
	end
end;

function MultiForceData(
	df::DataFrame;
	sampling_rate::Real,
	time_stamp_col::Union{Nothing, Symbol, String} = nothing,
	meta::Dict = Dict(),
)
	# dataframe to MultiForceData
	if time_stamp_col isa Symbol
		time_stamp_col = String(time_stamp_col)
	end
	ts = Int[]
	ids = String[]
	for n in names(df)
		if n == time_stamp_col
			ts = df[:, n]
		else
			push!(ids, n)
		end
	end
	return MultiForceData(Matrix(df[:, ids]), ts, sampling_rate, ids, meta)
end

Base.propertynames(::MultiForceData) = (:dat, :timestamps, :sampling_rate, :ids,
	:meta, :n_samples)

function Base.getproperty(x::MultiForceData, s::Symbol)
	if s === :timestamps
		return x.ts
	elseif s === :sampling_rate
		return x.sampling_rate
	elseif s === :n_samples
		return size(x.dat, 1)
	else
		return getfield(x, s)
	end
end

function force(fd::MultiForceData, id::Union{Int,String,Symbol})
	if id isa Int
		return fd.dat[:, id]
	else
		# returns force
		idx = findfirst(fd.ids .== Symbol(id))
		idx isa Integer || throw(ArgumentError(
			"Can not find force data labelled '$(id)'"))
		return force(fd, idx)
	end
end


function ForceData(fd::MultiForceData, id::Union{Integer, Symbol, String})
	return ForceData(force(fd, id), fd.ts, fd.sampling_rate, fd.meta)
end

function Base.show(io::IO, mime::MIME"text/plain", x::MultiForceData)
	println(io, "MultiForceData")
	print(io, "  $(x.n_samples) samples, sampling rate: $(x.sampling_rate)")
end;


############  ForceEpochs

"""
	ForceEpochs{T <: AbstractFloat}

TODO
"""
struct ForceEpochs{T <: AbstractFloat}
	dat::Matrix{T}
	sampling_rate::Float64
	design::DataFrame
	baseline::Vector{T}
	zero_sample::Int

	function ForceEpochs(force::Matrix{T}, sampling_rate::Real, design::DataFrame,
		baseline::Vector{T}, zero_sample::Int) where {T <: AbstractFloat}
		lf = size(force, 1)
		lb = length(baseline)
		lf == lb || throw(
			ArgumentError(
				"Number of rows of force ($(lf)) must match the length of baseline ($(lb)).",
			),
		)
		n = nrow(design)
		n == 0 || lf == n || throw(
			ArgumentError(
				"Number of rows of force ($(lf)) must match the number of rows ins the design ($(n)).",
			),
		)
		return new{T}(force, sampling_rate, design, baseline, zero_sample::Int)
	end
end;


Base.propertynames(::ForceEpochs) = (:dat, :sampling_rate, :design, :baseline,
	:zero_sample, :n_samples, :n_epochs)
function Base.getproperty(x::ForceEpochs, s::Symbol)
	if s === :n_epochs
		return size(x.dat, 1)
	elseif s === :n_samples
		return size(x.dat, 2)
	else
		return getfield(x, s)
	end
end

"""
	copy(fe::ForceEpochs)

TODO
"""
function Base.copy(fe::ForceEpochs)
	return ForceEpochs(copy(fe.dat), fe.sampling_rate, copy(fe.design),
		copy(fe.baseline), fe.zero_sample)
end

force(fe::ForceEpochs) = fe.dat


function Base.show(io::IO, mime::MIME"text/plain", x::ForceEpochs)
	println(io, "ForceEpochs")
	println(io, "  $(x.n_epochs) epochs")
	println(io, "  $(x.n_samples) samples, sampling rate: $(x.sampling_rate), zero sample: $(x.zero_sample)")
	print(io, "  Design: $(names(x.design))")
end


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
