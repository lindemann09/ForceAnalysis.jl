struct ForceData{T <: AbstractFloat}
	dat::Vector{T} # force data
	ts::Vector{Int} # timestamps
	sr::Int # sampling_rate
	meta::Dict

	function ForceData(dat::Vector{T}, ts::Vector{Int}, sr::Int, meta::Dict) where {T}
		check_timestamps(ts, length(dat))
		return new{T}(dat, ts, sr, meta)
	end
end;

Base.propertynames(::ForceData) = (:dat, :timestamps, :sampling_rate, :meta, :n_samples)
function Base.getproperty(x::ForceData, s::Symbol)
	if s === :timestamps
		return x.ts
	elseif s === :sampling_rate
		return x.sr
	elseif s === :n_samples
		return length(x.dat)
	else
		return getfield(x, s)
	end
end


struct MultiForceData{N, T <: AbstractFloat}
	dat::Matrix{T} # force data
	ts::Vector{Int} # timestamps
	sr::Int # sampling rate
	ids::NTuple{N, Symbol} # ids/labels
	meta::Dict

	function MultiForceData(dat::Matrix{T}, ts::Vector{Int}, sr::Int,
		ids::Vector, meta::Dict) where {T}
		check_timestamps(ts, size(dat, 1))
		N = size(dat, 2)
		ls = length(ids)
		if ls == 0
			ids = [Symbol(Char(x)) for x in 97:(97+N-1)]
		else
			ls == size(dat, 2) || throw(
				ArgumentError(
					"if ids are defined, they must one label for each sensor."),
			)
			ids = Symbol.(ids)
		end
		return new{N, T}(dat, ts, sr, NTuple{N, Symbol}(ids), meta)
	end
end;

function MultiForceData(
	df::DataFrame;
	sampling_rate::Integer,
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
		return x.sr
	elseif s === :n_samples
		return size(x.dat, 1)
	else
		return getfield(x, s)
	end
end


struct ForceProfiles{T <: AbstractFloat}
	dat::Matrix{T}
	sr::Int # sampling rate
	design::DataFrame
	baseline::Vector{T}
	zero_sample::Int

	function ForceProfiles(force::Matrix{T}, sr::Int, design::DataFrame,
		baseline::Vector{T}, zero_sample::Int) where {T <: AbstractFloat}
		lf = size(force, 1)
		lb = length(baseline)
		lf == lb || throw(
			ArgumentError(
				"Number of rows of force ($(lf)) must match the length of baseline ($(lb)).",
			),
		)
		return new{T}(force, sr, design, baseline, zero_sample::Int)
	end
end;


Base.propertynames(::ForceProfiles) = (:dat, :sampling_rate, :design, :baseline,
	:zero_sample, :n_samples, :n_profiles)
function Base.getproperty(x::ForceProfiles, s::Symbol)
	if s === :n_profiles
		return size(x.dat, 1)
	elseif s === :sampling_rate
		return x.sr
	elseif s === :n_samples
		return size(x.dat, 2)
	else
		return getfield(x, s)
	end
end


force(x::ForceData) = x.dat
force(x::ForceProfiles) = x.dat
force(x::MultiForceData, id::Int) = x.dat[:, id]
force(x::MultiForceData, id::String) = force(x, Symbol(id))

function ForceData(x::MultiForceData, id::Union{Integer, Symbol, String})
	return ForceData(force(x, id), x.ts, x.sr, x.meta)
end

function force(x::MultiForceData, id::Symbol)
	# returns force
	idx = findfirst(x.ids .== id)
	idx isa Integer || throw(ArgumentError(
		"Can not find force data labelled '$(id)'"))
	return force(x, idx)
end

function duration(n_samples::Int; sampling_rate::Int)
	# calculates n of samples into milliseconds duration
	return n_samples * (1000 / sampling_rate)
end


function Base.copy(fp::ForceProfiles)
	return ForceProfiles(copy(fp.dat), fp.sr, copy(fp.design),
		copy(fp.baseline), fp.zero_sample)
end


function Base.show(io::IO, mime::MIME"text/plain", x::ForceData)
	println(io, "ForceData")
	print(io, " $(x.n_samples) samples, sampling rate: $(x.sampling_rate)")
end;

function Base.show(io::IO, mime::MIME"text/plain", x::MultiForceData)
	println(io, "MultiForceData")
	print(io, "  $(x.n_samples) samples, sampling rate: $(x.sampling_rate)")
end;


function Base.show(io::IO, mime::MIME"text/plain", x::ForceProfiles)
	println(io, "ForceProfile")
	println(io, "  $(x.n_profiles) profiles")
	println(io, "  $(x.n_samples) samples, sampling rate: $(x.sampling_rate), zero sample: $(x.zero_sample)")
	print(io, "  Design: $(names(x.design))")
end

function load_force_profiles(filename::String)
	return convert(ForceProfiles, load_object(filename))
end;

function load_force_data(filename::String)
	return convert(ForceData, load_object(filename))
end;

function save_force(forces::Union{ForceProfiles, ForceData}, filename::String)
	if last(splitext(filename)) != ".jld2"
		filename = filename * ".jld2"
	end
	return jldsave(filename, true; forces)
end;

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
