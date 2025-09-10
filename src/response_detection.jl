"""
	OnsetCriterion{T <: AbstractFloat}

TODO
"""
struct OnsetCriterion{T <: AbstractFloat}
	min_increase::T  # min change between two samples
	min_increase_in_window::T
	window_size::Int
	function OnsetCriterion(min_increase::T, min_increase_in_window::T, window_size) where T <: AbstractFloat
		(min_increase >= 0 && min_increase_in_window >= 0
		 && window_size >= 0) || throw(ArgumentError(
			"No negative parameter allowed for OnsetCriterion"))
		return new{T}(min_increase, min_increase_in_window, window_size)
	end
end;

"""
	ForceResponse

TODO
"""
struct ForceResponse
	onset::Union{Missing, Int}
	offset::Union{Missing, Int}
	sampling_rate:: Real
	zero_sample::Int

	function ForceResponse(a::Union{Missing, Int}, b::Union{Missing, Int}, c::Float64, d::Int)
		(c > 0) || throw(ArgumentError("Sampling rate must be positive."))
		if !ismissing(a) && !ismissing(b) && a > b
			throw(ArgumentError("Onset must be before offset."))
		end
		return new(a, b, c, d)
	end
end
ForceResponse(onset::Union{Missing, Int}, offset::Union{Missing, Int},
			sampling_rate::Real) = ForceResponse(onset, offset, sampling_rate, 0)

function OnsetCriterion(;
	min_increase::Real,
	min_increase_in_window::Real,
	window_size::Integer,
)
	return OnsetCriterion(AbstractFloat(min_increase),
		AbstractFloat(min_increase_in_window), window_size)
end

function _response_onset(
	force_vector::AbstractVector{<:AbstractFloat},
	criterion::OnsetCriterion,
)::Int
	# returns the counter in the vector with the first onset

	# function for onset detection, finds first onset
	#
	# detects if  difference threshold (current-previous) is exceeded and force has
	# INCREASED after n trials for at least the increase_criterion
	#
	# returns sample of onset or -1

	force_diff = diff(force_vector)
	l = length(force_vector)
	if l == 0
		return -1
	end
	l - 1 == length(force_diff) || throw(ArgumentError(
		"length(force_difference) must be length(force_vector)-1 ."))

	for (c, d) in enumerate(force_diff)
		if !ismissing(d) && d > criterion.min_increase
			x = c + 1 # first difference refers to second force
			last = x + criterion.window_size
			if last > l
				last = l
			end
			init_force = force_vector[x]
			for f in skipmissing(force_vector[(x+1):last])
				if f > init_force + criterion.min_increase_in_window
					return x
				end
			end
		end
	end
	return -1
end

function _response_offset(
	force_vector::AbstractVector{<:AbstractFloat},
	criterion::OnsetCriterion;
	onset::Integer,
)::Int
	# returns the counter in the vector with the first offset

	if onset < 0
		return -1
	else
		# find approaching first crossing initial force
		l = length(force_vector)
		from = onset + criterion.window_size + 1
		initial_force = force_vector[onset]
		for x in range(from, l)
			f = force_vector[x]
			if !ismissing(f) && f <= initial_force
				return x
			end
		end
		# reverse data and find onset
		rev_onset = _response_onset(force_vector[l:-1:from], criterion)
		if rev_onset >= 0
			return l - rev_onset
		else
			return -1
		end
	end
end

"""
	response_detection(force_vector::AbstractVector{<:AbstractFloat}, sampling_rate:: Float64,
				criterion::OnsetCriterion; zero_sample::Integer)::ForceResponse
	response_detection(fe::BeForEpochs, sampling_rate:: Float64,
				criterion::OnsetCriterion)::Vector{ForceResponse}

Returns a `ForceResponse` or `Vector{ForceResponse}`
TODO
"""
function response_detection(
	force_vector::AbstractVector{<:AbstractFloat},
	sampling_rate::Float64,
	criterion::OnsetCriterion;
	zero_sample::Integer = 0
)::ForceResponse
	# returns ForceResponse relative to zero sample
	onset = _response_onset(force_vector, criterion)
	if onset >= 0
		offset = _response_offset(force_vector, criterion; onset)
		if offset < 0
			offset = missing
		end
		return ForceResponse(onset, offset, sampling_rate, zero_sample)
	else
		return ForceResponse(missing, missing, sampling_rate, zero_sample)
	end
end


function response_detection(
	fe::BeForEpochs,
	criterion::OnsetCriterion,
)::Vector{ForceResponse}
	# returns ForceResponse relative to zero sample
	zero_sample = fe.zero_sample
	sr = fe.sampling_rate
	return [response_detection(f, sr, criterion; zero_sample) for f in eachrow(fe.dat)]
end


"""
	duration
TODO
"""

function duration(rb::ForceResponse)
	# returns response latency in millisecond
	if ismissing(rb.onset) || ismissing(rb.offset)
		return missing
	end
	return _to_millisec(rb.offset - rb.onset, rb.sampling_rate)
end


"""
	latency
TODO
"""
function latency(rb::ForceResponse)
	# returns latency in millisecond
	if ismissing(rb.onset)
		return missing
	end
	return _to_millisec(rb.onset - rb.zero_sample, rb.sampling_rate)
end


"""
	peak_force(force_vector::AbstractVector{<:AbstractFloat}, rb::ForceResponse)

TODO
"""

function peak_force(
	force_vector::AbstractVector{<:AbstractFloat},
	rb::ForceResponse
)
	# returns (peak force, n sample to peak (reltive to onset))
	resp = _extract_response(force_vector, rb)
	if isnothing(resp)
		rtn = (missing, -1)
	else
		rtn = findmax(skipmissing(resp))
	end
	return (; peak = rtn[1], sample_to_peak = rtn[2])
end


"""
	impulse_size(force_vector::AbstractVector{<:AbstractFloat}, rb::ForceResponse)

TODO
"""
function impulse_size(
	force_vector::AbstractVector{<:AbstractFloat},
	rb::ForceResponse
)   # TODO not yet tested
	resp = _extract_response(force_vector, rb)
	if !isnothing(resp)
		return sum(skipmissing(resp .- resp[1]))
	end
end


# helper
function _extract_response(
	force_vector::AbstractVector{<:AbstractFloat},
	rb::ForceResponse,
)
	if ismissing(rb.offset) || ismissing(rb.offset)
		return nothing
	end
	return force_vector[(rb.onset):(rb.offset)]
end;

# samples to time
_to_millisec(nsamples::Int64, sr::Float64) = nsamples * (1000.0 / sr)
