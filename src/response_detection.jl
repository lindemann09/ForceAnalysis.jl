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

struct ForceResponse
	onset::Union{Missing, Int}
	offset::Union{Missing, Int}
	zero_sample::Int
	function ForceResponse(a::Union{Missing, Int}, b::Union{Missing, Int}, c::Int)
		if !ismissing(a) && !ismissing(b) && a > b
			throw(ArgumentError("Onset must be before offset."))
		end
		return new(a, b, c)
	end
end
ForceResponse(onset::Union{Missing, Int}, offset::Union{Missing, Int}) = ForceResponse(onset, offset, 0)

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

function response_detection(
	force_vector::AbstractVector{<:AbstractFloat},
	criterion::OnsetCriterion;
	zero_sample::Integer,
)::ForceResponse
	# returns ForceResponse relative to zero sample
	onset = _response_onset(force_vector, criterion)
	if onset >= 0
		offset = _response_offset(force_vector, criterion; onset)
		if offset < 0
			offset = missing
		end
		return ForceResponse(onset, offset, zero_sample)
	else
		return ForceResponse(missing, missing, 0)
	end
end


function response_detection(
	fe::ForceEpochs,
	criterion::OnsetCriterion,
)::Vector{ForceResponse}
	# returns ForceResponse relative to zero sample
	@unpack zero_sample = fe
	return [response_detection(f, criterion; zero_sample) for f in eachrow(fe.dat)]
end


function duration(
	rb::ForceResponse;
	sampling_rate::Real,
)
	# returns response latency in millisecond
	if ismissing(rb.onset) || ismissing(rb.offset)
		return missing
	end
	return _duration((rb.offset - rb.onset), sampling_rate)
end

function duration(
	vrb::AbstractVector{ForceResponse};
	sampling_rate::Real,
)
	return [_duration(rb, sampling_rate) for rb in vrb]
end

_duration(nsamples::Int64, sr::Float64) = nsamples * (1000.0 / sr)


function latency(
	rb::ForceResponse;
	sampling_rate::Real,
)
	# returns latency in millisecond
	if ismissing(rb.onset)
		return missing
	end
	return duration(rb.onset - rb.zero_sample; sampling_rate)
end

function latency(
	vrb::AbstractVector{ForceResponse};
	sampling_rate::Real,
)
	return [latency(rb; sampling_rate) for rb in vrb]
end


function peak_force(
	force_vector::AbstractVector{<:AbstractFloat},
	rb::ForceResponse,
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

function peak_force(
	fe::ForceEpochs,
	rb::AbstractVector{ForceResponse},
)
	# takes into account zero samples
	fe.n_epochs == length(rb) || throw(ArgumentError(
		"Number of epochs and ForceResponses don't match!"))
	return [peak_force(f, b) for (f, b) in zip(eachrow(fe.dat), rb)]
end


function impulse_size(
	force_vector::AbstractVector{<:AbstractFloat},
	rb::ForceResponse,
)   # TODO not yet tested
	resp = _extract_response(force_vector, rb)
	if !isnothing(resp)
		return sum(skipmissing(resp .- resp[1]))
	end
end

function impulse_size(
	fep::ForceEpochs,
	rb::AbstractVector{ForceResponse},
)
	fep.n_epochs == length(rb) || throw(ArgumentError(
		"Number of epochs and ForceResponse don't match!"))
	return [impulse_size(f, b) for (f, b) in zip(eachrow(fep.dat), rb)]
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
