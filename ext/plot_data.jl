const VecOrColorant = Union{Colorant, Base.AbstractVecOrTuple{Colorant}}
const row_ids = Union{Nothing, Integer, Base.AbstractVecOrTuple{Integer}}

function _response_marker(fe::BeForEpochs, oc::OnsetCriterion;
	mark_peak::Bool = true)
	rtn = Int[]
	responses = response_detection(fe, oc)
	for rb in responses
		if !ismissing(rb.onset)
			push!(rtn, rb.onset - rb.zero_sample)
		end
		if !ismissing(rb.offset)
			push!(rtn, rb.offset - rb.zero_sample)
		end
	end
	if mark_peak
		peaks = [peak_force(vec, resp) for (vec, resp) in zip(eachrow(fe.dat), responses)]
		for (r, p) in zip(responses, peaks)
			if p.sample_to_peak >= 0
				push!(rtn, p.sample_to_peak + r.onset - r.zero_sample)
			end
		end
	end
	return rtn
end

function _plot_force_matrix!(ax::Axis,
	force_mtx::Matrix{<:AbstractFloat};
	zero_sample::Integer = 0,
	ylims::UnitRange{Int} = -2000:2000,
	colors::VecOrColorant = RGBAf(0.2, 0.6, 0.2, 0.5),
	linewidth::Integer = 2,
	marker::Union{Nothing, AbstractVector{<:Integer}} = nothing,
	marker_color::Colorant = RGBAf(0.9, 0.4, 0.4, 0.8),
	marker_linewidth::Integer = 2,
	info_text::AbstractString = "",
	kwargs...,
)
	#epoch_lines
	xs = (1-zero_sample):(size(force_mtx, 2)-zero_sample)
	if colors isa Colorant
		colors = Iterators.cycle((colors,))
	end
	for (i, color) in zip(1:size(force_mtx, 1), colors)
		lines!(xs, force_mtx[i, :]; color, linewidth, kwargs...)
	end

	if !isnothing(marker)
		# write marker
		vlines!(ax, marker; linewidth = marker_linewidth, color = marker_color)
	end
	ylims!(ax, ylims.start, ylims.stop)
	if length(info_text) > 0
		text!(ax, 0, 1, text = info_text,
			align = (:left, :top), offset = (4, -2),
			space = :relative, fontsize = 18)
	end
	return ax
end

function ForceAnalysis.highlight_ranges!(ax::Axis, ranges::Base.AbstractVecOrTuple{UnitRange}, color::Any)
	start = [x.start for x in ranges]
	stop = [x.stop for x in ranges]
	vspan!(ax, start, stop; color)
end

ForceAnalysis.highlight_ranges!(ax::Axis, ranges::UnitRange, color::Any) = highlight_ranges!(ax, (ranges, ), color)

