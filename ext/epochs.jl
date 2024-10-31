
function Makie.plot!(ax::Axis, fe::BeForEpochs;
	rows::row_ids = nothing,
	ylims::UnitRange{Int} = -2000:2000,
	colors::VecOrColorant = RGBAf(0.2, 0.6, 0.2, 0.5),
	linewidth::Integer = 2,
	marker::Union{Nothing, AbstractVector{<:Integer}} = nothing,
	marker_color::Colorant = RGBAf(0.9, 0.4, 0.4, 0.8),
	marker_linewidth::Integer = 2,
	info_text::AbstractString = "",
	resp_criterion::Union{Nothing, OnsetCriterion} = nothing,
	mark_peak::Bool = true,
	kwargs...,
)
	if !isnothing(rows)
		fe = subset(fe, rows)
	end

	if !isnothing(resp_criterion)
		resp_mark = _response_marker(fe, resp_criterion; mark_peak)
		if isnothing(marker)
			marker = resp_mark
		else
			append!(marker, resp_mark)
		end
	end

	return _plot_force_matrix!(ax, force(fe); zero_sample = fe.zero_sample,
		ylims, colors, linewidth, marker, marker_color, marker_linewidth,
		info_text, kwargs...)
end

function Makie.plot!(fig::Figure, fe::BeForEpochs; kwargs...)
	return plot!(Axis(fig[1, 1]), fe; kwargs...)
end

function Makie.plot!(fig::Figure, epoch_mtx::Matrix; kwargs...)
	return _plot_force_matrix!(Axis(fig[1, 1]), epoch_mtx; kwargs...)
end

function ForceAnalysis.plot_good_bad!(ax::Axis, fe::BeForEpochs;
	rows::row_ids = nothing,
	ylims::UnitRange{Int} = -2000:2000,
	good_trials::Union{Nothing, AbstractVector{Bool}} = nothing,
	colors_good::Colorant = RGBAf(0.2, 0.6, 0.2, 0.3),
	color_bad::Colorant = RGBAf(0.9, 0.4, 0.4, 0.3),
	marker::Union{Nothing, AbstractVector{<:Integer}} = nothing,
	linewidth::Int = 2,
	info_text::AbstractString = "",
	kwargs...
)
	if isnothing(good_trials)
		# define sample range based on zero_sample
		colors = colors_good
	else
		colors = [x ? colors_good : color_bad for x in good_trials]
	end
	plot!(ax, fe; rows, ylims, colors, marker, linewidth, info_text,
				kwargs...)
	return ax
end

function ForceAnalysis.plot_good_bad!(fig::Figure, fe::BeForEpochs; kwargs...)
	return plot!(Axis(fig[1, 1]), fe; kwargs...)
end


function ForceAnalysis.plot_av_epoch!(ax::Axis, fe::BeForEpochs;
	condition::Symbol = :all,
	sd_err::Bool = true,
	colors::Union{<:ColorGradient, Vector{<:Colorant}, Nothing} = nothing,
	agg_fnc::Function = mean,
	linewidth::Real = 5,
	marker = Int64[],
	marker_color::Colorant = RGBAf(0.2, 0.2, 0.2, 0.9),
	highlight_ranges::Union{Nothing, Vector{UnitRange}} = nothing,
	highlight_color::Colorant = RGBAf(0.1, 0.9, 0.1, 0.3)
)
	# conditions is a variable with the conditions
	# has to have the same number of elemens as rows in froce

	xs = (1-fe.zero_sample):(fe.n_samples-fe.zero_sample)
	if length(marker) > 0
		vlines!(ax, marker, linewidth=1, color = marker_color)
	end

	if !isnothing(highlight_ranges)
		highlight_ranges!(ax, highlight_ranges, highlight_color)
	end

	if 	condition == :all
		agg_forces = aggregate(fe; agg_fnc = agg_fnc)
		cond = [:all]
		if sd_err
			sd_forces = aggregate(fe; agg_fnc = sderr)
		else
			sd_forces = nothing
		end
 	else
		agg_forces = aggregate(fe; condition, agg_fnc = agg_fnc)
		cond = agg_forces.design[:, condition]
		if sd_err
			sd_forces = aggregate(fe; condition, agg_fnc = sderr)
		else
			sd_forces = nothing
		end
	end

	if isnothing(colors)
		cols = cgrad(:roma, length(cond), categorical = true, rev=true)
	else
		cols = colors
	end
	for (i, c) in enumerate(sort(cond))
		val = vec(agg_forces.dat[cond.==c, :])
		lines!(xs, val; label = string(c), linewidth , color = cols[i])

		if !isnothing(sd_forces)
			err = vec(sd_forces.dat[cond.==c, :])
			band!(xs, val + err, val - err, color = (cols[i], 0.2))
		end
	end
	return ax
end

function ForceAnalysis.plot_av_epoch!(fig::Figure, fe::BeForEpochs; kwargs...)
	return plot_av_epoch!(Axis(fig[1, 1]), fe; kwargs...)
end
