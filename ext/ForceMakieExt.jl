module ForceMakieExt

using Makie
using Colors
using PlotUtils: ColorGradient

using ForceAnalysis

export 	plot!,
	plot_av_epoch!,
	plot_good_bad!,
	highlight_ranges!

include("plot_data.jl")
include("epochs.jl")

end # module
